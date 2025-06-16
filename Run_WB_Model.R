# ---------------------------------------------------------------------
# This script includes code for running and calibrating the water balance model. The model can be run with
# user-provided model paramters, or the parameters can be calibrated based on measured ET values from 
# the OpenET dataset. 
#
# EDITS IN PROGRESS
# figures for comparison with openET observations
# add code for future projections
# optimization with AET does not appear to work very well, at least at the monthly time step (NSE < 0)
# add functionality to run the model many times (like Joseph's wrapper script) - not sure what this should look like
# clean up, generally make user friendly
# add info on how to run to the user manual
# add nice visualizations (including some from Janelle/Connor's code?)
# ---------------------------------------------------------------------



#######################################################################
#######################################################################
### INTRODUCTION SECTION ###


#######################################################################
### Load libraries ###
library(sf); library(raster); library(ggplot2); library(dplyr); library(xts); library(geosphere)
library(lubridate); library(hydroGOF); library(stringr); library(terra); library(glue); library(tidyverse)
library(climateR); library(EGRET); library(daymetr); library(here); library(ggrepel); library(gridExtra); library(Kendall)
library(httr); library(jsonlite); library(sf); library(grid); library(GA); library(GGally); library(data.table)

### Source in function files ###
path <- here() 
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
### Set user-defined variables ###
PETMethod = "Oudin" 
optimization = TRUE 
delayStart = TRUE 
incompleteMonths = FALSE 
GridMET = TRUE
fillLeapDays = TRUE 
future_analysis = TRUE
runFutureWB = TRUE  # TRUE to re-run entire water balance model for future; FALSE to use pre-existing water balance projections from a Mike Tercek spreadsheet
userSetJTemp = FALSE 
make_plots = TRUE 
provide_coords = FALSE # if true, user provides lat/lon coords. if false, lat/long coords are pulled from centroid of watershed with given gage id
FolderName = "optim" 

### Define watershed ###
# centroid of watershed
# must be west of approximately -90 longitude to have openET data
SiteID = "Redwood Creek"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '11460151'                     #define stream gage location
if(provide_coords){
  lat = 37.9 
  lon = -122.59 
}

### Define time period for historical analysis ###
# for GridMET and stream gage; Daymet period starts one year after this period 
# openET data (monthly) is only available 2000 - present
startY = 2000; startM = 01; startD = 01 
endY = 2024; endM = 12; endD = 31


### Model names ###
# provide list of model names to generate plots of future streamflow with those models highlighted
individual_models = c('BNU-ESM.rcp45')


#######################################################################
### Set model variables ###

# default Water Balance variables
gw_add=0; vfm = 0.7555; jtemp = 1.982841; jrange = 3 ;hock = 4; hockros = 4; 
dro = 0; mondro = 0; aspect = 180; slope= 0; shade.coeff= 1; SWC.Max = 200

# Water Balance variables not to be optimized
Soil.Init = SWC.Max; Snowpack.Init = 0; T.Base = 0  

# water balance optimization lower and upper limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400)


#######################################################################
### Set other variables ###
# Optional scaling factors for GridMET time series: if no scaling, set slopes to 1 and bias to 0
tmmx_slope = 1; tmmx_bias = 0
tmmn_slope = 1; tmmn_bias = 0
p_slope = 1; p_bias = 0



#######################################################################
#######################################################################
### GET DATA ###

#######################################################################
### Establish variables, file paths, and names ###

# Set path variables
if(!dir.exists(here('Data', SiteID_FileName))) {dir.create(here('Data', SiteID_FileName))}; dataPath <- here('Data', SiteID_FileName)
if(!dir.exists(here('Output', SiteID_FileName))) {dir.create(here('Output', SiteID_FileName))}
if(!dir.exists(here('Output', SiteID_FileName, 'WaterBalance'))) {dir.create(here('Output', SiteID_FileName, 'WaterBalance'))}
if(!dir.exists(here('Output', SiteID_FileName, 'WaterBalance', FolderName))) {dir.create(here('Output', SiteID_FileName,'WaterBalance', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'WaterBalance', FolderName)

# Pull watershed shapefile from StreamStats database
# figure out format / how to assign variables
if(!provide_coords){
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon
}

# Get j_temp
if(!userSetJTemp){
  j.raster = raster(here('Data', "merged_jennings.tif"))
  jtemp = get_jtemp(lat = lat, lon= lon, j.raster = j.raster)}

#define lower and upper boundaries for jtemp
WB_lower = c(WB_lower, jtemp = jtemp-0.5) 
WB_upper = c(WB_upper, jtemp= jtemp+0.5)

#get elevation from Daymet data. This happens regardless or whether you use Daymet or GridMET climate data
elev = get_elev_daymet(lat, lon, startY, endY, SiteID_FileName)

# create start and end date objects of data collection. Daymet will start one year after the year listed here
startDate<- ymd(paste(startY, startM, startD)); endDate<-  ymd(paste(endY, endM, endD))



#######################################################################
### Scrape and clean openET data ###

if(optimization){
  MonthlyET <- get_et_point(startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', dataPath, api_key)
  if(startY < 2016){startY_daily <- 2016
  } else{startY_daily <- startY}
  DailyET <- get_et_point(startY_daily, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', dataPath, api_key)
}


#######################################################################
### Scrape and clean meteorological data (GridMET or Daymet) ### 

if(GridMET) {
  DailyClimData <- get_gridmet_data(SiteID_FileName, startY, endY, lat, lon, dataPath,
                                    tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope)
} else { 
  DailyClimData <- get_daymet_data(SiteID_FileName, startY, endY, lat, lon, dataPath)
}



#######################################################################
#######################################################################
### MODEL RUNNING ### 


results <- data.frame(SiteID = SiteID, start = startDate, end = endDate, PETMethod = PETMethod, optimization = optimization,
                      GridMET = GridMET, lon = lon, lat = lat,
                      startY = startY, startM = startM, startD = startD, endY = endY, endM = endM, endD = endD,
                      incompleteMonths = incompleteMonths)




#######################################################################
#######################################################################
### OPTIMIZATION ###

### First optimization: optimize water balance variables according to the NSE of AET/openET over historical period ###
# need new function
if(optimization){
  parms<- c(gw_add = gw_add, vfm = vfm, jrange = jrange, hock =  hock, hockros = hockros,dro = dro, mondro = mondro,
            aspect = aspect,slope= slope, shade.coeff= shade.coeff, SWC.Max = SWC.Max,
            jtemp = jtemp)
  
  #run the optimization routine
  strtTimeM <-Sys.time()
  set.seed(123) 
  WBcoeffs <- tibble()
  
  # try GA
  # returned value is labeled as l_par - not sure why???
  optMonth_init <- ga(type = "real-valued", fitness = function(x) 
    WB_Optim_AET(c(gw_add=x[1], vfm=x[2], jrange=x[3], hock=x[4], hockros=x[5], dro=x[6], mondro=x[7], aspect=x[8], slope=x[9], shade.coeff=x[10], SWC.Max=x[11], jtemp=x[12]), 
             Soil.Init = Soil.Init, Snowpack.Init = Snowpack.Init, T.Base = T.Base, PETMethod= PETMethod, 
             DailyClimData = DailyClimData, lat=lat,lon=lon, meas_aet_mon = DailyET), 
    lower=WB_lower, upper=WB_upper) #, maxiter=300)
  elpTimeM <- Sys.time() - strtTimeM
  
  # Define the water balance variables from the best run
  optValuesM <- data.frame(nseM = optMonth_init@fitnessValue, optMonth_init@solution)
  gw_add=optValuesM$gw_add; vfm=optValuesM$vfm; jrange=optValuesM$jrange; hock=optValuesM$hock
  hockros=optValuesM$hockros; dro=optValuesM$dro; mondro=optValuesM$mondro; aspect=optValuesM$aspect
  slope=optValuesM$slope; shade.coeff=optValuesM$shade.coeff; SWC.Max=optValuesM$SWC.Max; jtemp=optValuesM$jtemp
  
  # store and save results
  results = data.frame(results, optValuesM, elpTimeM = elpTimeM)
  saveRDS(WBcoeffs, file = paste0(outLocationPath, "/WBcoeffs.rds"))
  
  
  ### Rerun model with optimal variables ###
  DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect,
               slope, shade.coeff, jtemp, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
  
  
  ### PLOTS
  if(make_plots){
    if (dev.cur() != 1) dev.off()
    # Parallel coordinates plot of WB coeffs
    jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_ParallelCoords.jpg"), width=1600, height=800, res=100)
    print(ggparcoord(data=WBcoeffs, columns=1:12, groupColumn=13, scale="uniminmax") + 
            scale_color_gradient(low = "black", high = "gray90") + nps_theme()) # for red: "darkred" and "#fee5d9"
    dev.off()
    
    # Scatterplots of WB coeffs
    WBcoeffs_long <- reshape(WBcoeffs, varying = names(WBcoeffs)[1:12], v.names = "WBcoeffs", timevar = "Variable", 
                             times = names(WBcoeffs)[1:12], direction = "long")
    jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_Scatter.jpg"), width=800, height=500)
    print(ggplot(WBcoeffs_long, aes(x = WBcoeffs, y = nseM)) + geom_point() +
            facet_wrap(~ Variable, scales = 'free') + nps_theme() +
            labs(title = 'Water Balance Coefficients', x='', y = 'Monthly NSE'))
    dev.off()
    
    # Range plot of optimal WB coeffs
    wb_optim <- data.frame(var=c('Groundwater Addition', 'Volume Forcing Multiplier', 'Jennings Temperature Range','Hock','Hock Rain on Snow','Direct Runoff','Mondro','Aspect','Slope','Shade Coefficient','Max Soil Water Content','Jennings Temperature'),
                           value=c(gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope, shade.coeff, SWC.Max, jtemp),
                           lower=WB_lower, upper=WB_upper)
    jpeg(file=paste0(outLocationPath, "/", "WB_Optim_Coefficients.jpg"), width=600, height=400); par(mfrow = c(3, 4))
    for (i in 1:12){
      len <- ceiling(wb_optim$upper[i])-floor(wb_optim$lower[i])
      plot(c(wb_optim$lower[i], wb_optim$upper[i]), c(0, 0), type = "n", xlab = "", ylab = "",
           main = wb_optim$var[i], xlim = range(c(wb_optim$lower[i]-(len/5), wb_optim$upper[i]+(len/5))), ylim = c(-1, 1), xaxt = "n", yaxt = "n")
      segments(wb_optim$lower[i], 0, wb_optim$upper[i], 0, col = "black", lwd = 2)
      points(wb_optim$value[i], 0, col = "red", pch = 19, cex = 1.5)
      text(wb_optim$value[i], 0.3, sprintf("%.2f", wb_optim$value[i]), col='red')
      axis(1, at = seq(floor(wb_optim$lower[i]), ceiling(wb_optim$upper[i]), by=(len/5)))
    }
    dev.off()
  }
}





#######################################################################
### Run model without optimization ###
if(!optimization){
  # run model
  DailyWB <- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon)
  #MeasMod<- MeasModWB(DailyDrain, meas_flow_mon)
  #nseM = NSE(sim = MeasMod$Mod, obs = MeasMod$Meas)
  #results = data.frame(results, nseM =nseM, gw_add = gw_add, vfm =vfm,jrange =jrange,hock = hock, hockros =hockros,
  #                     dro =dro, mondro = mondro, aspect = aspect, slope = slope,
  #                     shade.coeff = shade.coeff, jtemp =jtemp, elpTimeM = NA)
  
  # store results
}


### save results as RDS files ###
saveRDS(results, file = paste0(outLocationPath, "/results.rds"))



#######################################################################
#######################################################################
##### MODEL PERFORMANCE ON HISTORICAL ET DATA #####
# EDIT

# Daily aggregation with measured and modeled ET
hist_et_daily <- merge(xts(with(DailyWB, cbind(AET)), order.by = as.Date(DailyWB$date)), meas_flow_daily_xts)
hist_et_daily <- hist_et_daily[complete.cases(hist_et_daily),]
colnames(hist_et_daily) <- c("AET", "Mod", "Meas")

# Monthly aggregation with measured and modeled ET
MeasMod <- MeasModWB(DailyDrain = DailyDrain, meas_et_mon = meas_et_mon, cutoffYear = cutoffYear)

# Annual aggregation with measured and modeled streamet
hist_et_ann <- as.data.frame(matrix(NA, nrow = nrow(apply.yearly(hist_et_daily[,"adj_runoff"], sum)),
                                      ncol = ncol(hist_et_daily), dimnames = list(c(), colnames(hist_et_daily))))
for(i in 1:ncol(hist_et_daily)){
  hist_et_ann[,i] <- apply.yearly(hist_et_daily[,i], sum)
}


