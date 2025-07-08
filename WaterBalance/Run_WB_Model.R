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
### Load libraries and source in function files ###
library('here'); lib_install <- FALSE
path <- here(); setwd(here('Code'))
sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
### Set user-defined variables ###
PETMethod = "Oudin" 
optimization = TRUE 
delayStart = TRUE 
incompleteMonths = FALSE 
GridMET = TRUE
fillLeapDays = TRUE 
future_analysis = TRUE
runFutureWB = TRUE
calcFutureWB = TRUE
userSetJTemp = FALSE 
make_plots = TRUE 
provide_coords = FALSE
point_location = FALSE
percent_skill_cutoff = 0.1 
FolderName = 'optim' 
filename = ''

### Define watershed ###
# centroid of watershed
# must be west of approximately -90 longitude to have OpenET data
SiteID = "Redwood Creek"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '11460151'        # USGS stream gage 
if(provide_coords){
  aoi <- ""
  lat = 37.9 
  lon = -122.59 
}

### Define time period for historical analysis ###
# for GridMET and stream gage; Daymet period starts one year after this period 
# openET data (monthly) is only available 2000 - present
# apparently right now monthly data is only available from 2010-2019? not sure why, check if this changes
startY = 2000; startM = 01; startD = 01 
endY = 2024; endM = 12; endD = 31


### Model names ###
# provide list of model names to generate plots of future streamflow with those models highlighted
individual_models = c('BNU-ESM.rcp45')


#######################################################################
### Set model variables ###

# Default water balance variables
gw_add=0; vfm = 0.7555; jtemp = 1.982841; jrange = 3; hock = 4; hockros = 4; 
dro = 0; mondro = 0; aspect = 180; slope= 0; shade.coeff= 1; SWC.Max = 200
aet_slope = 1; aet_bias = 0
Soil.Init = SWC.Max; Snowpack.Init = 0; T.Base = 0  

# Water balance optimization lower and upper limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10, aet_slope=-10, aet_bias=-5)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400, aet_slope=10, aet_bias=5)


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
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon; elev <- coords$elev; aoi <- coords$geometry
}
region <- get_region(lat,lon)

# Get j_temp
if(!userSetJTemp){
  j.raster = raster(here('Data', "merged_jennings.tif"))
  jtemp = get_jtemp(lat = lat, lon= lon, j.raster = j.raster)}

#define lower and upper boundaries for jtemp
WB_lower = c(WB_lower, jtemp = jtemp-0.5) 
WB_upper = c(WB_upper, jtemp= jtemp+0.5)

# create start and end date objects of data collection. Daymet will start one year after the year listed here
startDate<- ymd(paste(startY, startM, startD)); endDate<-  ymd(paste(endY, endM, endD))


#######################################################################
### Scrape and clean meteorological data (GridMET or Daymet) ### 

if(GridMET) {
  if(point_location){
    DailyClimData <- get_gridmet_point(SiteID_FileName, startY, endY, lat, lon, dataPath,
                                       tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope)
  } else {
    DailyClimData <- get_gridmet_area(SiteID_FileName, startY, endY, aoi, dataPath,
                                      tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope)
  }
} else { 
  if(point_location){
    DailyClimData <- get_daymet_point(SiteID_FileName, startY, endY, lat, lon, dataPath)
  } else{
    DailyClimData <- get_daymet_area(SiteID_FileName, startY, endY, aoi, dataPath)
  }
}



#######################################################################
### Scrape and clean openET data ###

MonthlyET <- get_et_point(startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ET', dataPath)
DailyET <- get_et_point(startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ET', dataPath)

MonthlyETo <- get_et_point(startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ETo', dataPath)
DailyETo <- get_et_point(startY_daily, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ETo', dataPath)



#######################################################################
#######################################################################
### OPTIMIZATION ###

### First optimization: optimize water balance variables according to the NSE of AET/openET over historical period ###
# need new function
if(optimization){
  # Create file for results
  results <- data.frame(SiteID = SiteID, start = startDate, end = endDate, PETMethod = PETMethod, optimization = optimization,
                        GridMET = GridMET, lon = lon, lat = lat,
                        startY = startY, startM = startM, startD = startD, endY = endY, endM = endM, endD = endD,
                        incompleteMonths = incompleteMonths)
  
  parms<- c(gw_add = gw_add, vfm = vfm, jrange = jrange, hock =  hock, hockros = hockros,dro = dro, mondro = mondro,
            aspect = aspect,slope= slope, shade.coeff= shade.coeff, SWC.Max = SWC.Max, aet_slope=aet_slope, aet_bias=aet_bias, jtemp = jtemp)
  
  # Run optimization routine
  strtTimeM <-Sys.time()
  #set.seed(123) 
  WBcoeffs <- tibble()
  
  # try GA - monthly time scale to smooth out errors
  # returned value is labeled as l_par - not sure why???
  optMonth_init <- ga(type = "real-valued", fitness = function(x)
    WB_Optim_AET(c(gw_add=x[1], vfm=x[2], jrange=x[3], hock=x[4], hockros=x[5], dro=x[6], mondro=x[7], aspect=x[8], slope=x[9], shade.coeff=x[10], SWC.Max=x[11], aet_slope=x[12], aet_bias=x[13], jtemp=x[14]),
             Soil.Init = Soil.Init, Snowpack.Init = Snowpack.Init, T.Base = T.Base, PETMethod= PETMethod,
             DailyClimData = DailyClimData, lat=lat, lon=lon, meas_aet = MonthlyET, interval='monthly'),
    lower=WB_lower, upper=WB_upper, popSize=200, maxiter=100, pmutation = 0.2, pcrossover = 0.8, elitism = 10, run = 50, monitor = TRUE )

  # see if DEoptim produces better optimization results
  #optMonth_init <- DEoptim(fn=function(x) -WB_Optim_AET(c(gw_add=x[1], vfm=x[2], jrange=x[3], hock=x[4], hockros=x[5], dro=x[6], mondro=x[7], aspect=x[8], slope=x[9], shade.coeff=x[10], SWC.Max=x[11], aet_slope=x[12], aet_bias=x[13], jtemp=x[14]), 
  #                                                                Soil.Init = Soil.Init, Snowpack.Init = Snowpack.Init, T.Base = T.Base, PETMethod= PETMethod, 
  #                                                                DailyClimData = DailyClimData, lat=lat, lon=lon, meas_aet = MonthlyET, interval='monthly'), 
  #                         lower=WB_lower, upper=WB_upper,  control = DEoptim.control(VTR=0.9, NP = 200, itermax = 100, trace = TRUE))
  elpTimeM <- Sys.time() - strtTimeM
  
  # Define the water balance variables from the best run
  optValuesM <- data.frame(nseM = optMonth_init@fitnessValue, optMonth_init@solution)
  #optValuesM <- as.data.frame(t(optMonth_init$optim$bestmem)); colnames(optValuesM) <- c('gw_add','vfm','jrange','hock','hockros','dro','mondro','aspect','slope','shade.coeff','SWC.Max','aet_slope','aet_bias','jtemp')
  #optValuesM <- optValuesM %>% mutate(nseM = optMonth_init$optim$bestval)
  gw_add=optValuesM$gw_add; vfm=optValuesM$vfm; jrange=optValuesM$jrange; hock=optValuesM$hock
  hockros=optValuesM$hockros; dro=optValuesM$dro; mondro=optValuesM$mondro; aspect=optValuesM$aspect
  slope=optValuesM$slope; shade.coeff=optValuesM$shade.coeff; SWC.Max=optValuesM$SWC.Max
  aet_slope=optValuesM$aet_slope; aet_bias=optValuesM$aet_bias; jtemp=optValuesM$jtemp
  
  # store and save results
  results = data.frame(results, optValuesM, elpTimeM = elpTimeM)
  saveRDS(WBcoeffs, file = paste0(outLocationPath, "/WBcoeffs.rds"))
  
  # print time
  print(paste("Time elapsed:", elpTimeM))
  
  ### PLOTS
  if(make_plots){
    if (dev.cur() != 1) dev.off()
    # Parallel coordinates plot of WB coeffs
    jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_ParallelCoords.jpg"), width=1600, height=800, res=100)
    print(ggparcoord(data=WBcoeffs, columns=1:14, groupColumn=15, scale="uniminmax") + 
            scale_color_gradient(low = "black", high = "gray90") + nps_theme()) # for red: "darkred" and "#fee5d9"
    dev.off()
    
    # Scatterplots of WB coeffs
    WBcoeffs_long <- reshape(WBcoeffs, varying = names(WBcoeffs)[1:14], v.names = "WBcoeffs", timevar = "Variable", 
                             times = names(WBcoeffs)[1:14], direction = "long")
    jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_Scatter.jpg"), width=800, height=500)
    print(ggplot(WBcoeffs_long, aes(x = WBcoeffs, y = nse)) + geom_point() +
            facet_wrap(~ Variable, scales = 'free') + nps_theme() +
            labs(title = 'Water Balance Coefficients', x='', y = 'Monthly NSE'))
    dev.off()
    
    # Range plot of optimal WB coeffs
    wb_optim <- data.frame(var=c('Groundwater Addition', 'Volume Forcing Multiplier', 'Jennings Temperature Range','Hock','Hock Rain on Snow','Direct Runoff','Mondro','Aspect','Slope','Shade Coefficient','Max Soil Water Content','AET Slope','AET Bias','Jennings Temperature'),
                           value=c(gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope, shade.coeff, SWC.Max, aet_slope, aet_bias, jtemp),
                           lower=WB_lower, upper=WB_upper)
    jpeg(file=paste0(outLocationPath, "/", "WB_Optim_Coefficients.jpg"), width=800, height=600); par(mfrow = c(3, 5))
    for (i in 1:14){
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
#######################################################################
### Run model  ###
DailyWB <- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, aet_slope, aet_bias, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon)



#######################################################################
#######################################################################
##### MODEL PERFORMANCE ON HISTORICAL ET DATA #####
# EDIT

# Daily, monthly, annual aggregation with measured and modeled ET
hist_et_daily <- merge(xts(with(DailyWB, cbind(AET)), order.by = as.Date(DailyWB$date)), DailyET); hist_et_daily <- hist_et_daily[complete.cases(hist_et_daily),]
colnames(hist_et_daily) <- c("Mod", "Meas")
hist_et_monthly <- apply.monthly(xts(with(DailyWB, cbind(AET)), order.by = as.Date(DailyWB$date)), function(x) {colSums(x, na.rm = TRUE)})
index(hist_et_monthly) <- as.Date(format(index(hist_et_monthly), "%Y-%m-01"))
hist_et_monthly <- merge(hist_et_monthly, MonthlyET); hist_et_monthly <- hist_et_monthly[complete.cases(hist_et_monthly),]; colnames(hist_et_monthly) <- c("Mod", "Meas")
hist_et_ann <- apply.yearly(hist_et_monthly, function(x) {colSums(x, na.rm = TRUE)})


#######################################################################
### Create summary plots ###

# scatterplot of Historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
# there are two trend lines in the scatter plot because the intercept is set to 0 in one and allowed to vary in the other
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_Scatter.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  meas <- coredata(hist_et_daily$Meas); mod <- coredata(hist_et_daily$Mod)
  plot(mod, meas, main =paste("Daily Average AET\nNSE:",round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # monthly
  mod <- coredata(hist_et_monthly$Mod); meas <- coredata(hist_et_monthly$Meas)
  plot(mod, meas, main =paste("Monthly Total AET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # annual
  mod <- coredata(hist_et_ann$Mod); meas <- coredata(hist_et_ann$Meas)
  plot(mod,meas, main =paste("Annual Total AET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  dev.off()
}

# time series plot of historical Measured vs Modeled AET for daily, monthly, annual aggregation
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_TimeSeries.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  plot(hist_et_daily[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Daily AET (mm)", main = "Daily", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # monthly
  plot(hist_et_monthly[,c("Mod", "Meas")], type = "l", lwd = 2, xlab = "Date", ylab = "Monthly Sum AET (mm)", main = "Monthly", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # annual
  plot(hist_et_ann[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Annual Sum AET (mm)", main = "Annual", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}





