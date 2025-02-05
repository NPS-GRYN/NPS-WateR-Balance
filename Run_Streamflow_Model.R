# ---------------------------------------------------------------------
# This script contains code for running and calibrating the streamflow model (which is comprised of
# the water balance model and the IHACRES rainfall-streamflow runoff model). The model can be run with
# user-provided model paramters, or the parameters can be calibrated based on measured streamflow values.
# The model is designed to be run for a USGS streamflow gage, to allow for calibration, but users can also
# provide their own coordinates if they simply want to run the model.
# See the User Manual for a more detailed description of the model.
#
# EDITS IN PROGRESS:
# add functionality to run the model many times (like Joseph's wrapper script) - not sure what this should look like
# clean up, generally make user friendly
# add info on how to run to the user manual
# ---------------------------------------------------------------------


#######################################################################
#######################################################################
### INTRODUCTION SECTION ###


#######################################################################
### Load libraries ###
library(sf); library(raster); library(ggplot2); library(dplyr); library(xts); library(geosphere)
library(lubridate); library(hydroGOF); library(stringr); library(terra); library(glue); library(tidyverse)
library(climateR); library(EGRET); library(daymetr); library(here); library(ggrepel); library(gridExtra); 
library(httr); library(jsonlite); library(sf); library(grid); library(GA); library(GGally)

### Source in function files ###
path <- here() 
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
### Set user-defined variables ###
PETMethod = "Oudin" 
optimization = TRUE 
delayStart = TRUE 
NonZeroDrainInitCoeff = FALSE
incompleteMonths = FALSE 
GridMET = TRUE
fillLeapDays = TRUE 
future_analysis = TRUE
runFutureWB = TRUE  # TRUE to re-run entire water balance model for future; FALSE to use pre-existing water balance projections from a Mike Tercek spreadsheet
userSetJTemp = FALSE 
make_plots = TRUE 
provide_coords = FALSE # if true, user provides lat/lon coords. if false, lat/lon coords are pulled from centroid of watershed with given gage id
flow_components = 3  # change the number of components that characterize the flow. can be 2 or 3. 2: flow has quick and slow components; 3: flow has quick, slow, and very slow components.
FolderName = "optim" 

### Define watershed ###
# centroid of watershed
SiteID = "Cataloochee"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '03460000'                  #define stream gage location (RWC: "11460151")
if(provide_coords) lat = 37.9; lon = -122.59 

### Define time period for historical analysis ###
# for GridMET and stream gage; Daymet period starts one year after this period 
startY = 1979; startM = 01; startD = 01 
endY = 2024; endM = 12; endD = 31


### Model names ###
# provide list of model names to generate plots of future streamflow with those models highlighted
individual_models = c('BNU-ESM.rcp45')


#######################################################################
### Set model variables ###
# Initial IHACRES flow coefficients
q0<-0; s0<-0; v0<-0
if(flow_components==3){
  qa<-0.62; qb<-0.22; sa<-0.58; sb<-0.06; va<-0.974; vb<-calc_vb(qa,qb,sa,sb,va)
  #qa<- 1; qb<-0; sa<-1; sb<-0; va<-1; vb<-calc_vb(qa,qb,sa,sb,va)
}else if(flow_components==2){
  qa<-0; qb<-0; sa<-0.58; sb<-0.06; va<-0.974; vb<-calc_vb(qa,qb,sa,sb,va) 
} else{print('invalid number of flow components')}

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
if(!dir.exists(here('Output', SiteID_FileName, 'IHACRES'))) {dir.create(here('Output', SiteID_FileName, 'IHACRES'))}
if(!dir.exists(here('Output', SiteID_FileName, 'IHACRES', FolderName))) {dir.create(here('Output', SiteID_FileName, 'IHACRES', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'IHACRES', FolderName)

# Pull watershed shapefile from StreamStats database
# figure out format / how to assign variables
if(!provide_coords){
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon
}

# Define variables that do not need to be defined outside of the function
# is this supposed to be 11 or 1???
if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 

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
### Scrape and clean USGS stream gage data ###

gage_data <- get_gage_data(GageSiteID, incompleteMonths, dataPath)
meas_flow_daily_xts <- gage_data$meas_flow_daily_xts; meas_flow_mon <- gage_data$meas_flow_mon


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

### Get initial flow conditions ###
if(NonZeroDrainInitCoeff){
  InitCond <- get_Init_Drain_Coef(DailyClimData, gw_add, vfm , jrange ,hock, hockros,
                                  dro, mondro , aspect, slope, shade.coeff, jtemp, SWC.Max, 
                                  Soil.Init, Snowpack.Init, T.Base, PETMethod, 
                                  q0, s0, v0, qa, qb, sa, sb, va, vb, lat, lon, cutoffYear)
  
  q0 = InitCond[["Quick"]]; s0= InitCond[["Slow"]]; v0 = InitCond[["Very_Slow"]]
}

results <- data.frame(SiteID = SiteID, start = startDate, end = endDate, PETMethod = PETMethod, optimization = optimization,
                      GridMET = GridMET, lon = lon, lat = lat,
                      startY = startY, startM = startM, startD = startD, endY = endY, endM = endM, endD = endD,
                      cutoffYear = cutoffYear, NonZeroDrainInitCoeff = NonZeroDrainInitCoeff, incompleteMonths = incompleteMonths)




#######################################################################
#######################################################################
### OPTIMIZATION ###
# look into optim(): increment jrange in whole numbers

### First optimization: optimize water balance variables according to the NSE of monthly summed streamflow over historical period ###
if(optimization){
  parms<- c(gw_add = gw_add, vfm = vfm, jrange = jrange, hock =  hock, hockros = hockros,dro = dro, mondro = mondro,
            aspect = aspect,slope= slope, shade.coeff= shade.coeff, SWC.Max = SWC.Max, jtemp = jtemp)
  
  #run the optimization routine
  strtTimeM <-Sys.time()
  set.seed(123) #this ensures reproducibility each time
  WBcoeffs <- tibble()
  
  # Use genetic algorithm (GA) for optimization
  optMonth_init <- ga(type = "real-valued", fitness = function(x) 
    WB_Optim(c(gw_add=x[1], vfm=x[2], jrange=x[3], hock=x[4], hockros=x[5], dro=x[6], mondro=x[7], aspect=x[8], slope=x[9], shade.coeff=x[10], SWC.Max=x[11], jtemp=x[12]), 
            meas_flow_daily_xts = meas_flow_daily_xts, cutoffYear = cutoffYear, q0=q0, s0=s0, v0=v0,qa=qa, qb=qb, sa=sa, sb=sb,va=va, vb=vb, Soil.Init = Soil.Init, 
            Snowpack.Init = Snowpack.Init, T.Base = T.Base, PETMethod= PETMethod, DailyClimData = DailyClimData, lat=lat,lon=lon, meas_flow_mon = meas_flow_mon), 
    lower=WB_lower, upper=WB_upper)
  elpTimeM <- Sys.time() - strtTimeM
  
  # Define the water balance variables from the best run
  optValuesM <- data.frame(nseM = optMonth_init@fitnessValue, optMonth_init@solution)
  gw_add=optValuesM$gw_add; vfm=optValuesM$vfm; jrange=optValuesM$jrange; hock=optValuesM$hock
  hockros=optValuesM$hockros; dro=optValuesM$dro; mondro=optValuesM$mondro; aspect=optValuesM$aspect
  slope=optValuesM$slope; shade.coeff=optValuesM$shade.coeff; SWC.Max=optValuesM$SWC.Max; jtemp=optValuesM$jtemp
  
  # store and save results
  results = data.frame(results, optValuesM, elpTimeM = elpTimeM)
  saveRDS(WBcoeffs, file = paste0(outLocationPath, "/WBcoeffs.rds"))
  
  
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
### Second optimization: optimize IHACRES A and B coefficients according to the NSE of daily streamflow over historical period ###
# check how it works w 2 flow components - used to be bad
if(optimization){
  #Re run the Water Balance model as an input for the second optimization with the optimized water balance variables
  DailyWB<- WB(DailyClimData, gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon)
  
  # Define variables for optimization
  IHACRES_lower <- c(l_qa=0, l_qb=0, l_sa=0, l_sb=0, l_va=0)
  IHACRES_upper <- c(u_qa=1, u_qb=1, u_sa=1, u_sb=1, u_va=1)
  strtTimeD <-Sys.time()
  set.seed(123)
  
  
  ### use GA (genetic algorithm) to explore parameter space and estimate IHACRES flow coefficients ###
  IHcoeffs <- tibble()
  optDaily_init <- ga(type = "real-valued", fitness = function(x) IHACRESFlow(c(qa=x[1], qb=x[2], sa=x[3], sb=x[4], va=x[5]), q0=q0, s0=s0, v0=v0, DailyWB=DailyWB, meas_flow_daily_xts= meas_flow_daily_xts, cutoffYear = cutoffYear), lower =IHACRES_lower, upper = IHACRES_upper)
  if(flow_components==3){
    qa <- optDaily_init@solution[1]; qb <- optDaily_init@solution[2]; sa <- optDaily_init@solution[3]; sb <- optDaily_init@solution[4]; va <- optDaily_init@solution[5]; vb<-calc_vb(qa, qb, sa, sa, va)
  } else if(flow_components==2){
    optSol <- IHcoeffs[which.max(IHcoeffs$nseD), ]
    qa<-optSol$qa; qb<-optSol$qb; sa<-optSol$sa; sb<-optSol$sb; va<-optSol$va; vb<-optSol$vb
  }
  saveRDS(IHcoeffs, file = paste0(outLocationPath, "/IHcoeffs_initial.rds"))
  
  # plot 
  if(make_plots){
    # Parallel coordinates plot of IH coeffs
    jpeg(file=paste0(outLocationPath, "/", "IHACRES_Coeffs_ParallelCoords.jpg"), width=600, height=400)
    print(ggparcoord(data=IHcoeffs, columns=1:6, groupColumn=7, scale="uniminmax") + 
      scale_color_gradient(low = "black", high = "gray90") + nps_theme()) # for red: "darkred" and "#fee5d9"
    dev.off()
    
    # Scatterplots of IH coeffs
    IHcoeffs_long <- reshape(IHcoeffs, varying = names(IHcoeffs)[1:6], v.names = "IHcoeffs", timevar = "Variable", 
                             times = names(IHcoeffs)[1:6], direction = "long")
    jpeg(file=paste0(outLocationPath, "/", "IHACRES_Coeffs_Scatter.jpg"), width=600, height=500)
    print(ggplot(IHcoeffs_long, aes(x = IHcoeffs, y = nseD)) + geom_point() +
      facet_wrap(~ Variable, scales = 'free') + nps_theme() +
      labs(title = 'IHACRES Coefficients', x='', y = 'Daily NSE'))
    dev.off()
  }
  
  
  ### Run optim() to find optimal IHACRES flow coefficients ###
  IHcoeffs <- tibble()
  IHACRES_parms <- c(qa=qa, qb=qb, sa=sa, sb=sb, va=va)
  optDaily <- optim(par = IHACRES_parms, fn = IHACRESFlow, method = "L-BFGS-B",
                    lower = IHACRES_lower, upper = IHACRES_upper, hessian=TRUE, control = list(fnscale = -1, factr = '1e-2'),
                    #these parameters are carried through to IHACRESFlow
                    q0=q0, s0=s0, v0=v0, DailyWB=DailyWB, meas_flow_daily_xts=meas_flow_daily_xts, cutoffYear=cutoffYear)
  
  elpTimeD <- Sys.time() - strtTimeD
  
  # redefine IHACRES variables from optimized run
  # EDIT THIS - unnecessarily complicated
  optValuesD <- data.frame(nseD = optDaily$value, t(as.matrix(optDaily$par)))#data.frame(t(data.frame(optDaily$par)))
  optValuesD$vb <- calc_vb(optValuesD$qa, optValuesD$qb, optValuesD$sa, optValuesD$sb, optValuesD$va)
  qa=optValuesD$qa; qb=optValuesD$qb; sa=optValuesD$sa; sb=optValuesD$sb; va=optValuesD$va; vb=optValuesD$vb
  
  # store and save results
  results<- data.frame(results, qa=optValuesD$qa, qb=optValuesD$qb, sa=optValuesD$sa, sb=optValuesD$sb, va=optValuesD$va, vb=optValuesD$vb, nseD=optDaily$value, elpTimeD=elpTimeD)
  saveRDS(IHcoeffs, file = paste0(outLocationPath, "/IHcoeffs_final.rds"))
  
  
  ### Plot optimized IHACRES values
  if(make_plots){
    jpeg(file=paste0(outLocationPath, "/", "IHACRES_Coeffs_Final.jpg")); par(mfrow = c(3,2))
    for (i in 2:7){
      val = optValuesD[1,i]
      plot(c(0, 1), c(0, 0), type = "n", xlab = "", ylab = "",
           main = colnames(optValuesD)[i], xlim = c(-0.2, 1.2), ylim = c(-1, 1))
      segments(0, 0, 1, 0, col = "black", lwd = 2)
      points(val, 0, col = "red", pch = 19, cex = 1.5)
      text(val, 0.3, sprintf("%.4f", val), col='red')
    }
    dev.off()
  }
  
  
  ### Rerun entire model with optimal variables ###
  DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect,
               slope, shade.coeff, jtemp, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
  DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
}






#######################################################################
### Run model without optimization ###
if(!optimization){
  # run model
  DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon)
  DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
  MeasMod<- MeasModWB(DailyDrain, meas_flow_mon, cutoffYear)
  nseM = NSE(sim = MeasMod$Mod, obs = MeasMod$Meas)
  results = data.frame(results, nseM =nseM, gw_add = gw_add, vfm =vfm,jrange =jrange,hock = hock, hockros =hockros,
                       dro =dro, mondro = mondro, aspect = aspect, slope = slope,
                       shade.coeff = shade.coeff, jtemp =jtemp, elpTimeM = NA)
  
  # store results
  IHcoeffs <- tibble()
  nseD <-  IHACRESFlow(c(qa=qa, qb=qb, sa=sa, sb=sb, va=va), q0, s0, v0, DailyWB, meas_flow_daily_xts, cutoffYear)
  results<- data.frame(results, IHcoeffs, elpTimeD=NA)
}


### save results as RDS files ###
saveRDS(results, file = paste0(outLocationPath, "/results.rds"))



#######################################################################
#######################################################################
##### MODEL PERFORMANCE ON HISTORICAL FLOW #####


#######################################################################
### Create dataframes for aggregations ###

# Daily aggregation with measured and modeled streamflow
hist_flow_daily <- merge(xts(with(DailyDrain, cbind(adj_runoff, Quick, Slow, Very_Slow, total)), order.by = as.Date(DailyDrain$date)), meas_flow_daily_xts)
hist_flow_daily <- hist_flow_daily[complete.cases(hist_flow_daily),]
colnames(hist_flow_daily) <- c("adj_runoff", "quick", "slow", "very slow", "Mod", "Meas")

# Monthly aggregation with measured and modeled streamflow
MeasMod <- MeasModWB(DailyDrain = DailyDrain, meas_flow_mon = meas_flow_mon, cutoffYear = cutoffYear)

# Annual aggregation with measured and modeled streamflow
hist_flow_ann <- as.data.frame(matrix(NA, nrow = nrow(apply.yearly(hist_flow_daily[,"adj_runoff"], sum)),
                                    ncol = ncol(hist_flow_daily), dimnames = list(c(), colnames(hist_flow_daily))))
for(i in 1:ncol(hist_flow_daily)){
  hist_flow_ann[,i] <- apply.yearly(hist_flow_daily[,i], sum)
}


#######################################################################
### Create plots ### - DOUBLE CHECK ALL OF THESE FOR ACCURACY/CONSISTENTCY

# scatterplot of Historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
# there are two trend lines in the scatter plot because the intercept is set to 0 in one and allowed to vary in the other
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_Scatter.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  meas <- coredata(hist_flow_daily$Meas); mod <- coredata(hist_flow_daily$Mod)
  plot(mod, meas, main ="Daily Average Streamflow", xlab = "Modeled Streamflow (mm)", ylab = "Measured Streamflow (mm)")
  text(0.25*par("usr")[2], 0.85*par("usr")[4], paste('NSE:',round(NSE(mod, meas),digits=2)), cex = 3)
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # monthly
  mod <- MeasMod$Mod; meas <- MeasMod$Meas
  plot(mod, meas, main ="Monthly Total Streamflow", xlab = "Modeled Streamflow (mm)", ylab = "Measured Streamflow (mm)")
  text(0.25*par("usr")[2], 0.85*par("usr")[4], paste('NSE:',round(NSE(mod, meas),digits=2)), cex = 3)
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # annual
  mod <- coredata(hist_flow_ann$Mod); meas <- coredata(hist_flow_ann$Meas)
  plot(mod,meas, main ="Annual Total Streamflow", xlab = "Modeled Streamflow (mm)", ylab = "Measured Streamflow (mm)")
  text(0.35*par("usr")[2], 0.85*par("usr")[4], paste('NSE:',round(NSE(mod, meas),digits=2)), cex = 3)
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  dev.off()
}

# time series plot of historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_TimeSeries.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  plot(hist_flow_daily[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Daily Streamflow (mm)", main = "Daily", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # monthly
  plot(xts(MeasMod[,c("Mod", "Meas")], order.by = ym(MeasMod$YrMon)), 
       type = "l", lwd = 2, xlab = "Date", ylab = "Monthly Sum Streamflow (mm)", main = "Monthly", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # annual
  plot(xts(hist_flow_ann[,c('Mod','Meas')], order.by=as.Date(index(apply.yearly(hist_flow_daily, sum)))), 
       type = "l", lwd = 2, xlab = "Date", ylab = "Annual Sum Streamflow (mm)", main = "Annual", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}

# log plot of historical measured vs modeled daily streamflow 
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_Daily_Log.jpg"), width=700, height=400)
  plot(hist_flow_daily[,c('Mod','Meas')], type = "l", log=TRUE, lwd = 2, xlab = "Date", ylab = "Streamflow (mm)", main = "Log Plot of Daily Historical Streamflow", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}



#######################################################################
#######################################################################
### FUTURE STREAMFLOW PROJECTIONS ###
if(future_analysis){
  source("Future_Analysis.R")
}

