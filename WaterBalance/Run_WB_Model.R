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
optimization_var = 'AET'
delayStart = TRUE 
incompleteMonths = FALSE 
GridMET = TRUE
fillLeapDays = TRUE 
future_analysis = TRUE
runFutureWB = TRUE
calcFutureWB = TRUE
userSetJTemp = FALSE 
make_plots = FALSE 
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

### Provide geographic data or pull shapefile from StreamStats database ###
if(provide_coords) {
  lat = 37.9; lon = -122.59; aoi <- st_read("")  # add path to shapefile
} else{
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon; elev <- coords$elev; aoi <- coords$geometry
}
region <- get_region(lat,lon)

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
k_c = 1; et_slope = 1; et_bias = 0
Soil.Init = SWC.Max; Snowpack.Init = 0; T.Base = 0  

# Get j_temp
if(!userSetJTemp){
  j.raster = raster(here('Data', "merged_jennings.tif"))
  jtemp = get_jtemp(lat = lat, lon= lon, j.raster = j.raster)}

# Water balance optimization lower and upper limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10,  jtemp = jtemp-0.5, k_c=0, et_slope=-5, et_bias=-5)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400,  jtemp = jtemp+0.5, k_c=5, et_slope=5, et_bias=5)


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
DailyETo <- get_et_point(startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ETo', dataPath)



#######################################################################
#######################################################################
### OPTIMIZATION ###

### Optimize water balance variables according to the NSE of AET/openET over historical period ###
# need new function
if(optimization) source('WaterBalance//WB_Optimization.R')


#######################################################################
#######################################################################
### Run model  ###
DailyWB <- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, k_c, et_slope, et_bias, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, optimization_var)



#######################################################################
#######################################################################
##### MODEL PERFORMANCE ###

#######################################################################
### HISTORICAL AET #####

# Daily, monthly, annual aggregation with measured and modeled ET
hist_aet_daily <- merge(xts(with(DailyWB, cbind(AET)), order.by = as.Date(DailyWB$date)), DailyET); hist_aet_daily <- hist_aet_daily[complete.cases(hist_aet_daily),]
colnames(hist_aet_daily) <- c("Mod", "Meas")
hist_aet_monthly <- apply.monthly(xts(with(DailyWB, cbind(AET)), order.by = as.Date(DailyWB$date)), function(x) {colSums(x, na.rm = TRUE)})
index(hist_aet_monthly) <- as.Date(format(index(hist_aet_monthly), "%Y-%m-01"))
hist_aet_monthly <- merge(hist_aet_monthly, MonthlyET); hist_aet_monthly <- hist_aet_monthly[complete.cases(hist_aet_monthly),]; colnames(hist_aet_monthly) <- c("Mod", "Meas")
hist_aet_ann <- apply.yearly(hist_aet_monthly, function(x) {colSums(x, na.rm = TRUE)})


#######################################################################
### Create summary plots ###

# scatterplot of Historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
# there are two trend lines in the scatter plot because the intercept is set to 0 in one and allowed to vary in the other
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_AET_Scatter.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  meas <- coredata(hist_aet_daily$Meas); mod <- coredata(hist_aet_daily$Mod)
  plot(mod, meas, main =paste("Daily Average AET\nNSE:",round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # monthly
  mod <- coredata(hist_aet_monthly$Mod); meas <- coredata(hist_aet_monthly$Meas)
  plot(mod, meas, main =paste("Monthly Total AET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # annual
  mod <- coredata(hist_aet_ann$Mod); meas <- coredata(hist_aet_ann$Meas)
  plot(mod,meas, main =paste("Annual Total AET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  dev.off()
}

# time series plot of historical Measured vs Modeled AET for daily, monthly, annual aggregation
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_AET_TimeSeries.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  plot(hist_aet_daily[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Daily AET (mm)", main = "Daily", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # monthly
  plot(hist_aet_monthly[,c("Mod", "Meas")], type = "l", lwd = 2, xlab = "Date", ylab = "Monthly Sum AET (mm)", main = "Monthly", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # annual
  plot(hist_aet_ann[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Annual Sum AET (mm)", main = "Annual", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}

#######################################################################
### HISTORICAL PET #####
# ETo from OpenET is for a grass reference crop, so isn't accurate in areas 
# with more complex vegetation cover. These results should not be considered
# unless the model calibration has been PET to ETo.

# Daily, monthly, annual aggregation with measured and modeled ET
hist_pet_daily <- merge(xts(with(DailyWB, cbind(PET)), order.by = as.Date(DailyWB$date)), DailyETo); hist_pet_daily <- hist_pet_daily[complete.cases(hist_pet_daily),]
colnames(hist_pet_daily) <- c("Mod", "Meas")
hist_pet_monthly <- apply.monthly(xts(with(DailyWB, cbind(PET)), order.by = as.Date(DailyWB$date)), function(x) {colSums(x, na.rm = TRUE)})
index(hist_pet_monthly) <- as.Date(format(index(hist_pet_monthly), "%Y-%m-01"))
hist_pet_monthly <- merge(hist_pet_monthly, MonthlyETo); hist_pet_monthly <- hist_pet_monthly[complete.cases(hist_pet_monthly),]; colnames(hist_pet_monthly) <- c("Mod", "Meas")
hist_pet_ann <- apply.yearly(hist_pet_monthly, function(x) {colSums(x, na.rm = TRUE)})


#######################################################################
### Create summary plots ###

# scatterplot of Historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
# there are two trend lines in the scatter plot because the intercept is set to 0 in one and allowed to vary in the other
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_PET_Scatter.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  meas <- coredata(hist_pet_daily$Meas); mod <- coredata(hist_pet_daily$Mod)
  plot(mod, meas, main =paste("Daily Average PET\nNSE:",round(NSE(mod, meas),digits=2)), xlab = "Modeled PET (mm)", ylab = "Measured PET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # monthly
  mod <- coredata(hist_pet_monthly$Mod); meas <- coredata(hist_pet_monthly$Meas)
  plot(mod, meas, main =paste("Monthly Total PET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled PET (mm)", ylab = "Measured PET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # annual
  mod <- coredata(hist_pet_ann$Mod); meas <- coredata(hist_pet_ann$Meas)
  plot(mod,meas, main =paste("Annual Total PET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled PET (mm)", ylab = "Measured PET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  dev.off()
}

# time series plot of historical Measured vs Modeled PET for daily, monthly, annual aggregation
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_PET_TimeSeries.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  plot(hist_pet_daily[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Daily PET (mm)", main = "Daily", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # monthly
  plot(hist_pet_monthly[,c("Mod", "Meas")], type = "l", lwd = 2, xlab = "Date", ylab = "Monthly Sum PET (mm)", main = "Monthly", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # annual
  plot(hist_pet_ann[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Annual Sum PET (mm)", main = "Annual", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}








