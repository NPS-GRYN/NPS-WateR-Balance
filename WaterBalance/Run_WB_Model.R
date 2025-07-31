# ---------------------------------------------------------------------
# This script includes code for running and calibrating the water balance model. The model can be run with
# user-provided model paramters, or the parameters can be calibrated based on measured ET values from 
# the OpenET dataset. 
#
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
calcFutureWB = FALSE
userSetJTemp = FALSE 
make_plots = TRUE 
provide_coords = FALSE
point_location = FALSE
percent_skill_cutoff = 0.1 
FolderName = 'optim' 

### Define watershed ###
# must be west of approximately -90 longitude to have OpenET data
# USGS gage number: https://waterdata.usgs.gov/nwis/rt
SiteID = "Redwood Creek"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '11460151'      

### Provide geographic data or pull shapefile from StreamStats database ###
if(provide_coords) {
  lat = 37.9; lon = -122.59; aoi <- st_read("")  # add path to shapefile
} else{
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon; elev <- coords$elev; aoi <- coords$geometry
}
region <- get_region(lat,lon)

### Define time period for historical analysis ###
startY = 2000; startM = 01; startD = 01 
endY = 2024; endM = 12; endD = 31


### Model names for future streamflow ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')


#######################################################################
### Set model variables ###

# Read in previously optimized variables from results file, if they exist
if(file.exists(paste0(outLocationPath, "/optim_results.rds"))){
  results <- readRDS(paste0(outLocationPath, "/optim_results.rds"))
  gw_add <- results$gw_add; vfm <- results$vfm; jrange <- results$jrange; hock <- results$hock; hockros <- results$hockros; dro <- results$dro; 
  mondro <- results$mondro; aspect <- results$aspect; slope <- results$slope; shade.coeff <- results$shade.coeff; SWC.Max <- results$SWC.Max; 
  k_c <- results$k_c; et_slope <- results$et_slope; et_bias <- results$et_bias
  
  # If not, manually set variables
} else {
  # Default water balance variables
  gw_add=0; vfm = 0.7555; jtemp = 1.982841; jrange = 3 ;hock = 4; hockros = 4; 
  dro = 0; mondro = 0; aspect = 180; slope= 0; shade.coeff= 1; SWC.Max = 200
  k_c = 1; et_slope = 1; et_bias = 0
}

# Non-optimized WB variables
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

MonthlyET <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ET', dataPath)
DailyET <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ET', dataPath)

MonthlyETo <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ETo', dataPath)
DailyETo <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ETo', dataPath)



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
### ANALYSIS ###

### MODEL PERFORMANCE ###
source('WB/WB_Model_Accuracy.R')


### FUTURE ANALYSIS ###
if(future_analysis) source("WB/WB_Future_Analysis.R")



