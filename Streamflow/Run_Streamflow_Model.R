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
# add code for watershed averages, not just centroids
# CONSISTENCY IN NAMING THINGS: projection and date are all lowercase
# ---------------------------------------------------------------------


#######################################################################
#######################################################################
### INTRODUCTION SECTION ###


#######################################################################
### Load libraries and function files ###
library('here'); lib_install <- FALSE
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
### Set user-defined variables ###
PETMethod = "Oudin" 
optimization = FALSE 
delayStart = FALSE 
NonZeroDrainInitCoeff = FALSE
incompleteMonths = FALSE 
GridMET = TRUE
fillLeapDays = TRUE 
historical_analysis = TRUE
future_analysis = TRUE
calcFutureWB = TRUE  
userSetJTemp = FALSE 
make_plots = TRUE 
provide_coords = FALSE
point_location = FALSE
flow_components = 3
percent_skill_cutoff = 0.1 
FolderName = "optim" 
#filename_future_wb = "\\Users\\mcburns\\OneDrive - DOI\\water-balance\\Data\\LittleRiver\\littleriver_water_balance_future.csv"

### Define watershed ###
SiteID = "Wet Beaver Creek"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '09505200'                  #define stream gage location (RWC: "11460151", LR: 03497300, C: 03460000)

### Provide geographic data or pull shapefile from StreamStats database ###
if(provide_coords) {
  lat = 37.9; lon = -122.59; aoi <- st_read("")  # add path to shapefile
} else{
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon; elev <- coords$elev; aoi <- coords$geometry
}
region <- get_region(lat,lon)

### Define time period for historical analysis ###
# GridMET begins in 1970, Daymet begins in 1980 
startY = 1979; startM = 01; startD = 01 
endY = 2023; endM = 12; endD = 31
if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 


### Model names for future streamflow ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')



### Set path variables ###
if(!dir.exists(here('Data', SiteID_FileName))) {dir.create(here('Data', SiteID_FileName))}; dataPath <- here('Data', SiteID_FileName)
if(!dir.exists(here('Output', SiteID_FileName))) {dir.create(here('Output', SiteID_FileName))}
if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow'))) {dir.create(here('Output', SiteID_FileName, 'Streamflow'))}
if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow', FolderName))) {dir.create(here('Output', SiteID_FileName, 'Streamflow', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'Streamflow', FolderName)


#######################################################################
### Set model variables ###

# Read in previously optimized variables from results file, if they exist
if(file.exists(paste0(outLocationPath, "/optim_results.rds"))){
  results <- readRDS(paste0(outLocationPath, "/optim_results.rds"))
  gw_add <- results$gw_add; vfm <- results$vfm; jrange <- results$jrange; hock <- results$hock; hockros <- results$hockros; dro <- results$dro; 
  mondro <- results$mondro; aspect <- results$aspect; slope <- results$slope; shade.coeff <- results$shade.coeff; SWC.Max <- results$SWC.Max; 
  jtemp <- results$jtemp; qa <- results$qa; qb <- results$qb; sa <- results$sa; sb <- results$sb; va <- results$va; vb <- results$vb

# If not, manually set variables
} else {
  # Default IHACRES flow coefficients
  if(flow_components==3){
    qa<-0.62; qb<-0.22; sa<-0.58; sb<-0.06; va<-0.974; vb<-calc_vb(qa,qb,sa,sb,va)
  }else if(flow_components==2){
    qa<-0; qb<-0; sa<-0.58; sb<-0.06; va<-0.974; vb<-calc_vb(qa,qb,sa,sb,va) 
  } else{print('invalid number of flow components')}
  
  # Default water balance variables
  gw_add=0; vfm = 0.7555; jtemp = 1.982841; jrange = 3 ;hock = 4; hockros = 4; 
  dro = 0; mondro = 0; aspect = 180; slope= 0; shade.coeff= 1; SWC.Max = 200
}

# Get j_temp
if(!userSetJTemp){
  j.raster = raster(here('Data', "merged_jennings.tif"))
  jtemp = get_jtemp(lat = lat, lon= lon, j.raster = j.raster)}

# Non-optimized WB variables
Soil.Init = SWC.Max; Snowpack.Init = 0; T.Base = 0 

# Water balance optimization limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10, jtemp = jtemp-0.5)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400, jtemp = jtemp+0.5)



#######################################################################
### Set other variables ###
# Optional scaling factors for GridMET time series: if no scaling, set slopes to 1 and bias to 0
tmmx_slope = 1; tmmx_bias = 0
tmmn_slope = 1; tmmn_bias = 0
p_slope = 1; p_bias = 0

if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 


#######################################################################
#######################################################################
### GET DATA ###

#######################################################################
### Scrape and clean USGS stream gage data ###

gage_data <- get_gage_data(GageSiteID, incompleteMonths, fillLeapDays, dataPath)
meas_flow_daily <- gage_data$meas_flow_daily; meas_flow_mon <- gage_data$meas_flow_mon


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
#######################################################################
### MODEL RUNNING ### 

### Get initial flow conditions or set all to 0 ###
if(NonZeroDrainInitCoeff){
  InitCond <- get_Init_Drain_Coef(DailyClimData, gw_add, vfm , jrange ,hock, hockros,
                                  dro, mondro , aspect, slope, shade.coeff, jtemp, SWC.Max, 
                                  Soil.Init, Snowpack.Init, T.Base, PETMethod, 
                                  q0, s0, v0, qa, qb, sa, sb, va, vb, lat, lon, cutoffYear)
  
  q0 = InitCond[["Quick"]]; s0= InitCond[["Slow"]]; v0 = InitCond[["Very_Slow"]]
} else{
  q0<-0; s0<-0; v0<-0
}


### Call optimization routine to get optimal variables ###
if(optimization){
  source('Streamflow//Streamflow_Optimization.R')
}


### Run model ###
DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
             shade.coeff, jtemp,SWC.Max, 1, 1, 0, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, "")
DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
MeasMod<- MeasModWB(DailyDrain, meas_flow_mon, cutoffYear)



#######################################################################
#######################################################################
### ANALYSIS ###

### MODEL PERFORMANCE ON HISTORICAL FLOW ###
source('Streamflow//Streamflow_Model_Accuracy.R')


### HISTORICAL STREAMFLOW ANALYSIS ###
if(historical_analysis){
  source('Streamflow//Streamflow_Historical_Analysis.R')
}


### FUTURE STREAMFLOW PROJECTIONS ###
if(future_analysis){
  source("Streamflow//Streamflow_Future_Analysis.R")
}

