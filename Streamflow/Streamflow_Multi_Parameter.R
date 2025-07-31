# ---------------------------------------------------------------------
# This script contains code for running the streamflow model
# for a single watershed with multiple different parameters,
# without additional user-provided input. This script is
# intended to work as a "manual" way to investigate the impacts of 
# changing individual parameters on model calibration, accuracy, and results. 
# ---------------------------------------------------------------------

### Load libraries and function files ###
library('here'); lib_install <- FALSE
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
#######################################################################
### DEFINE WATERSHED ###

### Watershed name and USGS gage number ###
SiteID = "Wet Beaver Creek"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '09505200'

### Model parameters ###
provide_coords = FALSE

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

### Model names for future analysis ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')


#######################################################################
#######################################################################
### SET MODEL PARAMETERS ###
make_plots = TRUE; historical_analysis = FALSE; future_analysis = FALSE

# Water balance optimization limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400)

# Store all parameters in tibble
v = tibble(
  # Folder name
  FolderName = "optim_Oudin",
  
  # General model variables
  optimization = TRUE, point_location = FALSE, percent_skill_cutoff = 0.1, 
  
  # Meteorological variables
  GridMET = TRUE, incompleteMonths = FALSE, fillLeapDays = TRUE, delayStart = FALSE,
  tmmx_slope = 1, tmmx_bias = 0, tmmn_slope = 1, tmmn_bias = 0, p_slope = 1, p_bias = 0,
  
  # Water balance variables
  PETMethod = "Oudin",  calcFutureWB = TRUE, userSetJTemp = FALSE,
  gw_add=0, vfm = 0.7555, jtemp = 1.982841, jrange = 3 ,hock = 4, hockros = 4, 
  dro = 0, mondro = 0, aspect = 180, slope= 0, shade.coeff= 1, SWC.Max = 200,
  Soil.Init = SWC.Max, Snowpack.Init = 0, T.Base = 0,
  
  # Streamflow variables
  flow_components = 3, NonZeroDrainInitCoeff = FALSE,
  qa=0.62, qb=0.22, sa=0.58, sb=0.06, va=0.974, vb=calc_vb(qa,qb,sa,sb,va)
)

### Create second row of variables and change as desired ###
v[2,] <- v[1,]; v[3,] <- v[1,]
v[2,"PETMethod"] = "Penman"; v[2,"FolderName"] = "optim_Penman" 
v[3,"PETMethod"] = "Hamon"; v[3,"FolderName"] = "optim_Hamon" 


#######################################################################
#######################################################################
### RUN MODEL FOR EACH ITERATION OF VARIABLES ###

# Create extra column to save model performance
v$nseM <- NA_real_

for(j in 1:nrow(v)){
  
  #######################################################################
  ### Define variables ###
  
  ### Read in from v tibble
  FolderName=v$FolderName[j]; 
  optimization=v$optimization[j]; point_location=v$point_location[j]; percent_skill_cutoff=v$percent_skill_cutoff[j] 
  GridMET=v$GridMET[j]; incompleteMonths=v$incompleteMonths[j]; fillLeapDays=v$fillLeapDays[j]; delayStart=v$delayStart[j]
  tmmx_slope=v$tmmx_slope[j]; tmmx_bias=v$tmmx_bias[j]; tmmn_slope=v$tmmn_slope[j]; tmmn_bias=v$tmmn_bias[j]; p_slope=v$p_slope[j]; p_bias=v$p_bias[j]
  PETMethod=v$PETMethod[j]; calcFutureWB=v$calcFutureWB[j]; userSetJTemp=v$userSetJTemp[j]
  gw_add=v$gw_add[j]; vfm=v$vfm[j]; jtemp=v$jtemp[j]; jrange=v$jrange[j]; hock=v$hock[j]; hockros=v$hockros[j] 
  dro=v$dro[j]; mondro=v$mondro[j]; aspect=v$aspect[j]; slope=v$slope[j]; shade.coeff=v$shade.coeff[j]; SWC.Max=v$SWC.Max[j]
  Soil.Init=v$Soil.Init[j]; Snowpack.Init=v$Snowpack.Init[j]; T.Base=v$T.Base
  flow_components=v$flow_components[j]; NonZeroDrainInitCoeff=v$NonZeroDrainInitCoeff[j]
  qa=v$qa[j]; qb=v$qb[j]; sa=v$sa[j]; sb=v$sb[j]; va=v$va[j]; vb=v$vb[j]
  
  ### Define other parameters
  if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 
  if(!userSetJTemp){
    j.raster = raster(here('Data', "merged_jennings.tif"))
    jtemp = get_jtemp(lat = lat, lon= lon, j.raster = j.raster)}
  WB_lower['jtemp'] <- jtemp-0.5; WB_upper['jtemp'] <- jtemp+0.5
  
  ### Set path variables
  if(!dir.exists(here('Data', SiteID_FileName))) {dir.create(here('Data', SiteID_FileName))}; dataPath <- here('Data', SiteID_FileName)
  if(!dir.exists(here('Output', SiteID_FileName))) {dir.create(here('Output', SiteID_FileName))}
  if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow'))) {dir.create(here('Output', SiteID_FileName, 'Streamflow'))}
  if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow', FolderName))) {dir.create(here('Output', SiteID_FileName, 'Streamflow', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'Streamflow', FolderName)
  
  #######################################################################
  ### Get data ###
  
  # USGS gage data
  gage_data <- get_gage_data(GageSiteID, incompleteMonths, fillLeapDays, dataPath)
  meas_flow_daily <- gage_data$meas_flow_daily; meas_flow_mon <- gage_data$meas_flow_mon
  
  # Historical meteorological data
  if(GridMET) {
    if(point_location){DailyClimData <- get_gridmet_point(SiteID_FileName, startY, endY, lat, lon, dataPath,
                                                          tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope)
    } else {DailyClimData <- get_gridmet_area(SiteID_FileName, startY, endY, aoi, dataPath,
                                              tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope)}
  } else { 
    if(point_location){DailyClimData <- get_daymet_point(SiteID_FileName, startY, endY, lat, lon, dataPath)
    } else{DailyClimData <- get_daymet_area(SiteID_FileName, startY, endY, aoi, dataPath)}
  }
  
  #######################################################################
  # Run model
  
  ### Get initial flow conditions or set all to 0
  if(NonZeroDrainInitCoeff){
    InitCond <- get_Init_Drain_Coef(DailyClimData, gw_add, vfm , jrange ,hock, hockros, dro, mondro , aspect, slope, shade.coeff, jtemp, SWC.Max, 
                                    Soil.Init, Snowpack.Init, T.Base, PETMethod,q0, s0, v0, qa, qb, sa, sb, va, vb, lat, lon, cutoffYear)
    q0 = InitCond[["Quick"]]; s0= InitCond[["Slow"]]; v0 = InitCond[["Very_Slow"]]
  } else{
    q0<-0; s0<-0; v0<-0
  }
  
  # Call optimization routine
  if(optimization){
    source('Streamflow//Streamflow_Optimization.R')
  }
  
  # Run model
  DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, 1, 1, 0, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, "")
  DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
  MeasMod<- MeasModWB(DailyDrain, meas_flow_mon, cutoffYear)
  
  # Calculate statistics
  v$nseM[j] <- NSE(MeasMod$Mod, MeasMod$Meas)
  
  #######################################################################
  # Analysis
  
  # Model accuracy
  source('Streamflow//Streamflow_Model_Accuracy.R')
  
  # Historical streamflow analysis
  if(historical_analysis){
    source('Streamflow//Streamflow_Historical_Analysis.R')
  }
  
  # Future streamflow projections
  if(future_analysis){
    source("Streamflow//Streamflow_Future_Analysis.R")
  }
}

# Save results of runs
saveRDS(v, file = paste0(outLocationPath, "/MultiParameter_Results.rds"))