# ---------------------------------------------------------------------
# This script contains code for running the water balance model
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
SiteID = "Redwood Creek"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '11460151'

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
startY = 2000; startM = 01; startD = 01 
endY = 2024; endM = 12; endD = 31

### Model names for future analysis ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')


#######################################################################
#######################################################################
### SET MODEL PARAMETERS ###
make_plots = TRUE; historical_analysis = FALSE; future_analysis = FALSE

# Water balance optimization limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10, k_c=0, et_slope=-5, et_bias=-5)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400, k_c=5, et_slope=5, et_bias=5)

# Store all parameters in tibble
v = tibble(
  # Folder name
  FolderName = "optim_Oudin",
  
  # General model variables
  optimization = TRUE, optimization_var='AET', point_location = TRUE, percent_skill_cutoff = 0.1, 
  
  # Meteorological variables
  GridMET = TRUE, incompleteMonths = FALSE, fillLeapDays = TRUE, delayStart = FALSE,
  tmmx_slope = 1, tmmx_bias = 0, tmmn_slope = 1, tmmn_bias = 0, p_slope = 1, p_bias = 0,
  
  # Water balance variables
  PETMethod = "Oudin",  calcFutureWB = TRUE, userSetJTemp = FALSE,
  gw_add=0, vfm = 0.7555, jtemp = 1.982841, jrange = 3 ,hock = 4, hockros = 4, 
  dro = 0, mondro = 0, aspect = 180, slope= 0, shade.coeff= 1, SWC.Max = 200,
  k_c = 1, et_slope = 1, et_bias = 0,
  Soil.Init = SWC.Max, Snowpack.Init = 0, T.Base = 0,
)

### Create second row of variables and change as desired ###
v[2,] <- v[1,]; v[3,] <- v[1,]
v[2,"PETMethod"] = "Penman"; v[2,"FolderName"] = "optim_Penman" 
v[3,"PETMethod"] = "Hamon"; v[3,"FolderName"] = "optim_Hamon" 


#######################################################################
#######################################################################
### RUN MODEL FOR EACH ITERATION OF VARIABLES ###

# Create extra column to save model performance
v$nseM <- NA_real_; v$nseD<- NA_real_; 

for(j in 1:nrow(v)){
  
  #######################################################################
  ### Define variables ###
  
  ### Read in from v tibble
  FolderName=v$FolderName[j]; 
  optimization=v$optimization[j]; optimization_var=v$optimization_var[j]; point_location=v$point_location[j]; percent_skill_cutoff=v$percent_skill_cutoff[j] 
  GridMET=v$GridMET[j]; incompleteMonths=v$incompleteMonths[j]; fillLeapDays=v$fillLeapDays[j]; delayStart=v$delayStart[j]
  tmmx_slope=v$tmmx_slope[j]; tmmx_bias=v$tmmx_bias[j]; tmmn_slope=v$tmmn_slope[j]; tmmn_bias=v$tmmn_bias[j]; p_slope=v$p_slope[j]; p_bias=v$p_bias[j]
  PETMethod=v$PETMethod[j]; calcFutureWB=v$calcFutureWB[j]; userSetJTemp=v$userSetJTemp[j]
  gw_add=v$gw_add[j]; vfm=v$vfm[j]; jtemp=v$jtemp[j]; jrange=v$jrange[j]; hock=v$hock[j]; hockros=v$hockros[j] 
  dro=v$dro[j]; mondro=v$mondro[j]; aspect=v$aspect[j]; slope=v$slope[j]; shade.coeff=v$shade.coeff[j]; SWC.Max=v$SWC.Max[j]
  k_c=v$k_c[j]; et_slope=v$et_slope[j]; et_bias=v$et_bias[j]
  Soil.Init=v$Soil.Init[j]; Snowpack.Init=v$Snowpack.Init[j]; T.Base=v$T.Base
  
  ### Define other parameters
  if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 
  if(!userSetJTemp){
    j.raster = raster(here('Data', "merged_jennings.tif"))
    jtemp = get_jtemp(lat = lat, lon= lon, j.raster = j.raster)}
  WB_lower['jtemp'] <- jtemp-0.5; WB_upper['jtemp'] <- jtemp+0.5
  
  ### Set path variables
  if(!dir.exists(here('Data', SiteID_FileName))) {dir.create(here('Data', SiteID_FileName))}; dataPath <- here('Data', SiteID_FileName)
  if(!dir.exists(here('Output', SiteID_FileName))) {dir.create(here('Output', SiteID_FileName))}
  if(!dir.exists(here('Output', SiteID_FileName, 'WaterBalance'))) {dir.create(here('Output', SiteID_FileName, 'WaterBalance'))}
  if(!dir.exists(here('Output', SiteID_FileName, 'WaterBalance', FolderName))) {dir.create(here('Output', SiteID_FileName,'WaterBalance', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'WaterBalance', FolderName)
  
  #######################################################################
  ### Get data ###
  
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
  
  # OpenET data
  if(point_location){
    MonthlyET <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ET', dataPath)
    DailyET <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ET', dataPath)
    
    MonthlyETo <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ETo', dataPath)
    DailyETo <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ETo', dataPath)
  } else{
    MonthlyET <- get_et_area(aoi, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ET', dataPath)
    DailyET <- get_et_area(aoi, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ET', dataPath)
    
    MonthlyETo <- get_et_area(aoi, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ETo', dataPath)
    DailyETo <- get_et_area(aoi, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ETo', dataPath)
    
  }
 
  
  #######################################################################
  # Run model
  
  # Call optimization routine
  if(optimization) source('WaterBalance//WB_Optimization.R')
  
  # Run model
  DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, k_c, et_slope, et_bias, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, "")
  
  # Calculate statistics (AET)
  Mod <- DailyWB %>% mutate(month = floor_date(as.Date(date), "month")) %>%  group_by(month) %>%
      summarise(Mod = sum(AET, na.rm = TRUE)) %>% rename(date = month)
  MeasMod <- dplyr::full_join(MonthlyET, Mod, by = join_by(date)); colnames(MeasMod)<- c("date", "Meas","Mod")
  MeasMod<- MeasMod[complete.cases(MeasMod),]
  v$nseM[j] <- NSE(MeasMod$Mod, MeasMod$Meas)
    
  Mod <- DailyWB %>% select(date, AET)
  MeasMod <- dplyr::full_join(DailyET, Mod, by = join_by(date)); colnames(MeasMod)<- c("date", "Meas","Mod")
  MeasMod<- MeasMod[complete.cases(MeasMod),]
  v$nseD[j] <- NSE(MeasMod$Mod, MeasMod$Meas)
  
  #######################################################################
  # Analysis
  
  # Model accuracy
  source('WaterBalance/WB_Model_Accuracy.R')
  
  # Future streamflow projections
  if(future_analysis) source("WaterBalance//WB_Future_Analysis.R")
  
}

# Save results of runs
saveRDS(v, file = paste0(outLocationPath, "/MultiParameter_Results.rds")) 