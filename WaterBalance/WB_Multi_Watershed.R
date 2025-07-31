# ---------------------------------------------------------------------
# This script contains code for running and calibrating the water balance model
# across multiple watersheds without additional user-provided input. 
#
# ---------------------------------------------------------------------

### Load libraries and function files ###
library('here'); lib_install <- FALSE
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
#######################################################################
### WATERSHED NAMES AND USGS GAGE NUMBERS ###
watershed_siteid <- c("Redwood Creek", "Wet Beaver Creek")
watershed_gageid <- c("11460151", "09505200")
watershed_foldername <- c('optim', 'optim')
num_watersheds <- length(watershed_siteid)

#######################################################################
#######################################################################
### SET PARAMETERS FOR ALL WATERSHEDS ###
# These parameters will be used for all watersheds
# To experiement with the effects of parameter changes on calibration, use
# the WB_Multi_Parameter.R script.
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

### Define time period for historical analysis ###
startY = 2000; startM = 01; startD = 01 
endY = 2024; endM = 12; endD = 31
if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 

### Model names for future streamflow ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')

### Default model parameters ###
# Water balance variables
gw_add=0; vfm = 0.7555; jtemp = 1.982841; jrange = 3 ;hock = 4; hockros = 4; 
dro = 0; mondro = 0; aspect = 180; slope= 0; shade.coeff= 1; SWC.Max = 200
k_c = 1; et_slope = 1; et_bias = 0

# Non-optimized WB variables
Soil.Init = SWC.Max; Snowpack.Init = 0; T.Base = 0 

# Water balance optimization limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10,  jtemp = jtemp-0.5, k_c=0, et_slope=-5, et_bias=-5)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400,  jtemp = jtemp+0.5, k_c=5, et_slope=5, et_bias=5)

### Optional scaling factors for GridMET ###
# if no scaling, set slopes to 1 and bias to 0
tmmx_slope = 1; tmmx_bias = 0
tmmn_slope = 1; tmmn_bias = 0
p_slope = 1; p_bias = 0

#######################################################################
#######################################################################
### RUN MODEL FOR EACH WATERSHED ###

for(i in 1:num_watersheds){
  #######################################################################
  ### Establish preliminary variables ###
  SiteID = watershed_siteid[i]; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
  GageSiteID <- watershed_gageid[i]
  FolderName <- watershed_foldername[i]
  
  # Set paths
  if(!dir.exists(here('Data', SiteID_FileName))) {dir.create(here('Data', SiteID_FileName))}; dataPath <- here('Data', SiteID_FileName)
  if(!dir.exists(here('Output', SiteID_FileName))) {dir.create(here('Output', SiteID_FileName))}
  if(!dir.exists(here('Output', SiteID_FileName, 'WaterBalance'))) {dir.create(here('Output', SiteID_FileName, 'WaterBalance'))}
  if(!dir.exists(here('Output', SiteID_FileName, 'WaterBalance', FolderName))) {dir.create(here('Output', SiteID_FileName,'WaterBalance', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'WaterBalance', FolderName)
  
  # Get watershed area and coordinates
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon; elev <- coords$elev; aoi <- coords$geometry
  region <- get_region(lat,lon)
  
  # Read in previously optimized variables from results file, if they exist
  if(file.exists(paste0(outLocationPath, "/optim_results.rds"))){
    results <- readRDS(paste0(outLocationPath, "/optim_results.rds"))
    gw_add <- results$gw_add; vfm <- results$vfm; jrange <- results$jrange; hock <- results$hock; hockros <- results$hockros; dro <- results$dro; 
    mondro <- results$mondro; aspect <- results$aspect; slope <- results$slope; shade.coeff <- results$shade.coeff; SWC.Max <- results$SWC.Max; 
    jtemp <- results$jtemp; qa <- results$qa; qb <- results$qb; sa <- results$sa; sb <- results$sb; va <- results$va; vb <- results$vb
  }
  
  # Get j_temp
  if(!userSetJTemp){
    j.raster = raster(here('Data', "merged_jennings.tif"))
    jtemp = get_jtemp(lat = lat, lon= lon, j.raster = j.raster)}
  
  
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
  MonthlyET <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ET', dataPath)
  DailyET <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ET', dataPath)
  
  MonthlyETo <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'monthly', 'ETo', dataPath)
  DailyETo <- get_et_point(lat, lon, startY, startM, startD, endY, endM, endD, SiteID_FileName, 'daily', 'ETo', dataPath)
  
  
  #######################################################################
  # Run model
  
  # Call optimization routine
  if(optimization) source('WaterBalance//WB_Optimization.R')
  
  # Run model
  DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
               shade.coeff, jtemp,SWC.Max, k_c, et_slope, et_bias, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, "")
  
  
  #######################################################################
  # Analysis
  
  # Model accuracy
  source('WaterBalance/WB_Model_Accuracy.R')
  
  # Future streamflow projections
  if(future_analysis) source("WaterBalance//WB_Future_Analysis.R")
}

