# ---------------------------------------------------------------------
# This script contains code for running and calibrating the streamflow model
# across multiple watersheds without additional user-provided input. 
# ---------------------------------------------------------------------

### Load libraries and function files ###
library('here'); lib_install <- FALSE
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
#######################################################################
### WATERSHED NAMES AND USGS GAGE NUMBERS ###
watershed_siteid <- c("Little River", "Cataloochee")
watershed_gageid <- c("03497300", "03460000")
watershed_foldername <- c('optim', 'optim')
num_watersheds <- length(watershed_siteid)


#######################################################################
#######################################################################
### SET PARAMETERS FOR ALL WATERSHEDS ###
# These parameters will be used for all watersheds
# To experiement with the effects of parameter changes on calibration, use
# the Streamflow_Multi_Parameter.R script.

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
point_location = FALSE
flow_components = 3
percent_skill_cutoff = 0.1 

### Define time period for historical analysis ###
# GridMET begins in 1970, Daymet begins in 1980 
startY = 1979; startM = 01; startD = 01 
endY = 2023; endM = 12; endD = 31
if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 

### Model names for future streamflow ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')

### Default model parameters ###
# IHACRES flow coefficients
if(flow_components==3){
  qa<-0.62; qb<-0.22; sa<-0.58; sb<-0.06; va<-0.974; vb<-calc_vb(qa,qb,sa,sb,va)
}else if(flow_components==2){
  qa<-0; qb<-0; sa<-0.58; sb<-0.06; va<-0.974; vb<-calc_vb(qa,qb,sa,sb,va) 
} else{print('invalid number of flow components')}

# Water balance variables
gw_add=0; vfm = 0.7555; jtemp = 1.982841; jrange = 3 ;hock = 4; hockros = 4; 
dro = 0; mondro = 0; aspect = 180; slope= 0; shade.coeff= 1; SWC.Max = 200

# Non-optimized WB variables
Soil.Init = SWC.Max; Snowpack.Init = 0; T.Base = 0 

# Water balance optimization limits
WB_lower = c(gw_add=0, vfm = 0.25, jrange = 1, hock = 0.25, hockros = 0.25, dro= 0, mondro = 0, aspect= 0, slope =  0, shade.coeff = 0.1, SWC.Max = 10, jtemp = jtemp-0.5)
WB_upper = c(gw_add = 1, vfm = 1, jrange = 5, hock = 8, hockros = 8, dro = 1, mondro = 1, aspect = 360, slope = 90, shade.coeff = 1, SWC.Max = 400, jtemp = jtemp+0.5)

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
  if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow'))) {dir.create(here('Output', SiteID_FileName, 'Streamflow'))}
  if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow', FolderName))) {dir.create(here('Output', SiteID_FileName, 'Streamflow', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'Streamflow', FolderName)
  
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

