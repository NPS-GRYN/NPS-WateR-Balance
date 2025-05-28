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
# ---------------------------------------------------------------------


#######################################################################
#######################################################################
### INTRODUCTION SECTION ###


#######################################################################
### Load libraries ###
library('here'); lib_install <- FALSE

### Source in function files ###
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
calcFutureWB = FALSE  # TRUE to re-run entire water balance model for future; FALSE to use pre-existing CONUS water balance projections from a Mike Tercek spreadsheet
userSetJTemp = FALSE 
make_plots = TRUE 
provide_coords = FALSE # if true, user provides lat/lon coords. if false, lat/lon coords are pulled from centroid of watershed with given gage id
flow_components = 3  # change the number of components that characterize the flow. can be 2 or 3. 2: flow has quick and slow components; 3: flow has quick, slow, and very slow components.
percent_skill_cutoff = 0.1 # only include a certain percentage of GCMs, ranked by skill. metric is 0-1
point_location = TRUE  # if TRUE, pull all meteorological data for a single point (center of watershed). if FALSE, pull data for entire watershed and take average (this takes longer)
FolderName = "optim" 

### Define watershed ###
# centroid of watershed
SiteID = "Cataloochee"; SiteID_FileName = gsub(pattern = " ", x = SiteID, replacement = "")
GageSiteID <- '03460000'                  #define stream gage location (RWC: "11460151")
if(provide_coords) lat = 37.9; lon = -122.59 

### Define time period for historical analysis ###
# for GridMET and stream gage; Daymet period starts one year after this period 
startY = 1979; startM = 01; startD = 01 
endY = 2023; endM = 12; endD = 31
if(delayStart){ cutoffYear = startY+11 }else{cutoffYear = startY} 


### Model names ###
# provide list of model names to generate plots of future streamflow with those models highlighted
individual_models = c('BNU-ESM.rcp45')


### Set path variables ###
if(!dir.exists(here('Data', SiteID_FileName))) {dir.create(here('Data', SiteID_FileName))}; dataPath <- here('Data', SiteID_FileName)
if(!dir.exists(here('Output', SiteID_FileName))) {dir.create(here('Output', SiteID_FileName))}
if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow'))) {dir.create(here('Output', SiteID_FileName, 'Streamflow'))}
if(!dir.exists(here('Output', SiteID_FileName, 'Streamflow', FolderName))) {dir.create(here('Output', SiteID_FileName, 'Streamflow', FolderName))}; outLocationPath = here('Output', SiteID_FileName, 'Streamflow', FolderName)


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



### Alternatively, read in previously optimized variables from results file
if(file.exists(paste0(outLocationPath, "/optim_results.rds"))){
  results <- readRDS(paste0(outLocationPath, "/optim_results.rds"))
  gw_add <- results$gw_add; vfm <- results$vfm; jrange <- results$jrange; hock <- results$hock; hockros <- results$hockros; dro <- results$dro; 
  mondro <- results$mondro; aspect <- results$aspect; slope <- results$slope; shade.coeff <- results$shade.coeff; SWC.Max <- results$SWC.Max; 
  jtemp <- results$jtemp; qa <- results$qa; qb <- results$qb; sa <- results$sa; sb <- results$sb; va <- results$va; vb <- results$vb
}


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
### Get variables ###

# Pull watershed shapefile from StreamStats database
if(!provide_coords){
  coords <- get_coords(SiteID_FileName, GageSiteID); lat <- coords$lat; lon <- coords$lon; aoi <- coords$geometry
}
region <- get_region(lat,lon)

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
# what is this used for?? make sure it's not and then delete
#elev = get_elev_daymet(lat, lon, aoi, startY, endY, SiteID_FileName)

# create start and end date objects of data collection. Daymet will start one year after the year listed here
startDate <- ymd(paste(startY, startM, startD)); endDate <-  ymd(paste(endY, endM, endD))



#######################################################################
### Scrape and clean USGS stream gage data ###

gage_data <- get_gage_data(GageSiteID, incompleteMonths, fillLeapDays, dataPath)
meas_flow_daily <- gage_data$meas_flow_daily; meas_flow_mon <- gage_data$meas_flow_mon


#######################################################################
### Scrape and clean meteorological data (GridMET or Daymet) ### 

if(GridMET) {
  DailyClimData <- get_gridmet_data(SiteID_FileName, startY, endY, lat, lon, aoi, dataPath,
                             tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope)
} else { 
  DailyClimData <- get_daymet_data(SiteID_FileName, startY, endY, lat, lon, aoi, dataPath)
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


### Call optimization routine to get optimal variables ###
if(optimization){
  source('Streamflow//Optimization.R')
}


### Run model ###
DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
             shade.coeff, jtemp,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon)
DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
MeasMod<- MeasModWB(DailyDrain, meas_flow_mon, cutoffYear)


# if(!optimization){
#   # run model
#   DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
#                shade.coeff, jtemp,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon)
#   DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
#   MeasMod<- MeasModWB(DailyDrain, meas_flow_mon, cutoffYear)
#   nseM = NSE(sim = MeasMod$Mod, obs = MeasMod$Meas)
#   results = data.frame(results, nseM =nseM, gw_add = gw_add, vfm =vfm,jrange =jrange,hock = hock, hockros =hockros,
#                        dro =dro, mondro = mondro, aspect = aspect, slope = slope,
#                        shade.coeff = shade.coeff, jtemp =jtemp, elpTimeM = NA)
#   
#   # store results
#   IHcoeffs <- tibble()
#   nseD <-  IHACRESFlow(c(qa=qa, qb=qb, sa=sa, sb=sb, va=va), q0, s0, v0, DailyWB, meas_flow_daily, cutoffYear)
#   results<- data.frame(results, IHcoeffs, elpTimeD=NA)
#   
#   # save results as RDS file
#   saveRDS(results, file = paste0(outLocationPath, "/non_optim_results.rds"))
# }



#######################################################################
#######################################################################
### ANALYSIS ###

### MODEL PERFORMANCE ON HISTORICAL FLOW ###
source('Streamflow//Model_Accuracy.R')


### HISTORICAL STREAMFLOW ANALYSIS ###
if(historical_analysis){
  source('Streamflow//Historical_Analysis.R')
}


### FUTURE STREAMFLOW PROJECTIONS ###
if(future_analysis){
  source("Streamflow//Future_Analysis.R")
}

