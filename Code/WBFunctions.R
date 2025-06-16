# ---------------------------------------------------------------------
# This script includes functions related to the running of the Water Balance model. 
# These functions are either reproduced directly from the WaterBalance R package, 
# updated versions of those included in the WaterBalance R package or newly created functions
# that expand the capabilities of the WaterBalance R package.
# The original WaterBalance R package can be found here: https://github.com/CCRP-Adaptation/WaterBalance. 
#
# EDITS IN PROGRESS
# improve documentation
# ---------------------------------------------------------------------

### From WaterBalance package  ###

# Rain: Calculates rainfall totals based on precipitation and freeze factor.
# Args:
#   ppt: A vector of precipitation values.
#   freeze: A vector of freeze factor values, calculated from average temperature and Jennings et al., 2018. Values are 0-1.
# Returns:
#   Vector of rainfall totals
get_rain = function(ppt, freeze){
  rain = ppt*freeze
  return(rain)
}

# Snow: Calculates snowfall totals based on precipitation and freeze factor.
# Args:
#   ppt: A vector of precipitation values.
#   freeze: A vector of freeze factor values, calculated from average temperature and Jennings et al., 2018. Values are 0-1.
# Returns
#   Vector of snowfall totals
get_snow = function(ppt, freeze){
  snow = (1 - freeze)*ppt
  return(snow)
}

# Modify PET: Modifies PET by heat load according to method by Lutz et al. (2010).
# Args:
#   pet: A time series vector of PET values.
#   slope: Slope of the site (in degrees).
#   aspect: Aspect of the site (in degrees).
#   lat: Latitude of the site (in degrees).
#   freeze: A vector of freeze factor values, calculated from Tmean. Values are 0-1.
#   shade.coeff: (optional) A shade coefficient from 0-1. Default is 1.
# Returns:
#   Time series vector of modified PET values
modify_PET = function(pet, slope, aspect, lat, freeze, shade.coeff=NULL){
  f.aspect = abs(180 - abs(aspect - 225))
  lat.rad = ifelse(lat > 66.7, (66.7/180)*pi, (lat/180)*pi)
  slope.rad = (slope/180)*pi
  aspect.rad = (f.aspect/180)*pi
  heat.load = 0.339+0.808*cos(lat.rad)*cos(slope.rad) - 0.196*sin(lat.rad)*sin(slope.rad) - 0.482*cos(aspect.rad)*sin(slope.rad)
  sc = ifelse(!is.null(shade.coeff), shade.coeff, 1)
  freeze = ifelse(freeze == 0,0,1)
  PET.Lutz = pet*heat.load*sc*freeze
  return(PET.Lutz)
}

# Water reaching soil surface as snow plus rain: Calculates water reaching soil surface using rain and snowmelt.
# Args:
#   rain: A vector of daily rain.
#   melt: A vector of daily snowmelt.
# Returns:
#   Vector of water reaching soil surface
get_w = function(rain, melt){
  w = (melt+rain)
  return(w)
}

# Water reaching soil surface minus PET: Calculates water reaching soil surface minues the PET.
# Args:
#   w: A time series vector of water reaching soil surface as snow plus rain.
#   pet: A time series vector of PET values.
# Returns:
#   Time series vector of water reaching soil surface, minus PET
get_w_pet = function(w, pet){
  w_pet = (w-pet)
  return(w_pet)
}


# Soil Water Content (SWC): Calculates soil water content from available water (rain + snowmelt), PET, max. water-holding capacity, and initial SWC.
# Args: 
#   w: A time series vector of available water for soil charging (rain + snowmelt).
#   pet: A time series vector of PET.
#   w_pet: A time series vector of the difference between w and pet.
#   swc.max: The maximum soil water-holding capacity of the soil layer being assessed.
#   swc.0: (optional) The initial soil water content value. Default is 0.
# Returns:
#   Time series vector of soil water content
get_soil = function(w, swc.0=NULL, pet, w_pet, swc.max){
  swc.i = ifelse(!is.null(swc.0), swc.0,0)
  soil=c()
  for(i in 1:length(pet)){
    soil[i] = ifelse(w[i]>pet[i], min((w_pet[i]+swc.i),swc.max), swc.i-swc.i*(1-exp(-(pet[i]-w[i])/swc.max)))
    swc.i=soil[i]
  }
  return(soil)  
}


# Daily change in Soil Water Content (SWC): Calculates daily change in soil water content.
# Args:
#   swc: A time series vector of soil water content.
#   swc.0: (optional) The initial soil water content value. Default is 0.
# Returns:
#   Time series vector of delta soil water content
get_d_soil=function(swc, swc.0=NULL){
  swc.0 = ifelse(!is.null(swc.0), swc.0, 0)
  d_soil = swc - lag(swc, default=swc.0)
  return(d_soil)
}

# Actual Evapotranspiration (AET): Calculates actual evapotranspiration (AET) from available water, PET, and soil water.
# Args:
#   w:  A time series vector of available water for soil charging (rain + snowmelt).
#   pet: A time series vector of PET.
#   swc: A time series vector of soil water content.
#   swc.0: (optional) The initial soil water content value. Default is 0.
# Returns:
#   Time series vector of AET
get_AET = function(w, pet, swc, swc.0=NULL){
  swc.i = ifelse(!is.null(swc.0), swc.0, 0)
  AET = numeric(length(w))
  for(i in 1:length(AET)){
    AET[i] = ifelse(w[i] > pet[i], pet[i], w[i]+swc.i-swc[i])
    swc.i = swc[i]
  }
  return(AET)
}

# Climatic water deficit: Calculates daily climatic water deficit, which is PET - AET.
# Args:
#   pet: A time series vector of PET.
#   AET: A time series vector of actual evapotranspiration.
# Returns:
#   Time series vector of climatic water deficit
get_deficit=function(pet, AET){
  deficit = pet-AET
  return(deficit)
}

# Growing Degree-Days: Calculates growing degree-days at daily time steps based on mean temperature and a threshold temperature.
# Args:
#   tmean A time series vector of daily mean temperatures (deg C).
#   tbase (optional) A threshold temperature, above which growing degree-days are calculated. Default is 0.
# Returns:
#   Time series vector of growing degree days
get_GDD = function(tmean, tbase=NULL){
  tb = ifelse(!is.null(tbase), tbase, 0)
  GDD = ifelse(tmean < tb, 0, tmean - tb)
  return(GDD)
}


### Modified or new water balance functions ###
# Joseph Crane added these functions to improve model accuracy and match the R model to the Excel model. 

# Calculate direct runoff 
# Update: WB package only had get runoff and not direct runoff; Excel had both
# Improve effiency (removed for loop)
# Args:
#   DailyWB: A dataframe of 
#   mondro:
#   dro:
#   tmean: A vector of daily mean temperatures (deg C).
#   jtemp: the Jennings temperature extracted from the raster based on latitude and longitude.
#   jrange:
#   precip:
#   month:
# Returns
#   raindro: 
get_dro = function(DailyWB, mondro, dro, tmean, jtemp, jrange, precip, month){
  high = jtemp+jrange
  ranks <- percent_rank(precip)
  
  raindro <- numeric(nrow(DailyWB))
  raindro[(tmean>high & precip>0) & (month>6 & month<11)] <- ranks[(tmean>high & precip>0) & (month>6 & month<11)] * mondro * precip[(tmean>high & precip>0) & (month>6 & month<11)]
  raindro[(tmean>high & precip>0) & !(month>6 & month<11)] <- ranks[(tmean>high & precip>0) & !(month>6 & month<11)] * dro * precip[(tmean>high & precip>0) & !(month>6 & month<11)]

  return(raindro)
}



# Calculate precip
# Update: Now that R code calculates direct runoff, it must be subtracted from precip 
# Args:
# Returns:
get_precip = function(precip, raindro){
  upprec = precip - raindro
  return(upprec)
}

# Freeze factor using Jennings et al., 2018 thresholds to partition rain and snow
# Update: jrange was formerly hard coded to 3, update allows variation. 
# Args:
#   jtemp: the Jennings temperature extracted from the raster based on latitude and longitude.
#   tmean: A vector of daily mean temperatures (deg C).
#   jrange:
# Returns:
#   Freeze factor from 0-1 based on a temperature threshold from Jennings et al., 2018 and average temperature
get_freeze = function (jtemp, tmean, jrange) 
{
  freeze = ifelse(tmean <= (jtemp - jrange), 0, ifelse(tmean >= 
  (jtemp + jrange), 1, (1/((jtemp + jrange) - (jtemp - jrange))) * (tmean - (jtemp - jrange))))
}

# Melt: Calculates the amount of snowmelt at time steps from snowpack, temperature, and Hock melt factor.
# Update: Includes updated final hock value (combination of 2) and adds jrange
# Args:
#   rain: A vector of rainfall totals
#   sp.0: (optional) Initial snowpack value. Default is 0.
#   hockros: A melt factor of daily snowmelt during rain on snow event (?)
#   hock: A melt factor of daily snowmelt when warm enough to melt.
#   tmean: A vector of daily mean temperatures (deg C).
#   jtemp: the Jennings temperature extracted from the raster based on latitude and longitude.
#   snow: A time series vector of snowfall values.
#   jrange:
# Returns:
get_melt = function (rain, sp.0, hockros, hock, tmean, jtemp, snow, jrange) {
  # Set initial values
  sp.0 = ifelse(!is.null(sp.0), sp.0, 0)
  low = jtemp-jrange
  
  # Initialize vectors
  finhock <- numeric(length(tmean))
  melt <- numeric(length(tmean))
  snowpack <- numeric(length(tmean))
    
  # Establish values for first day in time series
  finhock[1]<- if(rain[1]>0&sp.0>0) hockros else if (rain[1]==0 &sp.0>0) hock else 0
  melt[1] = ifelse(tmean[1] < low | sp.0 == 0, 0, 
                   ifelse((tmean[1] - low) * finhock[1] > sp.0, sp.0, (tmean[1] - low) * finhock[1]))
  snowpack[1] = sp.0 + snow[1] - melt[1]
  
  # Iterate through remainder of the time series
  for (i in 2:length(tmean)) {
    finhock[i]<- if (rain[i]>0&snowpack[i-1]>0) hockros else if (rain[i]==0 &snowpack[i-1]>0) hock else 0
    melt[i] = ifelse(tmean[i] < low| snowpack[i - 1] == 0, 0,
                     ifelse((tmean[i] - low) * finhock[i] > snowpack[i - 1], snowpack[i - 1], (tmean[i] - low) * finhock[i]))
    snowpack[i] = snowpack[i - 1] + snow[i] - melt[i]
  }
  return(melt)
}

# Snowpack: Calculates snowpack accumulation at time steps, from a time series of snowfall and melt.
# Args:
#   snow: A time series vector of snowfall values.
#   melt:  A time series vector of snowmelt.
#   sp.0: optional) Initial snowpack value. Default is 0.
# Returns:
#   A time series vector of snowpack accumulation
get_snowpack = function (snow, melt, sp.0 = NULL) 
{
  sp.i = ifelse(!is.null(sp.0), sp.0, 0)
  snowpack <- vector()
  for (i in 1:length(melt)) {
    snowpack[i] = sp.i + snow[i] - melt[i]
    sp.i = snowpack[i]
  }
  return(snowpack)
}

# Runoff or excess input greater than soil water holding capacity: Calculates runoff at daily timesteps based on water reaching soil surface, AET, change in soil moisture, and a runoff coefficient
# Args:
#   w: A time series vector of available water for soil charging (rain + snowmelt).
#   d_soil: A time series vector of change in soil moisture from previous day.
#   RainDRO:
# Returns:
#   A time series vector of runoff
get_runoff = function (w, d_soil, AET, RainDRO) 
{
  runoff = w - d_soil - AET + RainDRO
  return(runoff)
}

# Calculate adjusted runoff
# Update: use gw_add and vfm to adjust runoff
# Args:
#   orig: vector containing initial runoff values
#   gw_add: amount of water contributed to the streamflow by groundwater, also known as baseflow (mm)
#   vfm: volume forcing multiplier (1 corresponds to no change in volume)
# Returns:
#   Vector of runoff multiplied by vfm, with groundwater addition added
#I added this function for adjusting the runoff
get_adj_runoff= function(orig, gw_add, vfm){
  adjusted = (orig*vfm)+gw_add
  return(adjusted)
}


# Adjust precipitation and temperature data
# Update: function did not exist before
# Args:
#   orig: vector with original meteorological data (i.e. temperature or precip)
#   bias: addition component of bias correction
#   slopeadj: multiplication component of bias correction
# Returns:
#   Vector of bias-corrected meteorological data, multiplied by slope adjustment and with bias added
get_slope_bias_adj= function(orig, bias, slopeadj){
  adjusted = orig*slopeadj+bias
  return(adjusted)
}

# Convert Kelvin to Celsius
# Update: function did not exist before
# Args:
#   kelvin: vector or value of temperature in K
# Returns:
#   vector or value of temperature in C
kelvin_to_celcius <- function(kelvin) {
  return(kelvin- 273.15)
}

# Temperature threshold using Jennings et al., 2018 to partition rain and snow: Extracts the rain-snow temperature threshold from a raster.
# Update: moved call to raster outside function and fixed incorrect order of coordinates
# Args:
#   lat: Latitude of the site (degrees).
#   lon: Longitude of the site (degrees).
#   j.raster: location of merged_jennings.tif file containing Jennings coefficient data
# Returns:
#   Jennings coefficient for rain/snow partitioning
get_jtemp = function(lat, lon, j.raster){
  projection = sp::CRS("+init=epsg:4326")
  coords = cbind(lon, lat)
  sp = sp::SpatialPoints(coords, proj4string = projection)
  jtemp = raster::extract(j.raster, sp)
  return(jtemp)
}

# Extract elevation for given location using DayMet dataset
# Update: function did not exist before 
# Note, if the format of the DayMet data changes, this function may not perform correctly
# It relies on the elevation being in a specific row and column and having specific text around it
# not actually sure what the purpose of this function is
# Args:
#   lat: Latitude of the site (degrees).
#   lon: Longitude of the site (degrees).
#   aoi: Shapefile of the basin
#   startY: First year to get data for (cannot be earlier than 1980)
#   endY: Last year to get data for (cannot be later than the most recent completed calendar year)
#   SiteID_FileName:
# Returns:
#   Elevation of point, in m
get_elev_daymet = function(lat, lon, aoi, startY, endY, SiteID_FileName){
  if(!file.exists(file.path(dataPath, paste0(paste("DayMet", SiteID_FileName, startY+1,endY, sep = "_"), ".csv")))){
    get_daymet_data(SiteID_FileName, startY, endY, lat, lon, aoi, dataPath)
  } 
  daymet <- read.csv(file.path(dataPath, paste0(paste("DayMet", SiteID_FileName, startY+1,endY, sep = "_"), ".csv")), skip = 0, header = TRUE, sep = ",")
  chars <- "Elevation: "
  chars2 <- " meters"
  elev <- as.numeric(gsub(chars2,"" , gsub(chars, '', daymet[3,1])))
  return(elev)
}
