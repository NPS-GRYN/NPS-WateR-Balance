# ---------------------------------------------------------------------
# This script includes functions related to the calculation of Water Balance model parameters. 
# These functions are either reproduced directly from the WaterBalance R package, 
# updated versions of those included in the WaterBalance R package or newly created functions
# that expand the capabilities of the WaterBalance R package.
# The original WaterBalance R package can be found here: https://github.com/CCRP-Adaptation/WaterBalance. 
# For many of the functions, units are not prescribed (but must match with the other inputs). When units 
# are prescribed, the function requires them.
#
# ---------------------------------------------------------------------

#######################################################################
### Functions to calculate water balance parameters ###

# Rain: Calculates rainfall totals based on precipitation and freeze factor.
# Args:
#   ppt: Time series vector of precipitation
#   freeze: Time series vector of freeze factors (0-1), calculated from average temperature and Jennings et al., 2018
# Returns:
#   Vector of rainfall totals
get_rain = function(ppt, freeze){
  rain = ppt*freeze
  return(rain)
}

# Snow: Calculates snowfall totals based on precipitation and freeze factor.
# Args:
#   ppt: Time series vector of precipitation
#   freeze: Time series vector of freeze factors (0-1), calculated from average temperature and Jennings et al., 2018
# Returns
#   Time series vector of snowfall totals
get_snow = function(ppt, freeze){
  snow = (1 - freeze)*ppt
  return(snow)
}

# Modify PET: Modifies PET by heat load according to method by Lutz et al. (2010).
# Args:
#   pet: Time series vector of PET (mm)
#   slope: Slope of the site (degrees)
#   aspect: Aspect of the site (degrees)
#   lat: Latitude of the site (degrees)
#   freeze: Time series vector of freeze factors (0-1)
#   shade.coeff: (optional) Shade coefficient from 0-1. Default is 1
# Returns:
#   Time series vector of modified PET (mm)
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
#   rain: Time series vector of daily rain
#   melt: Time series vector of daily snowmelt
# Returns:
#   Time series vector of water reaching soil surface
get_w = function(rain, melt){
  w = (melt+rain)
  return(w)
}

# Water reaching soil surface: Calculates water reaching soil surface, minus PET.
# Args:
#   w: Time series vector of water reaching soil surface as snow plus rain
#   pet: Time series vector of PET
# Returns:
#   Time series vector of water reaching soil surface, minus PET
get_w_pet = function(w, pet){
  w_pet = (w-pet)
  return(w_pet)
}


# Soil Water Content (SWC): Calculates soil water content from available water (rain + snowmelt), PET, max water-holding capacity, and initial SWC.
# Args: 
#   w: Time series vector of available water for soil charging (rain + snowmelt)
#   pet: Time series vector of PET
#   w_pet: Time series vector of the difference between w and pet
#   swc.max: Maximum soil water-holding capacity of soil layer
#   swc.0: (optional) The initial soil water content value. Default is 0
# Returns:
#   Time series vector of soil water content
get_soil = function(w, swc.0=NULL, pet, w_pet, swc.max){
  swc.i = ifelse(!is.null(swc.0), swc.0,0)
  soil=c()
  for(i in 1:length(pet)){
    soil[i] = ifelse(w[i]>pet[i], min((w_pet[i]+swc.i), swc.max), swc.i-swc.i*(1-exp(-(pet[i]-w[i])/swc.max)))
    swc.i=soil[i]
  }
  return(soil)  
}


# Daily change in Soil Water Content (SWC): Calculates daily change in soil water content.
# Args:
#   swc: Time series vector of soil water content
#   swc.0: (optional) Initial soil water content value. Default is 0
# Returns:
#   Time series vector of change in soil water content
get_d_soil=function(swc, swc.0=NULL){
  swc.0 = ifelse(!is.null(swc.0), swc.0, 0)
  d_soil = swc - lag(swc, default=swc.0)
  return(d_soil)
}

# Actual Evapotranspiration (AET): Calculates actual evapotranspiration (AET) from available water, PET, and soil water.
# Update: cap AET so it cannot be above PET, replacing ifelse(w[i] > pet[i], pet[i], w[i]+swc.i-swc[i])
# Args:
#   w:  Time series vector of available water for soil charging (rain + snowmelt)
#   pet: Time series vector of PET
#   swc: Time series vector of soil water content
#   swc.0: (optional) The initial soil water content value. Default is 0
# Returns:
#   Time series vector of AET
get_AET = function(w, pet, swc, swc.0=NULL){
  swc.i = ifelse(!is.null(swc.0), swc.0, 0)
  AET = numeric(length(w))
  for(i in 1:length(AET)){
    AET[i] = min(pet[i], w[i]+swc.i-swc[i])
    swc.i = swc[i]
  }
  return(AET)
}

# Climatic water deficit: Calculates daily climatic water deficit, which is PET - AET.
# Args:
#   pet: Time series vector of potential evapotranspiration
#   AET: Time series vector of actual evapotranspiration
# Returns:
#   Time series vector of climatic water deficit
get_deficit=function(pet, AET){
  deficit = pet-AET
  return(deficit)
}

# Growing Degree-Days: Calculates growing degree-days at daily time steps based on mean temperature and a threshold temperature.
# Args:
#   tmean Time series vector of daily mean temperatures (deg C)
#   tbase (optional) A threshold temperature, above which growing degree-days are calculated. Default is 0 (deg C)
# Returns:
#   Time series vector of growing degree days (deg C)
get_GDD = function(tmean, tbase=NULL){
  tb = ifelse(!is.null(tbase), tbase, 0)
  GDD = ifelse(tmean < tb, 0, tmean - tb)
  return(GDD)
}


# Calculate direct runoff 
# Update: WB package only had get runoff and not direct runoff; Excel had both
# Args:
#   DailyWB: A dataframe of daily water balance calculations
#   mondro: Direct runoff coefficient for monsoon season (July-October)
#   dro: Direct runoff coefficient for non-monsoon season (November-June)
#   tmean: Time series vector of daily mean temperatures (deg C)
#   jtemp: Jennings temperature (deg C)
#   jrange: Range of uncertainty for Jennings temperature (deg C)
#   precip: Time series vector of precipitation
#   month: Month of year (numerical)
# Returns: 
#   raindro: Time series vector of direct runoff from rain.
get_dro = function(DailyWB, mondro, dro, tmean, jtemp, jrange, precip, month){
  high = jtemp+jrange
  ranks <- percent_rank(precip)
  
  raindro <- numeric(nrow(DailyWB))
  raindro[(tmean>high & precip>0) & (month>6 & month<11)] <- ranks[(tmean>high & precip>0) & (month>6 & month<11)] * mondro * precip[(tmean>high & precip>0) & (month>6 & month<11)]
  raindro[(tmean>high & precip>0) & !(month>6 & month<11)] <- ranks[(tmean>high & precip>0) & !(month>6 & month<11)] * dro * precip[(tmean>high & precip>0) & !(month>6 & month<11)]

  return(raindro)
}



# Calculate effective precip (amount available after direct runoff)
# Update: Now that R code calculates direct runoff, it must be subtracted from precip 
# Args:
#   precip: Time series vector of precipitation
#   raindro: Time series vector of direct runoff from rain
# Returns:
#   Time series vector of precipitation available for infiltration, snow accumulation, or ET
get_precip = function(precip, raindro){
  upprec = precip - raindro
  return(upprec)
}


# Freeze factor using Jennings et al., 2018 thresholds to partition rain and snow
# Update: jrange was formerly hard coded to 3, update allows variation. 
# Args:
#   jtemp: Jennings temperature (deg C)
#   tmean: Time series vector of daily mean temperatures (deg C)
#   jrange: Range of uncertainty for Jennings temperature (deg C)
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
#   rain: Time series vector of rainfall
#   sp.0: (optional) Initial snowpack value. Default is 0
#   hockros: A melt factor of daily snowmelt when warm enought to melt and it is raining
#   hock: A melt factor of daily snowmelt when warm enough to melt
#   tmean: Time series vector of daily mean temperatures (deg C)
#   jtemp: Jennings temperature (deg C)
#   snow: Time series vector of snowfall values
#   jrange: Range of uncertainty for Jennings temperature (deg C)
# Returns:
#   Time series vector of snowpack melt
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


# Snowpack: Calculates snowpack accumulation at time steps from snowfall and melt.
# Args:
#   snow: Time series vector of snowfall
#   melt:  Time series vector of snowmelt
#   sp.0: optional) Initial snowpack value. Default is 0
# Returns:
#   Time series vector of snowpack accumulation
get_snowpack = function (snow, melt, sp.0 = NULL) 
{
  sp.i = ifelse(!is.null(sp.0), sp.0, 0)
  snowpack <- vector()
  for (i in 1:length(melt)) {
    snowpack[i] = min(sp.i + snow[i] - melt[i], 0)
    sp.i = snowpack[i]
  }
  return(snowpack)
}


# Runoff: Calculates runoff (excess input greater than soil water holding capacity) from water reaching soil surface, AET, soil moisture, and a runoff coefficient
# Args:
#   w: Time series vector of available water for soil charging (rain + snowmelt)
#   d_soil: Time series vector of change in soil moisture from previous day
#   RainDRO: Time series vector of direct runoff from rain
# Returns:
#   Time series vector of runoff
get_runoff = function (w, d_soil, AET, RainDRO) 
{
  runoff = w - d_soil - AET + RainDRO
  return(runoff)
}


# Calculate adjusted runoff
# Update: use gw_add and vfm to adjust runoff
# Args:
#   orig: Time series vector containing initial runoff values
#   gw_add: amount of water contributed to the streamflow by groundwater, also known as baseflow (intercept)
#   vfm: volume forcing multiplier (1 corresponds to no change in volume) (slope)
# Returns:
#   Time series vector of runoff multiplied by vfm, with groundwater addition added
get_adj_runoff= function(orig, gw_add, vfm){
  adjusted = (orig*vfm)+gw_add
  return(adjusted)
}

# Calculate adjusted ET (PET or AET) 
# Args:
#   et: Time series vector of ET (PET or AET)
#   et_slope: multiplication component of bias correction
#   et_bias: addition component of bias correction
# Returns:
#   Time series vector of bias-corrected ET data, multiplied by slope adjustment and with bias added
get_adj_et = function(aet, pet, et_bias, et_slope){
  adjusted = pmin(aet * et_slope + et_bias, pet)
  return(adjusted)
}

# Temperature threshold using Jennings et al., 2018 to partition rain and snow: Extracts rain-snow temperature threshold from a raster.
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

#######################################################################
### Supporting functions for ET calculation ###

# Daylength: Returns daylength in hours for a series of dates, based on latitude. Calls the 'geosphere' package.
# Args:
#   dates: A series of dates containing year, month, and day
#   lat: Latitude (degrees)
# Returns:
#   Daylength in hours for series of dates
get_daylength = function(dates, lat){
  yday = as.numeric(strftime(dates, "%j"))
  dayl_h = geosphere::daylength(lat, yday)
  return(dayl_h)
}

# Saturation Vapor Pressure: Calculates mean saturation vapor pressure of air based on temperature.
# Args:
#   temp: Temperature (deg C)
# Returns:
#   Saturation vapor pressure (kPa)
get_svp = function(temp){
  svp = 0.6108*exp((17.27*temp)/(temp + 237.3))
  return(svp)
}

# Relative Humidity: Calculates relative humidity from atmospheric vapor pressure and temperature
# Args:
#   vp: Vapor pressure (kPa)
#   temp: Temperature (deg C)
# Returns:
#   Relative humidity (%)
get_rh = function(vp, temp){
  svp = get_svp(temp)
  rh = vp/svp
  return(rh)
}

# Actual Vapor Pressure: Calculates actual vapor pressure of air based on maximum and minimum relative humidity and maximum and minimum temperature.
# Args:
#   rhmax: Daily maximum relative humidity (%).
#   rhmin: Daily minimum relative humidity (%).
#   tmax: Daily maximum temperature (deg C).
#   tmin: Daily minimum temperature (deg C).
# Returns:
#   Actual vapor pressure (kPa)
actual_vp = function(rhmax, rhmin, tmax, tmin){
  e.tmax = get_svp(tmax)
  e.tmin = get_svp(tmin)
  e.a = (e.tmin*(rhmax/100) + e.tmax*(rhmin/100))/2
  return(e.a)
}

# Slope of Saturation Vapor Curve: Calculates the slope of the saturation vapor curve for a given temperature.
# Args:
#   temp: A time series vector or single value of temperatures (deg C).
# Returns
#   Slope of saturation vapor curve
vapor_curve = function(temp){
  vap.curve = 4098*(0.6108*exp((17.27*temp/(temp+237.3)))/(temp+273.3)^2)
  return(vap.curve)
}

# Atmospheric Pressure: Estimates atmospheric pressure (kPa) at a given elevation.
# Args:
#   elev: Elevation (m).
# Returns:
#   Atmospheric pressure (kPa)
atm_press = function(elev){
  atm.press = 101.3*((293 - 0.0065*elev)/293)^5.26
  return(atm.press)
}

# Psychrometric Constant: Calculates the psychrometric constant relating partial pressure of water in air to the air temperature, based on atmospheric pressure. Calls the atm_press() function to estimate atmospheric pressure from elevation.
# Args:
#   elev: Elevation (m)
# Returns:
#   Psychrometric constant
psyc_constant = function(elev){
  atm.press = atm_press(elev)
  psyc.const = 0.000665*atm.press
  return(psyc.const)
}

# Clear Sky Radiation: Calculates incoming clear-sky radiation based on day-of-year, latitude, and elevation
# Args:
#   doy: Day-of-year (Julian date)
#   lat: Latitude (degrees)
#   elev: Elevation (m)
# Returns:
#   Clear sky radiation (MJ m^-2 day^-1)
clear_sky_rad = function(doy, lat, elev){
  d.r = 1 + 0.033*cos(((2*pi)/365)*doy)
  declin = 0.409*sin((((2*pi)/365)*doy)-1.39)
  lat.rad = (pi/180)*lat
  sunset.ang = acos(-tan(lat.rad)*tan(declin))
  R.a = ((24*60)/pi)*0.0820*d.r*(sunset.ang*sin(lat.rad)*sin(declin) + cos(lat.rad)*cos(declin)*sin(sunset.ang))
  R.so = (0.75 + 2e-5*elev)*R.a
  return(R.so)
}

# Outgoing Radiation: Calculates outgoing radiation based on daily Tmax, Tmin, incoming radiation, actual vapor pressure, and clear-sky radiation.
# Args:
#   tmax: Daily maximum temperatures (deg C).
#   tmin Daily minimum temperatures (deg C).
#   R.s Incoming solar radiation (MJ m^-2 day^-1).
#   e.a Actual vapor pressure (kPa).
#   R.so Clear-sky radiation (MJ m^-2 day^-1).
# Returns:
#   Outgoing radiation (MJ m^-2 day^-1)
outgoing_rad = function(tmax, tmin, R.s, e.a, R.so){
  R.nl = 4.903e-09*(((tmax + 273.16)^4 + (tmin + 273.16))/2)*(0.34-0.14*sqrt(e.a))*(1.35*(R.s/R.so) - 0.35)
  return(R.nl)
}

#######################################################################
### ET Calculation Methods ###

# Hamon Daily PET: Calculates Hamon PET from a daily time series of Tmean and daylength.
# Args:
#   tmax: a daily resolution array containing max temperature (deg C); tmin: a daily resolution array containing min temperature (deg C); date: a daily resolution array containing dates; lat: latitude
# Returns:
#   Daily time series of PET calculated according to Hamon method
ET_Hamon_daily = function(tmax, tmin, date, lat){
  tmean = (tmax + tmin)/2
  daylength = get_daylength(date, lat)
  
  et.hamon = 0.1651*(daylength/12)*(216.7*(6.108*exp((17.26*tmean)/(tmean+273.3))))/(tmean+273.3)
  
  return(et.hamon)
}

# Thornthwaite Monthly PET: Calculates PET from monthly Tmean and daylength, according to the Thornthwaite method.
# Args:
#   x: A monthly time series data frame containing Date, tmean_C (deg C), and daylength (hours)
# Returns:
#   Monthly time series of PET accordint to Thornthwaite method
ET_Thorn_monthly = function(x){
  x$month = strftime(x$Date, "%m")
  N = lubridate::days_in_month(as.numeric(x$month))
  e.s = get_svp(x$tmean_C)
  et.thorn = ifelse(x$tmean_C > 0, 29.8*N*x$daylength*(e.s/(x$tmean_C+273.2)), 0)
  return(et.thorn)
}

# Penman-Monteith Daily PET: Calculates PET (mm) from daily Tmax, Tmin, solar radiation, elevation, and latitude, according to the Penman-Monteith method. May also use daily maximum and minimum relative humidity, atmospheric vapor pressure, and wind speeds.
# EDIT INPUTS
# Args:
#   x: A daily time series data frame containing Date (date object), tmax_C (deg C), tmin_C (deg C), srad (MJ m^-2 day^-1). Optionally contains RHmax (percent), RHmin (percent), vp (kPa), and wind (m/s).
#   elev: Elevation of the site (m).
#   lat: Latitude of the site (degrees).
#   wind: (optional) An estimated value for daily average wind speeds (m/s). Use if input data frame does not contain daily wind speed values.
# Returns:
#   Daily time series of PET according to Penman-Monteith Method
ET_PenmanMonteith_daily = function(date, tmax, tmin, srad, vpd, vs, elev, lat){
  # Calculate inputs
  tmean = (tmax + tmin)/2
  doy = as.numeric(strftime(date, "%j"))
  #rh.max = x$RHmax
  #rh.min = x$RHmin
  R.s = srad * 0.0864 # convert from W/m2 to MJ/m2
  psyc.const = psyc_constant(elev)
  vap.curve = vapor_curve(tmean)
  
  #Auxilary calculations for wind terms
  DT = vap.curve/(vap.curve + psyc.const*(1+0.34*vs))
  PT = psyc.const/(vap.curve + (psyc.const*(1+0.34*vs)))
  TT = (900/(tmean + 273))*vs
  
  #Saturation vapor pressure
  e.tmax = get_svp(tmax)
  e.tmin = get_svp(tmin)
  e.s = (e.tmax + e.tmin)/2
  
  #Actual vapor pressure
  if(is.null(vpd) == TRUE){
    if(is.null(rh.max) == TRUE){
      e.a = e.tmin
    } else {
      e.a = actual_vp(rh.max, rh.min)
    }
  } else {
    e.a = vpd
  }
  
  #Solar angle and radiation calculations
  R.ns = (1 - 0.23)*R.s
  R.so = clear_sky_rad(doy, lat, elev)
  R.nl = outgoing_rad(tmax, tmin, R.s, e.a, R.so)
  R.n = R.ns - R.nl
  R.ng = 0.408*R.n
  
  #ET from radiation
  ET.rad = DT*R.ng
  #ET from wind
  ET.wind = PT*TT*(e.s - e.a)
  #Total ET
  ET.o = ET.rad + ET.wind
  return(ET.o)
}

# Oudin Daily PET: Calculates PET (mm) based on temperature, latitude, and solar radiation 
# Updates: removed heatload and terrain correction from within this function
# Args:
#   doy: Day-of-year (Julian date)
#   lat: Latitude of the site (degrees).
#   snowpack: A time series vector of snowpack accumulation values.
#   tmean: A vector of daily mean temperatures (deg C).
#   slope: Slope of the site (in degrees).
#   aspect: Aspect of the site (in degrees).
#   shade.coeff: (optional) A shade coefficient from 0-1. Default is 1.
# Returns:
#   Daily time series of PET according to Oudin method
get_OudinPET = function(doy, lat, snowpack, tmean, slope, aspect, shade.coeff=NULL){
  d.r = 1 + 0.033*cos((2*pi/365)*doy)
  declin = 0.409*sin((((2*pi)/365)*doy)-1.39)
  lat.rad = (pi/180)*lat
  sunset.ang = acos(-tan(lat.rad)*tan(declin))
  R.a = ((24*60)/pi)*0.082*d.r*((sunset.ang*sin(lat.rad)*sin(declin)) + (cos(lat.rad)*cos(declin)*sin(sunset.ang)))
  Oudin = ifelse(snowpack>2, 0, ifelse(tmean>-5, (R.a*(tmean+5)*0.408)/100, 0))
  #Folded_aspect = abs(180-abs((aspect)-225))
  #Heatload = (0.339+0.808*cos(lat*(pi/180))*cos(slope*(pi/180)))-(0.196*sin(lat.rad)*sin(slope*(pi/180)))-(0.482*cos(Folded_aspect*(pi/180))*sin(slope*(pi/180)))
  #sc = ifelse(!is.null(shade.coeff), shade.coeff, 1)
  #OudinPET = Oudin * Heatload * sc
  return(Oudin)
}


#######################################################################
### Supporting Functions ###

# Adjust precipitation and temperature data
# Args:
#   orig: Time series vector with original meteorological data (i.e. temperature or precip)
#   bias: addition component of bias correction (intercept)
#   slopeadj: multiplication component of bias correction (slope)
# Returns:
#   Time series vector of bias-corrected meteorological data, multiplied by slope adjustment and with bias added
get_slope_bias_adj= function(orig, bias, slopeadj){
  adjusted = orig*slopeadj+bias
  return(adjusted)
}

# Convert Kelvin to Celsius
# Args:
#   kelvin: vector or value of temperature (K)
# Returns:
#   vector or value of temperature (deg C)
kelvin_to_celcius <- function(kelvin) {
  return(kelvin- 273.15)
}
