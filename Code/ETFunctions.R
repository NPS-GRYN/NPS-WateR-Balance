# ---------------------------------------------------------------------
# This script includes functions to calculate potential evapotranspiration from input climate data,
# in support of the water balance model. Functions are based on Dave Thoma's water balance Excel 
# spreadsheet model.
# 
# EDITS IN PROGRESS
# test all ET functions (except for Oudin, which is verified to be working currently)
# improve documentation
# ---------------------------------------------------------------------


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

################################ ET Calculation Methods ##########################################

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
  Oudin = ifelse(snowpack>2,0,ifelse(tmean>-5,(R.a*(tmean+5)*0.408)/100,0))
  Folded_aspect = abs(180-abs((aspect)-225))
  Heatload = (0.339+0.808*cos(lat*(pi/180))*cos(slope*(pi/180)))-(0.196*sin(lat.rad)*sin(slope*(pi/180)))-(0.482*cos(Folded_aspect*(pi/180))*sin(slope*(pi/180)))
  sc = ifelse(!is.null(shade.coeff), shade.coeff, 1)
  OudinPET = Oudin * Heatload * sc
  return(OudinPET)
}