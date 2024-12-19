####these are functions that did exist in the WB package and I edited or that didn't exist and I created

#I edited this function to match the calculation in Excel
#The WB package only had get runoff and not direct runoff
get_dro = function(DailyWB, mondro = mondro, dro = dro, tmean, jtemp, jrange, precip, month){
  high = jtemp+jrange
  raindro<- c()
  for(i in 1:nrow(DailyWB)){
    raindro[i]<- if (tmean[i]>high&precip[i]>0) {
      if (month[i]>6&month[i]<11) {
        percent_rank(precip)[i]*mondro*precip[i]
      }else(percent_rank(precip)[i]*dro*precip[i])
    }else(0)
  }
  return(raindro)
}

#I made this to match Excel because the WB package only had runoff and not direct runoff
#This meant they did not subtract the direct runoff from precip
get_precip = function(precip, raindro){
  upprec = precip- raindro
  return(upprec)
}

#I edited this function to allow for jrange, It was hard coded to 3
get_freeze = function (jtemp, tmean, jrange=jrange) 
{
  freeze = ifelse(tmean <= (jtemp - jrange), 0, ifelse(tmean >= 
  (jtemp + jrange), 1, (1/((jtemp + jrange) - (jtemp - jrange))) * (tmean - (jtemp - jrange))))
}

#I edited this to include Final Hock which did not Exist in R, it just had one hock value
# I also added jrange to this function
#This function has two parts; the first sets the first values. 
#The for loop iterates starting at 2 because the first one has been set

get_melt = function (rain, sp.0, hockros, hock, tmean, jtemp, snow, jrange = jrange) 
{
  sp.0 = ifelse(!is.null(sp.0), sp.0, 0)
  finhock<- c()
  finhock[1]<- if (rain[1]>0&sp.0>0) {
    hockros
  } else if (rain[1]==0 &sp.0>0) {
    hock
  } else {
    0
  }
  low = jtemp-jrange
  melt <- vector()
  melt[1] = ifelse(tmean[1] < low | sp.0 == 0, 0, 
                   ifelse((tmean[1] - low) * finhock[1] > sp.0, sp.0, (tmean[1] - low) * finhock[1]))
  snowpack <- vector()
  snowpack[1] = sp.0 + snow[1] - melt[1]
  
  
  for (i in 2:length(tmean)) {
    finhock[i]<- if (rain[i]>0&snowpack[i-1]>0) {
      hockros
    } else if (rain[i]==0 &snowpack[i-1]>0) {
      hock
    } else {
      0
    }
    melt[i] = ifelse(tmean[i] < low| snowpack[i - 1] == 0, 0,
                     ifelse((tmean[i] - low) * finhock[i] > snowpack[i - 1], snowpack[i - 1], (tmean[i] - low) * finhock[i]))
    snowpack[i] = snowpack[i - 1] + snow[i] - melt[i]
  }
  return(melt)
}

#I removed a line of unnecessary code which defined low_thresh_temp and never called it again
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

#I edited this to match the calculations in Excel
get_runoff = function (w, d_soil, AET, RainDRO) 
{
  runoff = w - d_soil - AET + RainDRO
  return(runoff)
}

#I added this function for adjusting the runoff
get_adj_runoff= function(orig, gw_add, vfm){
  adjusted = (orig*vfm)+gw_add
  return(adjusted)
}

#I added this function for adjusting the Gridmet Precip and temperature data
get_slope_bias_adj= function(orig, bias, slopeadj){
  adjusted = orig*slopeadj+bias
  return(adjusted)
}

# I made this very complex function to convert Kelvin to Celcius
kelvin_to_celcius <- function(kelvin) {
  return(kelvin- 273.15)
}

#I edited this fixing the incorrect order of lat lon in coords and moved the call of the raster to outside of the function
get_jtemp = function(lat, lon, j.raster){
  projection = sp::CRS("+init=epsg:4326")
  coords = cbind(lon, lat)
  sp = sp::SpatialPoints(coords, proj4string = projection)
  jtemp = raster::extract(j.raster, sp)
  return(jtemp)
}

#I created this function. It scrapes the DayMet data and extracts the elevation.
#Note, if the format of the DayMet data changes, this function may not perform correctly
#It relies on the elevation being in a specific row and column and having specific text around it
get_elev_DayMet = function(lat=lat, lon=lon, startY =startY, endY = endY, ClimateSiteID_FileName){
  if(!file.exists(file.path(dataPath, paste0(paste("DayMet", ClimateSiteID_FileName, startY+1,endY, sep = "_"), ".csv")))){
    sites<- data.frame(site = paste("DayMet", ClimateSiteID_FileName, sep = "_"), latitude = lat,longitude = lon)
    write.csv(sites, file =file.path(dataPath,paste0("SiteFile", ClimateSiteID_FileName,".csv")), row.names = F)
    daymetr::download_daymet_batch(file_location = file.path(dataPath,paste0("SiteFile", ClimateSiteID_FileName, ".csv")),
                          start = startY+1,
                          end = endY,
                          internal = FALSE,
                          path=dataPath)
  }
  daymet<- read.csv(file.path(dataPath, paste0(paste("DayMet", ClimateSiteID_FileName, startY+1,endY, sep = "_"), ".csv")), skip = 0, header = TRUE, sep = ",")
  chars<- "Elevation: "
  chars2<- " meters"
  elev<- as.numeric(gsub(chars2,"" , gsub(chars, '', daymet[3,1])))
  return(elev)
}