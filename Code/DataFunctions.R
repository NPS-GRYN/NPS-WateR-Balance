# ---------------------------------------------------------------------
# This script includes functions for pulling, cleaning, and managing data required to run
# and calibrate the water balance and IHACRES streamflow models.
#
# EDITS IN PROGRESS
# function for pulling Mike Tercek's water balance data (website currently appears to be down?)
# add functions for pulling data for aoi (not point): MACA, gridmet/daymet, openET
# once both of those are done, figure out how to pull Mike data and calculate watershed average (if that's feasible)
# finish documentation
# ---------------------------------------------------------------------


# Load libraries
if(!lib_install){
  library(sf); library(raster); library(ggplot2); library(dplyr); library(xts); library(geosphere); library(quantreg); library(orca)
  library(lubridate); library(hydroGOF); library(stringr); library(terra); library(glue); library(tidyverse); library(RColorBrewer)
  library(climateR); library(EGRET); library(daymetr); library(here); library(ggrepel); library(gridExtra); library(Kendall)
  library(httr); library(jsonlite); library(sf); library(grid); library(GA); library(GGally); library(data.table); library(plotly)
  library(tseries); library(dgof); library(wql)
  lib_install <- TRUE
}


# Create folders for data and output
if(!dir.exists(here('Data'))) {dir.create(here('Data'))}
if(!dir.exists(here('Output'))) {dir.create(here('Output'))}

# list of GCMs
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')


# Get latitude and longitude coordinates of watershed centroid using StreamStats database
# Args:
#   SiteID_FileName:
#   GageSiteID:
# Returns:
#   Latitude and longitude coordiantes of watershed centroid
get_coords <- function(SiteID_FileName, GageSiteID){
  if(!file.exists(here('Data', SiteID_FileName, 'downloaded_shapefile', 'Layers', 'globalwatershed.shp'))){
    # get lat/lon data for watershed
    INFO <- readNWISInfo(siteNumber=GageSiteID, parameterCd = "", interactive=FALSE)
    
    # produce workspace ID for given watershed
    get_workspace = GET(paste0('https://streamstats.usgs.gov/streamstatsservices/watershed.geojson?rcode=',tail(strsplit(INFO$station_nm," ")[[1]],n=1),'&xlocation=',INFO$dec_long_va,"&ylocation=",INFO$dec_lat_va,"&crs=4326&includeparameters=true&includeflowtypes=false&includefeatures=true&simplify=true"))
    workspaceID <- fromJSON(content(get_workspace, as = "text", encoding = "UTF-8"))$workspaceID
    
    # use workspace ID to download shapefile
    get_download <- GET(paste0("https://streamstats.usgs.gov/streamstatsservices/download?workspaceID=",workspaceID,"&format=SHAPE"))
    writeBin(content(get_download, "raw"), here("Data",SiteID_FileName,"downloaded_shapefile.zip"))
    unzip(here("Data",SiteID_FileName,'downloaded_shapefile.zip'), exdir=here('Data',SiteID_FileName,'downloaded_shapefile')); unlink(here("Data",SiteID_FileName,'downloaded_shapefile.zip'), recursive = TRUE)
    
    # move shapefile so it can be accessed without workspace ID
    source_folder <- here('Data',SiteID_FileName,'downloaded_shapefile',workspaceID)
    destination_folder <- here('Data',SiteID_FileName,'downloaded_shapefile')
    files <- list.files(source_folder, full.names = TRUE)
    file.rename(files, file.path(destination_folder, basename(files)))
    unlink(source_folder, recursive = TRUE)
  }
  # Load watershed shapefile and get coordinates of centroid
  aoi <<- st_read(here('Data', SiteID_FileName, 'downloaded_shapefile', 'Layers', 'globalwatershed.shp'))
  centroid <- st_centroid(aoi)
  centroid_coords <- st_coordinates(st_transform(centroid, crs = 4326))
  lat <<- centroid_coords[1,'Y']; lon <<- centroid_coords[1,'X']
  return(data.frame(lat=lat, lon=lon, aoi=st_geometry(aoi)))
}



# Define regions (PWR, SER, SWR) based on latitude and longitude
# estimates based on: https://www.fs.usda.gov/rm/pubs_series/rmrs/gtr/rmrs_gtr413.pdf 
get_region <- function(lat, lon){
  if(lat>=41.5 & lat <=49.5 & lon >=-124.5 & lon<=-110.5){
    region='PWR'
  } else if(lat>=24.5 & lat<=40 & lon>=-100 & lon<=-70){
    region='SER'
  } else if(lat<=41 & lat>=31 & lon>=-125 & lon<=-107){
    region='SWR'
  } else{region='other'}
  return(region)
}



# Scrape data from USGS stream gage and aggregate
# Args:
#   GageSiteID
#   incompleteMonths:
#   dataPath
# Returns:
#   xts object of daily streamflow measurements and dataframe of monthly streamflow measurements
get_gage_data <- function(GageSiteID, incompleteMonths, fillLeapDays, dataPath){
  # Scrape data and save
  if(!file.exists(file.path(dataPath, paste0(paste("USGS_Gage",GageSiteID, sep = "_"), ".csv")))){
    DailyStream <- EGRET::readNWISDaily(siteNumber = GageSiteID, parameterCd = "00060") |>
      dplyr::filter(grepl('A', Qualifier)) |> # this filters for any Qualifier that has an A (A and A:E)
      dplyr::mutate(CFS = Q*35.314666212661) # convert Q from m3/s to cfs
    write.csv(DailyStream, file.path(dataPath, paste0(paste("USGS_Gage",GageSiteID, sep = "_"), ".csv")))
  }else{DailyStream<- read.csv(file.path(dataPath, paste0(paste("USGS_Gage",GageSiteID, sep = "_"), ".csv")))}
  DailyStream$date <- as.Date(DailyStream$Date); DailyStream <- DailyStream %>% select(-Date)
  
  # Generate dataframe with complete dates, fill in missing data with NA
  DailyStream <- left_join(data.frame(date = seq(min(DailyStream$date), max(DailyStream$date), by = "day")), 
                           DailyStream, by = "date")
  DailyStream$waterYear <- sapply(DailyStream$date, get_water_year)
  
  # Remove leap days from streamflow data according to user input
  # check this code
  if(!fillLeapDays){
    DailyStream$date<- ymd(DailyStream$date)
    DailyStream <- DailyStream[!(format(DailyStream$date,"%m") == "02" & format(DailyStream$date, "%d") == "29"), , drop = FALSE]
  }
  
  # Extract square mileage of the watershed from the EGRET package
  obj = readNWISInfo(siteNumber = GageSiteID, parameterCd = "00060", interactive = FALSE)
  sqmi <<- obj$drain_area_va
  
  # Aggregate gage discharge data daily and convert from cfs to mm 
  meas_flow_daily <- data.frame(date = DailyStream$date, MeasMM = DailyStream$CFS*28316847*86400/(2590000000000 * sqmi))
  meas_flow_daily_xts <- xts(meas_flow_daily$MeasMM, order.by = ymd(meas_flow_daily$date))
  
  # Aggregate gage discharge data monthly
  meas_flow_daily$YrMon<- format(as.Date(meas_flow_daily$date, format="%Y-%m-%d"),"%Y-%m")
  if(incompleteMonths){
    #sum all months of Measured Discharge, including incomplete months
    meas_flow_mon <- aggregate(meas_flow_daily$MeasMM, by=list(meas_flow_daily$YrMon), FUN=sum)
    colnames(meas_flow_mon)<- c("YrMon", "MeasMM")
  }else{
    #aggregate monthly by anyNA()and sum
    MeasInCompleteMonths<- aggregate(meas_flow_daily$MeasMM, by=list(meas_flow_daily$YrMon), FUN=anyNA)
    meas_flow_mon <- aggregate(meas_flow_daily$MeasMM, by=list(meas_flow_daily$YrMon), FUN=sum)
    
    #merge the incomplete months with the summed months and subset by complete months
    MeasCM <- merge(MeasInCompleteMonths, meas_flow_mon, by = "Group.1", all = TRUE)
    meas_flow_mon = subset(MeasCM, MeasCM$x.x == FALSE)
    
    # clean up
    meas_flow_mon = meas_flow_mon[, c("Group.1","x.y")]
    colnames(meas_flow_mon)<- c("YrMon", "MeasMM")
  }
  
  return(list(meas_flow_daily=meas_flow_daily_xts, meas_flow_mon=meas_flow_mon, DailyStream=DailyStream))
}



# Scrape GridMET meteorological data and clean
# Args:
# Returns: 
#   Dataframe with meteorological data at daily time scale
get_gridmet_data <- function(SiteID_FileName, startY, endY, lat, lon, aoi, dataPath,
                             tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope){
  # Scrape data and save
  if(!file.exists(file.path(dataPath, paste0(paste("GridMET", SiteID_FileName, startY, endY, sep = "_" ), '.csv')))){
    if(point_location) aoi <- data.frame(lon = lon, lat = lat) %>% vect(geom = c("lon", "lat"), crs = "EPSG:4326")
    GridMET_vars <- c("pr", "srad","tmmn", "tmmx", "vpd", "vs")
    DailyClimData <- getGridMET(aoi, varname = GridMET_vars,startDate = startDate, endDate = endDate,verbose = TRUE)
    write.csv(DailyClimData, file.path(dataPath, paste0(paste("GridMET",SiteID_FileName,startY, endY, sep = "_" ), ".csv")), row.names = FALSE) 
  }else{
    DailyClimData = read.csv(file.path(dataPath, paste0(paste("GridMET",SiteID_FileName,startY, endY, sep = "_" ), ".csv")))
  }
  
  # Convert temperature units
  DailyClimData$tmmn<- kelvin_to_celcius(DailyClimData$tmmn); DailyClimData$tmmx<- kelvin_to_celcius(DailyClimData$tmmx)
  
  # Bias adjustment for temperature 
  DailyClimData$tmmn = get_slope_bias_adj(orig = DailyClimData$tmmn, bias = tmmn_bias, slopeadj = tmmn_slope)
  DailyClimData$tmmx = get_slope_bias_adj(orig = DailyClimData$tmmx, bias = tmmx_bias, slopeadj = tmmx_slope)
  
  # Bias adjustment for precip 
  for(i in 1:nrow(DailyClimData)){
    DailyClimData$pr[i] = if (DailyClimData$pr[i] == 0) {
      0
    }else{DailyClimData$pr[i] = get_slope_bias_adj(orig = DailyClimData$pr[i], bias = p_bias, slopeadj = p_slope)}
  }
  
  DailyClimData$date <- as.Date(DailyClimData$date)
  return(DailyClimData)
}



# Scrape Daymet meteorological data and clean
# Args:
# Returns:
#   Dataframe with meteorological data at daily time scale
get_daymet_data <- function(SiteID_FileName, startY, endY, lat, lon, aoi, dataPath){
  # Scrape data and save
  if(!file.exists(file.path(dataPath, paste0(paste("Daymet", SiteID_FileName, startY+1,endY, sep = "_"), ".csv")))){
    if(point_location) aoi <- data.frame(lon = lon, lat = lat) %>% vect(geom = c("lon", "lat"), crs = "EPSG:4326")
    DailyClimData <- getDaymet(aoi, startDate = startDate, endDate = endDate,verbose = TRUE)
    write.csv(DailyClimData, file.path(dataPath, paste0(paste("Daymet",SiteID_FileName,startY, endY, sep = "_" ), ".csv")), row.names = FALSE) 
  } else{
    DailyClimData<- read.csv(file.path(dataPath, paste0(paste("Daymet", SiteID_FileName, startY+1,endY, sep = "_"), ".csv")), skip = 6, header = TRUE, sep = ",")
  }
  
  # Fill leap days according to user input
  # check this code - seems really unwieldy/possibly could be condensed
  HasLeapDays <- data.frame(date = as.Date(seq(0, (nrow(DailyClimData)-1), 1),
                                           origin = ymd(paste(startY+1, startM, startD))))
  NoLeapDays<- HasLeapDays[!(format(HasLeapDays$date,"%m") == "02" & format(HasLeapDays$date, "%d") == "29"), , drop = FALSE]
  DifRows = nrow(HasLeapDays)-nrow(NoLeapDays)
  LastDate = NoLeapDays[nrow(NoLeapDays),]
  LostDates <- data.frame(date = as.Date(seq(1, DifRows, 1), origin = LastDate))
  NoLeapDays<- rbind(NoLeapDays, LostDates)
  row.names(NoLeapDays)<- NULL
  DailyClimData$date<- ymd(NoLeapDays$Date)
  
  if(fillLeapDays){
    # Fill with data from previous day
    DateSeq <- rbind(HasLeapDays, LostDates)
    colnames(DateSeq)<- "date"
    DailyClimData = dplyr::full_join(DateSeq, DailyClimData, by = join_by("date"))
    na_rows <- as.numeric(rownames(DailyClimData[!complete.cases(DailyClimData), ]))
    dateColNumber = which(colnames(DailyClimData)=="date")
    for(i in na_rows){
      DailyClimData[i,-dateColNumber]<- DailyClimData[i-1,-dateColNumber]
    }
  }
  
  # Match format of GridMET data
  DailyClimData$month<- as.numeric(format(as.Date(DailyClimData$date, format="%Y-%m-%d"),"%m"))
  DailyClimData<- subset(DailyClimData, DailyClimData$date<=endDate)
  DailyClimData$year<- NULL
  DailyClimData$yday<- NULL
  DailyClimData$swe..kg.m.2.<- NULL
  if(fillLeapDays){ #order of columns was changed because of merging 
    colnames(DailyClimData)<- c("Date", "dayl..s." , "pr", "srad" ,"tmmx", "tmmn", "vp..Pa.", "month")
  }else{colnames(DailyClimData)<- c("dayl..s." , "pr","srad" ,"tmmx", "tmmn", "vp..Pa." ,"Date", "month")}
  
  return(DailyClimData)
}



# Pull MACA projections for a single point and clean data
# Args:
# Returns:
get_maca_point <- function(lat, lon, SiteID_FileName){
  # Pull future meteorological data
  if(!file.exists(here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100_point.csv', sep='_')))){
    # Pull data
    point <- data.frame(lon = lon, lat = lat) %>% vect(geom = c("lon", "lat"), crs = "EPSG:4326")
    future_climate_data <- getMACA(point, c('tasmin','tasmax','pr','rsds','vpd','vas','uas'), timeRes='day', model=gcm_list, scenario=c('rcp45','rcp85'), 
                                   startDate = '2023-01-01', endDate = '2099-12-31')
    
    # Clean and compile data
    precip<-NULL; tasmin<-NULL; tasmax<-NULL; rsds<-NULL; vpd<-NULL; vas<-NULL; uas<-NULL 
    for(i in 2:length(colnames(future_climate_data))){
      split_colnames <- strsplit(colnames(future_climate_data[i]), "_")
      combine <- data.frame(date=future_climate_data[,1], GCM=split_colnames[[1]][2], RCP=split_colnames[[1]][4],
                            var=future_climate_data[,i])
      if(split_colnames[[1]][1] == 'pr') {precip <- rbind(precip, combine)}
      if(split_colnames[[1]][1] == 'tasmin') {tasmin <- rbind(tasmin, combine)}
      if(split_colnames[[1]][1] == 'tasmax') {tasmax <- rbind(tasmax, combine)}
      if(split_colnames[[1]][1] == 'rsds') {rsds <- rbind(rsds, combine)}
      if(split_colnames[[1]][1] == 'vpd') {vpd <- rbind(vpd, combine)}
      if(split_colnames[[1]][1] == 'vas') {vas <- rbind(vas, combine)}
      if(split_colnames[[1]][1] == 'uas') {uas <- rbind(uas, combine)}
    }
    colnames(precip) <- c('date','GCM','RCP','pr'); colnames(tasmin) <- c('date','GCM','RCP','tmmn'); colnames(tasmax) <- c('date','GCM','RCP','tmmx')
    colnames(rsds) <- c('date','GCM','RCP','srad'); colnames(vpd) <- c('date','GCM','RCP','vpd'); colnames(vas) <- c('date','GCM','RCP','vas'); colnames(uas) <- c('date','GCM','RCP','uas')
    
    future_climate <- Reduce(function(x, y) merge(x, y, by = c('date', 'GCM', 'RCP'), all = TRUE), list(precip, tasmin, tasmax, rsds, vpd, vas, uas))
    
    future_climate$tmmn <- future_climate$tmmn - 273.15; future_climate$tmmx <- future_climate$tmmx - 273.15
    future_climate$vs <- sqrt(future_climate$vas^2 + future_climate$uas^2)
    future_climate <- future_climate %>% mutate(projection = paste0(GCM, '.', RCP))
    future_climate <- future_climate %>% select(-vas, -uas, -GCM, -RCP)
    future_climate <- future_climate[,c('projection','date','pr','srad','tmmn','tmmx','vs','vpd')]
    future_climate$date <- as.Date(future_climate$date)
    
    # Save
    write.csv(future_climate, file = here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100_point.csv', sep='_')), row.names = FALSE)
  } else {
    future_climate <- read.csv(here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100_point.csv', sep='_')))
    future_climate$date <- as.Date(future_climate$date)
  }
  return(future_climate)
}



# Pull MACA projections for an area of interest and clean data
# Args:
# Returns:
get_maca_data_area <- function(aoi, SiteID_FileName){
  if(!file.exists(here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100_area.csv', sep='_')))){
    # Pull data
    future_climate_data <- getMACA(aoi, c('pr','rsds','vpd','vas','uas'), timeRes='day', model=gcm_list, scenario=c('rcp45','rcp85'), 
                                   startDate = '2023-01-01', endDate = '2099-12-31')
    future_climate_data_tasmin <- getMACA(aoi, 'tasmin', timeRes='day', model=gcm_list, scenario=c('rcp45','rcp85'), 
                                         startDate = '2023-01-01', endDate = '2099-12-31')
    future_climate_data_tasmax <- getMACA(aoi, 'tasmax', timeRes='day', model=gcm_list, scenario=c('rcp45','rcp85'), 
                                          startDate = '2023-01-01', endDate = '2099-12-31')
    future_climate_data$tasmin <- future_climate_data_tasmin$air_temperature; future_climate_data$tasmax <- future_climate_data_tasmax$air_temperature
    
    # Check all the data is the correct length and re-query if not
    dict <- list(precipitation='pr', tasmin='tasmin', tasmax='tasmax', surface_downwelling_shortwave_flux_in_air='rsds', vpd='vpd', vas='northward_wind', uas='eastward_wind')
    for(var in c("precipitation","tasmin","tasmax","surface_downwelling_shortwave_flux_in_air","vpd","northward_wind","eastward_wind")){
      while(nlyr(future_climate_data[[var]]) != 731224){
        print(var)
        new_future_climate_data <- getMACA(aoi, dict[[var]], timeRes='day', model=gcm_list, scenario=c('rcp45','rcp85'), 
                       startDate = '2023-01-01', endDate = '2099-12-31')
        if(var=='tasmax' | var=='tasmin'){
          future_climate_data[[var]] <- new_future_climate_data$air_temperature
        }
        future_climate_data[[var]] <- new_future_climate_data[[var]]
      }
    }
    
    # Pull out each meteorological variable
    precip <- global(future_climate_data$precipitation, fun = "mean", na.rm = TRUE); colnames(precip) <- c('pr')
    precip <- precip %>% rownames_to_column("rowname") %>% separate(rowname, into = c("Variable", "date","GCM","Run","RCP"), sep = "_") %>% select(-Variable, -Run)
    tmmn <- global(future_climate_data$tasmin, fun = "mean", na.rm = TRUE); colnames(tmmn) <- c('tmmn')
    tmmn <- tmmn %>% rownames_to_column("rowname") %>% separate(rowname, into = c("Variable", "date","GCM","Run","RCP"), sep = "_") %>% select(-Variable, -Run)
    tmmx <- global(future_climate_data$tasmax, fun = "mean", na.rm = TRUE); colnames(tmmx) <- c('tmmx')
    tmmx <- tmmx %>% rownames_to_column("rowname") %>% separate(rowname, into = c("Variable", "date","GCM","Run","RCP"), sep = "_") %>% select(-Variable, -Run)
    tmmx$tmmx <- tmmx$tmmx - 273.15; tmmn$tmmn <- tmmn$tmmn - 273.15
    srad <- global(future_climate_data$surface_downwelling_shortwave_flux_in_air, fun = "mean", na.rm = TRUE); colnames(srad) <- c('srad')
    srad <- srad %>% rownames_to_column("rowname") %>% separate(rowname, into = c("Variable", "date","GCM","Run","RCP"), sep = "_") %>% select(-Variable, -Run)
    vpd <- global(future_climate_data$vpd, fun = "mean", na.rm = TRUE); colnames(vpd) <- c('vpd')
    vpd <- vpd %>% rownames_to_column("rowname") %>% separate(rowname, into = c("Variable", "date","GCM","Run","RCP"), sep = "_") %>% select(-Variable, -Run)
    vas <- global(future_climate_data$northward_wind, fun = "mean", na.rm = TRUE); colnames(vas) <- c('vas')
    vas <- vas %>% rownames_to_column("rowname") %>% separate(rowname, into = c("Variable", "date","GCM","Run","RCP"), sep = "_") %>% select(-Variable, -Run)
    uas <- global(future_climate_data$eastward_wind, fun = "mean", na.rm = TRUE); colnames(uas) <- c('uas')
    uas <- uas %>% rownames_to_column("rowname") %>% separate(rowname, into = c("Variable", "date","GCM","Run","RCP"), sep = "_") %>% select(-Variable, -Run)
    
    # Merge and forcibly correct discrepancies in order
    future_climate <- precip %>% select(-pr) %>% arrange(date, tolower(GCM), RCP)
    future_climate$pr <- precip$pr; future_climate$tmmx <- tmmx$tmmx; future_climate$tmmn <- tmmn$tmmn; future_climate$srad <- srad$srad; future_climate$vpd <- vpd$vpd; future_climate$vas <- vas$vas; future_climate$uas <- uas$uas

    # Clean data
    future_climate$vs <- sqrt(future_climate$vas^2 + future_climate$uas^2)
    future_climate <- future_climate %>% mutate(projection = paste0(GCM, '.', RCP))
    future_climate <- future_climate %>% select(-vas, -uas, -GCM, -RCP)
    future_climate <- future_climate[,c('projection','date','pr','srad','tmmn','tmmx','vs','vpd')]
    future_climate$date <- as.Date(future_climate$date)
    
    # Save
    write.csv(future_climate, file = here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100_area.csv', sep='_')), row.names = FALSE)
  } else {
    future_climate <- read.csv(here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100_area.csv', sep='_')))
    future_climate$date <- as.Date(future_climate$date)
  }
  return(future_climate)
}



# Pull gridded water balance data for a single point from CONUS model
# More details about the gridded product can be found:
# https://www.yellowstoneecology.com/research/Gridded_Water_Balance_Model_Version_2_User_Manual.pdf
# Args:
#   SiteID_FileName
#   lat, lon: Latitude and longitude of site, in degrees
#   startY_future, endY_future: start and end years of future projection period
# Returns:
#   Gridded water balance data for CONUS, which is saved as a csv file
# NOTES: 
# Mike Tercek's website is not consistently up and running so this will likely not work
# "agdd", fix this - return when agdd is back on the website; make sure to put agdd BEFORE AET - will mess up code if agdd is last
# MODIFIED from Janelle/Connor code
get_conus_wb <- function(SiteID_FileName, lat, lon, startY_future, endY_future){
  # Return file if it exists
  if(file.exists(file.path(dataPath, paste("WB_conus",SiteID_FileName,"2023_2100.csv", sep = "_")))){
    future_wb <- read.csv(file.path(dataPath, paste("WB_conus",SiteID_FileName,"2023_2100.csv", sep = "_")))
    future_wb$date <- as.Date(future_wb$date)
    # not sure if I need this
    #future_wb$adj_runoff<- get_adj_runoff(future_wb$runoff, gw_add = gw_add, vfm = vfm)
    return(future_wb)
  }
  
  # Pull data from online for each GCM, RCP, and year
  future_wb <- NULL
  for(GCM in gcm_list){
    for(RCP in c("rcp85", "rcp45")){
      print(paste("downloading", GCM, RCP))
      model_holder <- NULL
      for(yr in c(startY_future:endY_future)){ 
        leap <- lubridate::leap_year(yr) 
        if(leap == TRUE){enddate <- paste(yr,"12-31", sep = "-")
        } else {enddate <- paste(yr + 1, "01-01", sep = "-")}
        var_holder <- NULL
        for(climvar in c("soil_water", "runoff", "rain","accumswe", "PET", "Deficit", "AET")){
          holder <- NULL 
          data_url<-paste("http://www.yellowstone.solutions/thredds/ncss/daily_or_monthly/gcm/",RCP,"/",GCM,"/V_1_5_",yr,"_",GCM,"_",RCP,"_",climvar,".nc4?var=",climvar,"&latitude=",lat,"&longitude=",lon,"&time_start=",yr,"-01-01T12%3A00%3A00Z&time_end=",enddate,"T12%3A00%3A00Z&accept=csv_file",sep ="") 
          
          # test url: http://www.yellowstone.solutions/thredds/ncss/daily_or_monthly/gcm/rcp85/inmcm4/V_1_5_2099_inmcm4_rcp85_soil_water_monthly.nc4?var=soil_water&latitude=45&longitude=-111&time_start=2099-01-16T05%3A14%3A31.916Z&time_end=2099-12-17T00%3A34%3A14.059Z&accept=csv_file
          
          # Catch error if Mike's thredds server is not running
          holder <- result <- tryCatch(
            {
              holder <- data.frame(suppressWarnings(fread(data_url, verbose=FALSE, showProgress = FALSE)))
              return(holder)
            }, error = function(e) {
              calcFutureWB <<- TRUE
              message('The water balance server is not currently running. Please contact Mike Tercek (miketercek@yahoo.com) for more information.')
              return(NA)
              }
          )
          if(is.na(holder)){return(holder)}
           
          colnames(holder) <- c("date", "latitude", "longitude", paste(climvar))
          holder$GCM <- GCM; holder$RCP <- RCP
          if(is.null(var_holder)){
            var_holder <- holder
          } else{
            var_holder <- merge(var_holder, holder)}
        }
        model_holder <- rbind(model_holder, var_holder)
      }
      future_wb <- rbind(future_wb, model_holder)
    }
    # Clean data
    future_wb <- future_wb[future_wb$soil_water != -32767, ]
    future_wb$time <- as.Date(future_wb$time)
    future_wb$projection <- paste(future_wb$GCM, future_wb$RCP, sep='.')
    future_wb <- subset(future_wb, select = -c(latitude, longitude, GCM, RCP))
    
    # Convert from mm to in 
    future_wb <- future_wb %>% dplyr::mutate(across(where(is.numeric), ~ . / 25.4))
    
    # Get adjusted runoff
    future_wb$adj_runoff<- get_adj_runoff(future_wb$runoff, gw_add = gw_add, vfm = vfm)
  }
  write.csv(future_wb, file.path(dataPath, paste("WB_conus",SiteID_FileName,"2023_2100.csv", sep = "_")))
  return(future_wb)
}



# Pull gridded water balance data from a file provided by Mike Tercek
get_conus_wb_direct <- function(SiteID_FileName, dataPath, filename){
  if(file.exists(filename)){
    future_wb_conus <- read.csv(filename)
    
    # Clean: fix date, convert to mm, fix column names
    future_wb_conus <- future_wb_conus %>% rename(date = Date); future_wb_conus$date <- as.Date(future_wb_conus$date)
    future_wb_conus<-cbind(future_wb_conus[,c("date","GCM")], 25.4*(future_wb_conus[,c(which(colnames(future_wb_conus)=="Deficit.in"):ncol(future_wb_conus))]))
    colnames(future_wb_conus)<- c("date", "projection", "Deficit", "AET", "soil_water", "runoff", "rain", "accumswe", "PET")
    #future_wb<-subset(future_wb, projection !="MIROC-ESM-CHEM.rcp85")   # drop "MIROC-ESM-CHEM.rcp85" because it doesn't have an associated RCP 4.5  
    
    # adjust for ground water addition and volume forcing multiplier
    future_wb_conus$adj_runoff<- get_adj_runoff(future_wb_conus$runoff, gw_add = gw_add, vfm = vfm)
    
    # save 
    write.csv(future_wb_conus, file.path(dataPath, paste("WB_conus",SiteID_FileName,"2023_2100.csv", sep = "_")), row.names=FALSE)
    return(future_wb_conus)
  }
  return(NA)
}



# Pull OpenET data for a single point
get_et_point <- function(startY, startM, startD, endY, endM, endD, siteID_FileName, interval, dataPath, api_key){
  file_path <- here(dataPath, paste0(paste("OpenET", interval, SiteID_FileName, startY, endY, sep = "_" ), '.csv'))
  if(!file.exists(file_path)){
    header <- add_headers(accept = 'application/json', Authorization = api_key, content_type = 'application/json')
    
    # get data
    args <- list(date_range = c(paste(startY, sprintf("%02d", startM), sprintf("%02d", startD), sep='-'), 
                                paste(endY, sprintf("%02d", endM), sprintf("%02d", endD), sep='-')),
                 interval = interval,geometry = c(lon, lat), model = "Ensemble", variable = "ET", reference_et = "gridMET",
                 units = "mm", file_format = "JSON")
    response <- POST(url = "https://openet-api.org/raster/timeseries/point", header, body = args, encode = "json")
    
    # Check if request was successful
    if(response$status_code == 200){
      ET <- data.frame(fromJSON(content(response, as = "text", encoding = "UTF-8")))
      colnames(ET) <- c('date', 'Meas ET')
      write.csv(ET, file_path)
    } else{
      print(paste('No', interval, 'OpenET data for that region or time period. Optimization cannot occur.')); stop()
    }
  } else {
    ET <- read.csv(file_path, '.csv')
  }
  return(ET)
}



# Function to calculate pseudo R-squared
calculate_pseudo_r2 <- function(model, newdata, tau){
  preds <- predict(model, newdata = newdata)
  tau_preds <- preds[, which(model$tau == tau)]
  residuals <- newdata$Meas - tau_preds
  rss <- sum(residuals^2)
  tss <- sum((newdata$Meas - mean(newdata$Meas))^2)
  r2 <- 1 - (rss / tss)
  return(r2)
}



# Function to return water year
get_water_year <- function(date){
  if(month(date) < 10){
    return(year(date))
  } else {
    return(year(date) + 1)
  }
}



# Identify whether redundant climate futures were selected
ID.redundant.gcm <- function(PCA){
  #ID redundant diagonal
  redundant.diag=count(PCA,diagonals)$diagonals[which(count(PCA,diagonals)$n==1)] 
  
  #ID which PC has the redundant diagonal
  PC.foul = PCA$PC[which(PCA$diagonals == redundant.diag)] 
  
  #ID projection that is in both the  redundant diagonal and the duplicative PC
  PCA$projection[which(PCA$PC == PC.foul & PCA$projection != PCA$projection[which(PCA$diagonals == redundant.diag)])] 
}

  