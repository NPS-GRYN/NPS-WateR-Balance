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


# Get latitude and longitude coordinates of watershed centroid using StreamStats database
# Args:
#   SiteID_FileName:
#   GageSiteID:
# Returns:
#   Latitude and longitude coordiantes of watershed centroid
get_coords <- function(SiteID_FileName, GageSiteID){
  if(!file.exists(here('Data', SiteID_FileName, 'downloaded_shapefile', 'Layers', 'globalwatershed.shp'))){
    # get lat/lon data for watershed
    INFO <- readNWISInfo(siteNumber=GageSiteID, parameterCd = "")
    
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
  return(data.frame(lat=lat, lon=lon))
}



# Scrape data from USGS stream gage and aggregate
# Args:
#   GageSiteID
#   incompleteMonths:
#   dataPath
# Returns:
#   xts object of daily streamflow measurements and dataframe of monthly streamflow measurements
get_gage_data <- function(GageSiteID, incompleteMonths, dataPath){
  # Scrape data and save
  if(!file.exists(file.path(dataPath, paste0(paste("USGS_Gage",GageSiteID, startY+1,endY, sep = "_"), ".csv")))){
    DailyStream <- EGRET::readNWISDaily(siteNumber = GageSiteID, parameterCd = "00060", 
                                        startDate = paste(startY, startM, startD, sep='-'), endDate = paste(endY, endM, endD, sep='-')) |>
      dplyr::filter(grepl('A', Qualifier)) |> #this filters for any Qualifier that has an A. It will return A and A:E
      dplyr::mutate(CFS = Q*35.314666212661) #converts Q to flow cfs
    write.csv(DailyStream, file.path(dataPath, paste0(paste("USGS_Gage",GageSiteID, startY+1,endY, sep = "_"), ".csv")))
  }else{DailyStream<- read.csv(file.path(dataPath, paste0(paste("USGS_Gage",GageSiteID, startY+1,endY, sep = "_"), ".csv")))}
  DailyStream$Date <- as.Date(DailyStream$Date)
  
  # Extract square mileage of the watershed from the EGRET package
  obj = readNWISInfo(siteNumber = GageSiteID, parameterCd = "00060", interactive = FALSE)
  sqmi = obj$drain_area_va
  
  # Aggregate gage discharge data daily and convert from cfs to mm 
  meas_flow_daily <- data.frame(Date = DailyStream$Date, MeasMM = DailyStream$CFS*28316847*86400/(2590000000000 * sqmi))
  meas_flow_daily_xts <- xts(meas_flow_daily$MeasMM, order.by = ymd(meas_flow_daily$Date))
  
  # Aggregate gage discharge data monthly
  meas_flow_daily$YrMon<- format(as.Date(meas_flow_daily$Date, format="%Y-%m-%d"),"%Y-%m")
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
  
  # Remove leap days from streamflow data according to user input
  # check this code
  if(!fillLeapDays){
    DailyStream$Date<- ymd(DailyStream$Date)
    DailyStream <- DailyStream[!(format(DailyStream$Date,"%m") == "02" & format(DailyStream$Date, "%d") == "29"), , drop = FALSE]
  }
  
  return(list(meas_flow_daily_xts=meas_flow_daily_xts, meas_flow_mon=meas_flow_mon))
}



# Scrape GridMET meteorological data and clean
# Args:
# Returns: 
#   Dataframe with meteorological data at daily time scale
get_gridmet_data <- function(SiteID_FileName, startY, endY, lat, lon, dataPath,
                             tmmn_bias, tmmn_slope, tmmx_bias, tmmx_slope, p_bias, p_slope){
  # Scrape data and save
  if(!file.exists(file.path(dataPath, paste0(paste("GridMET", SiteID_FileName, startY, endY, sep = "_" ), '.csv')))){
    point <- data.frame(lon = lon, lat = lat) %>% vect(geom = c("lon", "lat"), crs = "EPSG:4326")
    GridMET_vars <- c("pr", "srad","tmmn", "tmmx", "vpd", "vs")
    DailyClimData <- getGridMET(point,varname = GridMET_vars,startDate = startDate, endDate = endDate,verbose = TRUE)
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
  return(DailyClimData)
}



# Scrape Daymet meteorological data and clean
# Args:
# Returns:
#   Dataframe with meteorological data at daily time scale
get_daymet_data <- function(SiteID_FileName, startY, endY, lat, lon, dataPath){
  # Scrape data and save
  if(!file.exists(file.path(dataPath, paste0(paste("Daymet", SiteID_FileName, startY+1,endY, sep = "_"), ".csv")))){
    point <- data.frame(lon = lon, lat = lat) %>% vect(geom = c("lon", "lat"), crs = "EPSG:4326")
    DailyClimData <- getDaymet(point, startDate = startDate, endDate = endDate,verbose = TRUE)
    write.csv(DailyClimData, file.path(dataPath, paste0(paste("Daymet",SiteID_FileName,startY, endY, sep = "_" ), ".csv")), row.names = FALSE) 
  } else{
    DailyClimData<- read.csv(file.path(dataPath, paste0(paste("Daymet", SiteID_FileName, startY+1,endY, sep = "_"), ".csv")), skip = 6, header = TRUE, sep = ",")
  }
  
  # Fill leap days according to user input
  # check this code - seems really unwieldy/possibly could be condensed
  HasLeapDays <- data.frame(Date = as.Date(seq(0, (nrow(DailyClimData)-1), 1),
                                           origin = ymd(paste(startY+1, startM, startD))))
  NoLeapDays<- HasLeapDays[!(format(HasLeapDays$Date,"%m") == "02" & format(HasLeapDays$Date, "%d") == "29"), , drop = FALSE]
  DifRows = nrow(HasLeapDays)-nrow(NoLeapDays)
  LastDate = NoLeapDays[nrow(NoLeapDays),]
  LostDates <- data.frame(Date = as.Date(seq(1, DifRows, 1), origin = LastDate))
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
    colnames(DailyClimData)<- c("date", "dayl..s." , "pr", "srad" ,"tmmx", "tmmn", "vp..Pa.", "month")
  }else{colnames(DailyClimData)<- c("dayl..s." , "pr","srad" ,"tmmx", "tmmn", "vp..Pa." ,"date", "month")}
  
  return(DailyClimData)
}



# Pull MACA projections for a single point and clean data
# Args:
# Returns:
get_maca_point <- function(){
  # Pull future meteorological data
  if(!file.exists(here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100.csv', sep='_')))){
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
    
    # Save
    write.csv(future_climate, file = here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100.csv', sep='_')), row.names = FALSE)
  } else {
    future_climate <- read.csv(here('Data', SiteID_FileName, paste('MACA', SiteID_FileName, endY, '2100.csv', sep='_')))
  }
  return(future_climate)
}



# Pull MACA projections for an area of interest and clean data
# ADD this
get_maca_data_area <- function(){
  
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
# "agdd", fix this - return when agdd is back on the website; make sure to put agdd BEFORE AET - will mess up code if agdd is last
# MODIFIED from Janelle/Connor code
# STILL EDITING: mike tercek's website appears to be down. next step: TEST
get_conus_wb <- function(SiteID_FileName, lat, lon, startY_future, endY_future){
  future_wb <- NULL
  # Loop through and pull data for each GCM, RCP, year, and variable
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
          
          holder <-data.frame(fread(data_url, verbose=FALSE, showProgress = FALSE,)) 
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
    
    # convert from mm to in 
    future_wb <- future_wb %>% dplyr::mutate(across(where(is.numeric), ~ . / 25.4))
  }
  write.csv(future_wb, file.path(dataPath, paste("WB_conus",SiteID_FileName,"2023_2100.csv", sep = "_")))
  return(future_wb)
}



# Pull OpenET data for a single point
get_et_point <- function(startY, startM, startD, endY, endM, endD, siteID_FileName, interval, dataPath){
  file_path <- here(dataPath, paste0(paste("OpenET", interval, SiteID_FileName, startY, endY, sep = "_" ), '.csv'))
  if(!file.exists(file_path)){
    api_key <- 'ZZjI9EAHEFsVhFf8WVgBD2J6ks14IbJZJgHYR1iBPO82EcYO2XxeDJAcwAN9'
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
      colnames(ET) <- c('Date', 'Meas ET')
      write.csv(ET, file_path)
    } else{
      print(paste('No', interval, 'OpenET data for that region or time period. Optimization cannot occur.')); stop()
    }
  } else {
    ET <- read.csv(file_path, '.csv')
  }
  return(ET)
}

  