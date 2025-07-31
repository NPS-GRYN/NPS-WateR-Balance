# ---------------------------------------------------------------------
# This script includes code to generate future projections of streamflow using the IHACRES 
# rainfall-streamflow methodology. This requires future projections of adjusted runoff,
# which are either calculated using pre-calibrated coefficients or pulled from a pre-generated
# gridded water balance model produced/maintained by Mike Tercek. This code also provides 
# preliminary visualizations of these future streamflow projections. 
# 
# ---------------------------------------------------------------------



#######################################################################
### GENERATE FUTURE WATER BALANCE MODELS ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')

# Remove low skill models using list from Rupp et al. 2016
low_skill_models = read.delim('./Data/GCM_skill_by_region.txt', header=TRUE) %>% 
  filter(Region == ifelse(region %in% Region, region, "mean")) %>% top_n(n=round(length(gcm_list)*percent_skill_cutoff), wt=Rank)
gcm_list <- gcm_list[!gcm_list %in% low_skill_models$GCM]


### Use Mike Tercek's pre-generated gridded CONUS water balance model for future projections ###
if(!calcFutureWB){
  # Pull directly from website
  future_wb_conus <- get_conus_wb(SiteID_FileName, lat, lon, endY, 2099)
  
  # Uncomment to access file (obtained directly from Mike Tercek)
  # Add name of path
  #filename_future_wb = ""  
  #if(file.exists(filename_future_wb)) future_wb_conus <- get_conus_wb_direct(SiteID_FileName, dataPath, filename_future_wb)
  
  # If neither version can be accessed, calculate water balance
  if(!exists("future_wb_conus")){
    calcFutureWB <- TRUE
  }
}



### Re-run water balance model to generate future projections ###
# for entire MACA period of record (2006-2100)
if(calcFutureWB){
  if(!file.exists(here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, "2006_2100.csv", sep='_')))){
    # Get future climate data
    if(point_location) {future_climate <- get_maca_point(lat, lon, 2006, 2100, SiteID_FileName, gcm_list)
    } else future_climate <- get_maca_area(aoi, 2006, 2100, SiteID_FileName, gcm_list)
    
    # Run water balance code for each future projection
    wb_list <- vector("list", length(unique(future_climate$projection)))
    i <- 1
    for(proj in unique(future_climate$projection)){
      print(proj)
      ClimData <- future_climate %>% filter(projection==proj) %>% select(-projection)
      DailyWB_future <- WB(ClimData, gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope,
                           shade.coeff, jtemp, SWC.Max, 1, 1, 0, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon, "")
      wb_list[[i]] <- cbind(proj, DailyWB_future)
      i <- i + 1
    }
    future_wb_calc <- do.call(rbind, wb_list); future_wb_calc <- future_wb_calc %>% rename(projection=proj)
    
    
    # Save calculated WB
    write.csv(future_wb_calc, here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, "2006_2100.csv", sep='_')), row.names=FALSE)
  } else{
    # Read in calculated WB 
    future_wb_calc <- read.csv(here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, "2006_2100.csv", sep='_')))
    future_wb_calc$date <- as.Date(future_wb_calc$date)
  }
}



### Compare the two future water balance projections, just for fun ###
# plot AET, deficit, adj runoff for a sample year and a sample model
# edit if statement
if(exists("future_wb_conus") & exists("future_wb_calc")){
  model_run = 'HadGEM2-CC365.rcp45'; yr = 2060
  if(make_plots){
    plot_aet <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=AET), col='black')+
      geom_line(data=future_wb_calc %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=AET), col='red')+
      labs(x='Date',y='AET [mm]', title='Actual Evapotranspiration') +
      theme(legend.position = "none") + nps_theme()
    plot_d <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=Deficit), col='black')+
      geom_line(data=future_wb_calc %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=D), col='red')+
      labs(x='Date',y='Deficit [mm]', title='Deficit') +
      theme(legend.position = "none") + nps_theme()
    plot_run <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=adj_runoff, col='Gridded WB'))+
      geom_line(data=future_wb_calc %>% filter(projection==model_run& year(date)==yr), aes(x=date, y=adj_runoff, col='Calculated WB'))+
      labs(x='Date',y='Adjusted Runoff [mm]', title='Adjusted Runoff') +  
      scale_color_manual(values = c("Gridded WB"="black", "Calculated WB"="red"), name="WB Projections") + nps_theme()
    
    nameReduce = gsub(pattern = " ",replacement = "_", x = paste(SiteID, "Future WB Projection Comparison"))
    jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=2000, height=600)
    grid.arrange(plot_aet, plot_d, plot_run, ncol = 3, widths=c(1,1,1.3), top = textGrob(paste('WB Projection Comparisons for', model_run, ':', yr),gp=gpar(fontsize=30)))
    dev.off() 
  }
}



#######################################################################
### GENERATE FUTURE STREAMFLOW PROJECTIONS ###

# Define which future water balance projection to use
if(calcFutureWB){
  future_wb <- future_wb_calc
  if(!dir.exists(file.path(outLocationPath, 'WB_Calc'))) {dir.create(file.path(outLocationPath, 'WB_Calc'))}; outLocationPathFuture = file.path(outLocationPath, 'WB_Calc')
} else{
  future_wb <- future_wb_conus
  if(!dir.exists(file.path(outLocationPath, 'WB_CONUS'))) {dir.create(file.path(outLocationPath, 'WB_CONUS'))}; outLocationPathFuture = file.path(outLocationPath, 'WB_CONUS')
}

# Run IHACRES model for each future projection
#if(!file.exists(paste(outLocationPathFuture,'Future_Streamflow_Projections',SiteID_FileName, endY, "2100.csv", sep='_'))){
gcms<-unique(future_wb$projection)
futures <- NULL
for (j in 1:length(gcms)){
  # Subset one model
  fut_ro <- subset(future_wb, projection == gcms[j]); print(gcms[j])
  data <- data.frame(fut_ro$adj_runoff)
  colnames(data) <- c("adj_runoff")
  
  # Run IHACRES model
  DailyDrainFuture <- Drain(data, q0, s0, v0, qa, qb, sa, sb, va, vb)
  
  # Save streamflow projection to futures dataframe
  drainage_qsvt <- cbind(gcms[j], fut_ro$date, DailyDrainFuture)
  colnames(drainage_qsvt)[] <- c("projection", "date", "adj_runoff", "quick", "slow", "veryslow", "total")
  futures <-rbind(futures, drainage_qsvt)
}

# Extract GCM and RCP
futures$gcm <- sapply(X = strsplit(futures$projection, split=".rcp"), FUN = "[", 1) 
futures$rcp <- as.numeric(x = sapply(strsplit(futures$projection, split=".rcp"), FUN = "[", 2)) 

# Save future streamflow projections
write.csv(futures, here('Data',SiteID_FileName,paste('Future_Streamflow_Projections',SiteID_FileName, endY, "2100.csv", sep='_')), row.names=FALSE)
#} else{
#  futures <- read.csv(here('Data',SiteID_FileName,paste('Future_Streamflow_Projections',SiteID_FileName, endY, "2100.csv", sep='_')))
#  futures$date <- as.Date(futures$date)
#}

#######################################################################
### COMPILE FUTURE AND HISTORICAL PROJECTIONS ###

# Get historical model projections
hist_flow_mod <- data.frame(date=DailyDrain$date, adj_runoff=DailyDrain$adj_runoff, quick=DailyDrain$Quick,
                            slow=DailyDrain$Slow, veryslow=DailyDrain$Very_Slow, total=DailyDrain$total, 
                            projection = rep('Historical', length(DailyDrain$date)), gcm=rep('Historical', length(DailyDrain$date)), rcp=rep('Hist', length(DailyDrain$date)))


# Combine historical and future model projections
daily_df<-rbind(futures, hist_flow_mod)
daily_df$date <- as.Date(daily_df$date); daily_df$yr<-as.numeric(format(daily_df$date,"%Y")); daily_df$mo<-format(daily_df$date,"%m"); daily_df$yr_mo<-format(daily_df$date,"%Y-%m"); daily_df$day <- factor(yday(daily_df$date), levels = unique(yday(daily_df$date)))
daily_df$water_year <- sapply(daily_df$date, get_water_year); 
daily_df <- daily_df %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
daily_df<-daily_df[,c("date", "projection", "gcm", "rcp", "yr", "mo", "yr_mo", "water_year", "adj_runoff", "quick", "slow", "veryslow", "total")]
daily_df$Period<-ifelse(daily_df$yr<=2005,"Historical",ifelse (daily_df$yr>=2006 & daily_df$yr<=2050,"Early",
                                                               ifelse (daily_df$yr>=2051 & daily_df$yr<=2070,"Middle", ifelse (daily_df$yr>=2071, "Late","NA"))))


# Aggregate data to annual
annual_df <- as.data.frame(daily_df %>% group_by(gcm, rcp, yr) %>%
  dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), water_year=first(water_year),
                   adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                   slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE), total = sum(total, na.rm = TRUE), 
                   Period=first(Period)))
annual_df$date<-as.Date(paste(annual_df$yr,"-01", "-01",sep=""))

# Aggregate data to monthly 
monthly_df <- as.data.frame(daily_df %>% group_by(gcm,rcp, yr_mo) %>%
  dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), water_year=first(water_year),
                   adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                   slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE) ,total = sum(total, na.rm = TRUE),
                   Period=first(Period)))
monthly_df$date<-as.Date(paste(monthly_df$yr_mo,"-01",sep=""))

# Mean of future daily streamflow projections
mean_daily_df <- daily_df %>% filter(projection!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total))



#######################################################################
### SELECT DIVERGENT CLIMATE FUTURES ###
#scenario_names <- c("Warm Wet", "Hot Wet", "Warm Dry", "Hot Dry")
#color_names <- c("#12045C","#E10720","#9A9EE5","#F3D3CB")
color_names <- c("#4E5BA6","#E10720","#B4B8F0","#E79C9C")
#color_names <- c('#7570B3','#E69F00','#92A8D1','#D55E00')
future_means <- select_climate_futures(color_names)


### Identify models in format for plotting ###
model_names <- (future_means %>% filter(!is.na(pca)))$projection; scenario_names <- (future_means %>% filter(!is.na(pca)))$pca
num_models <- length(model_names)

# identify warm wet/hot dry [1:2] or warm dry/hot wet [3:4] or all (comment out both lines)
#model_names <- model_names[1:2]; scenario_names <- scenario_names[1:2]; color_names <- color_names[1:2]
#model_names <- model_names[3:4]; scenario_names <- scenario_names[3:4]; color_names <- color_names[3:4]



#######################################################################
### GET MACA DATA FOR CLIMATE FUTURES ###
# edit so only 4 climate futures are selected
if(point_location){
  if(!file.exists(here('Data',SiteID_FileName,paste('Historical_Streamflow_Projections',SiteID_FileName, "1960_2005_point.csv", sep='_')))){
    maca_hist <- get_maca_hist_point(lat, lon, 1960, 2005, SiteID_FileName, model_names)
    maca_hist$gcm <- maca_hist$projection; maca_hist$projection <- paste0(maca_hist$projection, ".Hist")
    
    
    # run streamflow model
    hist <- NULL
    for(mod in model_names){
      hist_clim <- maca_hist %>% filter(grepl(strsplit(mod, "\\.")[[1]][1], gcm))
      
      hist_wb <- WB(hist_clim, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
                    shade.coeff, jtemp,SWC.Max, 1, 1, 0, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, "")
      hist_wb <- Drain(hist_wb, q0, s0, v0, qa, qb, sa, sb, va, vb)
      hist_wb <- hist_wb %>% select(projection, date, adj_runoff, Quick, Slow, Very_Slow, total, gcm)
      hist_wb <- hist_wb %>% rename(quick=Quick, slow=Slow, veryslow=Very_Slow)
      hist_wb$rcp <- "Historical"
      
      hist <- rbind(hist, hist_wb)
    }
    hist <- as.data.frame(hist)
    
    write.csv(hist, here('Data',SiteID_FileName,paste('Historical_Streamflow_Projections',SiteID_FileName, "1960_2005_point.csv", sep='_')), row.names=FALSE)
  } else{
    hist <- read.csv(here('Data',SiteID_FileName,paste('Historical_Streamflow_Projections',SiteID_FileName, "1960_2005_point.csv", sep='_')))
    hist$date <- as.Date(hist$date)
  }
} else{
  if(!file.exists(here('Data',SiteID_FileName,paste('Historical_Streamflow_Projections',SiteID_FileName, "1960_2005_area.csv", sep='_')))){
    maca_hist <- get_maca_hist_area(aoi, 1960, 2005, SiteID_FileName, model_names)
    maca_hist$gcm <- maca_hist$projection; maca_hist$projection <- paste0(maca_hist$projection, ".Hist")
    
    # pull only the four identified models
    #maca_hist <- filter(grepl(paste(sapply(strsplit(model_names, "\\."), `[`, 1), collapse = "|"), gcm))
    
    # run streamflow model
    hist <- NULL
    for(mod in gcm_list){
      hist_clim <- maca_hist %>% filter(grepl(strsplit(mod, "\\.")[[1]][1], gcm))
      
      hist_wb <- WB(hist_clim, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
                    shade.coeff, jtemp,SWC.Max, 1, 1, 0, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, "")
      hist_wb <- Drain(hist_wb, q0, s0, v0, qa, qb, sa, sb, va, vb)
      hist_wb <- hist_wb %>% select(projection, date, adj_runoff, Quick, Slow, Very_Slow, total, gcm)
      hist_wb <- hist_wb %>% rename(quick=Quick, slow=Slow, veryslow=Very_Slow)
      hist_wb$rcp <- "Historical"
      
      hist <- rbind(hist, hist_wb)
    }
    hist <- as.data.frame(hist)
    
    write.csv(hist, here('Data',SiteID_FileName,paste('Historical_Streamflow_Projections',SiteID_FileName, "1960_2005_area.csv", sep='_')), row.names=FALSE)
  } else{
    hist <- read.csv(here('Data',SiteID_FileName,paste('Historical_Streamflow_Projections',SiteID_FileName, "1960_2005_area.csv", sep='_')))
    hist$date <- as.Date(hist$date)
  }
}


# Create dataframe with only the MACA data
futures$rcp <- as.character(futures$rcp)
model_df <- rbind(hist, futures) 

# Add useful data
model_df$yr<-as.numeric(format(model_df$date,"%Y")); model_df$mo<-format(model_df$date,"%m"); model_df$yr_mo<-format(model_df$date,"%Y-%m"); model_df$day <- factor(yday(model_df$date), levels = unique(yday(model_df$date)))
model_df$water_year <- sapply(model_df$date, get_water_year); 
model_df <- model_df %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
model_df<-model_df[,c("date", "projection", "gcm", "rcp", "yr", "mo", "yr_mo", "water_year", "adj_runoff", "quick", "slow", "veryslow", "total")]
model_df$Period<-ifelse(model_df$yr<=2005,"Historical",ifelse (model_df$yr>=2006 & model_df$yr<=2050,"Early",
                                                               ifelse (model_df$yr>=2051 & model_df$yr<=2070,"Middle", ifelse (model_df$yr>=2071, "Late","NA"))))


# Aggregate data to annual
model_annual_df <- as.data.frame(model_df %>% group_by(projection, yr) %>%
                             dplyr::summarize(gcm=first(gcm), rcp=first(rcp), water_year=first(water_year),
                                              adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                                              slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE), total = sum(total, na.rm = TRUE), 
                                              Period=first(Period)))
model_annual_df$date<-as.Date(paste(model_annual_df$yr,"-01", "-01",sep=""))

# Aggregate data to monthly 
model_monthly_df <- as.data.frame(model_df %>% group_by(projection, yr_mo) %>%
                              dplyr::summarize(gcm=first(gcm), rcp=first(rcp), water_year=first(water_year),
                                               adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                                               slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE) ,total = sum(total, na.rm = TRUE),
                                               Period=first(Period)))
model_monthly_df$date<-as.Date(paste(model_monthly_df$yr_mo,"-01",sep=""))




#######################################################################
#######################################################################
### PLOTS WITH ALL MODELS (MACA HISTORICAL) ###

if(make_plots){
  # Daily streamflow projections for all models
  jpeg(file=paste0(outLocationPathFuture, "/Daily_Streamflow_MACA.jpg"), width=1500, height=500)
  plot<- ggplot() + geom_line(data=model_df, aes(x = date, y = total, color=rcp)) + 
    facet_wrap(~gcm, ncol=7)+ labs(title="Daily Streamflow", y="Streamflow (mm)", x="Year", color="RCP") + nps_theme() +
    scale_color_manual(values = c('Historical'='black', "45"="orange", "85"='red'), labels = c("Historical"="Historical", "45"="RCP 4.5", "85"="RCP 8.5"))
  print(plot)
  dev.off()
  
  # Monthly streamflow projections for all models
  jpeg(file=paste0(outLocationPathFuture, "/Monthly_Streamflow_MACA.jpg"), width=1500, height=500)
  plot<- ggplot() + geom_line(data=model_monthly_df, aes(x = date, y = total, color = rcp)) + 
    facet_wrap(~gcm, ncol=7)+ labs(title="Monthly Streamflow", y="Streamflow (mm)", x="Year", color="RCP") + nps_theme() +
    scale_color_manual(values = c('Historical'='black', "45"="orange", "85"='red'), labels = c("Historical"="Historical", "45"="RCP 4.5", "85"="RCP 8.5"))
  print(plot)
  dev.off()
  
  # Annual streamflow projections for all models
  jpeg(file=paste0(outLocationPathFuture, "/Annual_Streamflow_MACA.jpg"), width=1500, height=500)
  plot<- ggplot() + geom_line(data=model_annual_df, aes(x = date, y = total, color = rcp)) + 
    facet_wrap(~gcm, ncol=7)+ labs(title="Annual Streamflow", y="Streamflow (mm)", x="Year", color="RCP") + nps_theme() +
    scale_color_manual(values = c('Historical'='black', "45"="orange", "85"='red'), labels = c("Historical"="Historical", "45"="RCP 4.5", "85"="RCP 8.5"))
  print(plot)
  dev.off()
  

  # Time series of daily, monthly, annual streamflow
  plot_daily <- ggplot() + labs(title='Daily Modeled Streamflow', x='Year', y='Streamflow [mm]', color='Model') + nps_theme() + 
    geom_line(data=(model_df %>% filter(rcp!='Historical')), aes(x=date, y=total), col='gray', alpha = 0.7) + 
    geom_line(data=(model_df %>% filter(rcp=='Historical')), aes(x=date, y=total), col='gray',alpha=0.7) +
    geom_line(data=model_df %>% filter(rcp!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total), col='black', alpha=1, linewidth=1.5) + 
    geom_line(data=model_df %>% filter(rcp=='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total), col='blue', alpha=1, linewidth=1.5) +
    guides(alpha = "none") + theme(legend.position="none")
  plot_mon <-  ggplot() + labs(title='Monthly Modeled Streamflow', x='Year', y='Streamflow [mm]', color='Model') + nps_theme() + 
    geom_line(data=(model_monthly_df %>% filter(projection!='Historical')), aes(x=date, y=total), col='gray', alpha = 0.7) + 
    geom_line(data=(model_monthly_df %>% filter(projection=='Historical')), aes(x=date, y=total), col='gray', alpha=1, linewidth=1.5) +
    geom_line(data=model_monthly_df %>% filter(projection!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total), col='black',alpha=1, linewidth=1.5) + 
    geom_line(data=model_monthly_df %>% filter(rcp=='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total), col='blue', alpha=1, linewidth=1.5) +
    guides(alpha = "none") + theme(legend.position="none")
  plot_ann <-  ggplot() + labs(title='Annual Modeled Streamflow', x='Year', y='Streamflow [mm]', color='Model') + nps_theme() + 
    geom_line(data=(model_annual_df %>% filter(projection!='Historical')), aes(x=date, y=total, col='Individual Model Projections'), alpha = 0.7) + 
    geom_line(data=(model_annual_df %>% filter(projection=='Historical')), aes(x=date, y=total, col='Mean Model Projections'), alpha=1, linewidth=1.5) +
    geom_line(data=model_annual_df %>% filter(projection!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total, col='Mean Model Projections'), alpha=1, linewidth=1.5) + 
    geom_line(data=model_annual_df %>% filter(rcp=='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total, col='Historical'), alpha=1, linewidth=1.5) +
    scale_color_manual(values = c('Historical'='blue', "Individual Model Projections" = "gray", 'Mean Model Projections'="black"),
                       labels= c(expression('Historical', "Individual\nModel Projections"), expression('Mean Model\nProjections'))) + 
    guides(alpha = "none") 
  jpeg(file=paste0(outLocationPathFuture, "/Streamflow_Projections_Time_Series_MACA.jpg"), width=2000, height=600)
  grid.arrange(plot_daily, plot_mon, plot_ann, ncol = 3, widths=c(1,1,1.3))
  dev.off()
  
  
  # Climate future scatterplot: annual magnitude vs daily standard deviation
  delta_plot_ann <- (model_annual_df %>% filter(rcp!="Historical")) %>% 
    left_join((model_annual_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(mean_total=mean(total, na.rm=TRUE)), by = "gcm") %>%
    mutate(delta_annual_mm=total-mean_total)
  delta_plot_sd <- (model_df %>% filter(rcp!="Historical")) %>% group_by(gcm, rcp, yr) %>% summarize(Period=first(Period), fut_sd=sd(total, na.rm=TRUE)) %>%
    left_join((model_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(hist_sd=sd(total, na.rm=TRUE)), by = "gcm") %>% 
    mutate(delta_daily_sd=fut_sd-hist_sd)
  delta_plot <- delta_plot_sd %>% left_join(delta_plot_ann %>% select(gcm, rcp, yr, delta_annual_mm, Period), 
                                         by = c("gcm", "rcp", "yr","Period"))
  delta_plot <- delta_plot %>% filter(gcm!='Historical') %>% group_by(gcm,rcp,Period) %>% dplyr::summarise(delta_daily_sd=mean(delta_daily_sd), delta_annual_mm=mean(delta_annual_mm))
  delta_plot$projection <- paste0(delta_plot$gcm, ".rcp", delta_plot$rcp) 
  
  annual_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                                xintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_annual_mm, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_annual_mm, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_annual_mm, 0.5)))
  sd_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                            yintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_daily_sd, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_daily_sd, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_daily_sd, 0.5)))
  annual_zero <- data.frame(Period = c("Early", "Middle", "Late"), xintercept=c(0,0,0)); sd_zero <- data.frame(Period = c("Early", "Middle", "Late"), yintercept = c(0,0,0))
  
  # plot with RCPs
  jpeg(file=paste0(outLocationPathFuture, "/Streamflow_Climate_Future_Scatter_RCP.jpg"), width=1000, height=400)
  plot <- ggplot(data=delta_plot, aes(x=delta_annual_mm,y=delta_daily_sd,color=rcp)) + geom_point() +
    geom_text_repel(aes(label = gcm), color = 'black', max.overlaps=Inf) +
    geom_hline(data = sd_quantile, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_quantile, aes(xintercept = xintercept), color = "black") +
    #geom_hline(data = sd_zero, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_zero, aes(xintercept = xintercept), color = "black") +
    facet_wrap(~factor(Period, levels = c('Early', 'Middle', 'Late'))) + labs(title=paste('Changes in streamflow at',SiteID), x='Change in annual magnitude [mm]', y='Change in daily standard deviation [mm]', color='RCP') + 
    scale_color_manual(values = c("45" = "orange", "85" = "red")) + nps_theme()
  print(plot)
  dev.off()
  
  # plot with selected climate futures
  jpeg(file=paste0(outLocationPathFuture, "/Streamflow_Climate_Future.jpg"), width=1500, height=600)
  plot <- ggplot(data=delta_plot) + 
    geom_point(data=delta_plot %>% filter(!(projection %in% model_names)), aes(x=delta_annual_mm,y=delta_daily_sd,color='Other')) +
    geom_point(data=delta_plot %>% filter(projection %in% model_names), aes(x=delta_annual_mm,y=delta_daily_sd,color=projection)) + 
    geom_text_repel(data=delta_plot %>% filter(!(projection %in% model_names)), aes(x=delta_annual_mm,y=delta_daily_sd, label=projection), color = 'black', size=2, max.overlaps=Inf) +
    geom_text_repel(data=delta_plot %>% filter(projection %in% model_names), aes(x=delta_annual_mm,y=delta_daily_sd, label=projection, color=projection), size=4,  fontface = "bold", max.overlaps=Inf) + 
    geom_hline(data = sd_quantile, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_quantile, aes(xintercept = xintercept), color = "black") +
    #geom_hline(data = sd_zero, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_zero, aes(xintercept = xintercept), color = "black") +
    facet_wrap(~factor(Period, levels = c('Early', 'Middle', 'Late'))) + labs(title=paste('Changes in streamflow at',SiteID), x='Change in annual magnitude [mm]', y='Change in daily standard deviation [mm]', color='RCP') + 
    scale_color_manual(values = c("Other" = "black", setNames(color_names, model_names)), labels=c("Other" = "black", setNames(scenario_names, model_names))) + nps_theme()
  print(plot)
  dev.off()
}


#######################################################################
#######################################################################
### PLOTS WITH MACA HISTORICAL ###


### Heatmap - daily
model_df <- model_df %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))

plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
  plot_list[[i]] <- ggplot(analysis_df, aes(factor(water_day), water_year, fill=total)) + geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), transform='log', breaks=c(min(model_df$total), 0.01, 0.1, 10, 100, max(model_df$total)), 
                         labels=c(sprintf('%.3f',  min(model_df$total)),'.01','0.1','10', '100', sprintf('%.0f', max(model_df$total))),
                         limits = c(min(model_df$total, na.rm = TRUE), max(model_df$total, na.rm = TRUE))) +
    labs(title = paste0(scenario, " Daily Streamflow (", proj, ")"), x='Month', y='Water Year', fill = "Streamflow [mm]") +
    scale_x_discrete(breaks=c("1", "32", "63", "94", "120", "151", "181", "212", "243", "274", "305", "335"), 
                     labels = c("1" = "October", "32" = "November", "63" = "December", "94" = "January", "120" = "February", "151" = "March", 
                                "181" = "April", "212" = "May", "243" = "June", "274" = "July", "305" = "August", "335" = "September")) +
    theme(axis.text.x = element_text(angle = 90))
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_MACA_Daily_Heatmap.jpg"), width=200*num_models, height=150*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()


### Heatmap - monthly 
# order by water year?
model_monthly_df$month <- month(model_monthly_df$date); model_monthly_df$year <- year(model_monthly_df$date)

plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- model_monthly_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
  
  plot_list[[i]] <- ggplot(analysis_df, aes(factor(month), year, fill = total)) + geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), transform='log', breaks=c(0.1, 1, 10, 100, 1000, max(model_monthly_df$total)), 
                         labels=c('0.1', '1', '10', '100', '1000', sprintf('%.0f', max(model_monthly_df$total))),
                         limits = c(min(model_monthly_df$total, na.rm = TRUE), max(model_monthly_df$total, na.rm = TRUE))) + 
    labs(title = paste0(scenario, " Monthly Streamflow (", proj, ")"), x='Month', y='Water Year', fill = "Streamflow [mm]") +
    scale_x_discrete(labels = c("1" = "January", "2" = "February", "3" = "March", 
                                "4" = "April", "5" = "May", "6" = "June", 
                                "7" = "July", "8" = "August", "9" = "September", 
                                "10" = "October", "11" = "November", "12" = "December")) +
    nps_theme() + theme(axis.text.x = element_text(angle = 90))
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_MACA_Monthly_Heatmap.jpg"), width=400*num_models, height=300*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()


#######################################################################
### PLOT ECDF HISTORICAL VS FUTURE ###
### KS test to evaluate distributions ###
# daily
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
    
    daily_ks <- ks.test((analysis_df %>% filter(rcp=='Historical'))$total, (analysis_df %>% filter(rcp!='Historical'))$total)
    if(daily_ks$p.value <= 0.05){label <- sprintf('Difference: Significant \n p-value: %.2f', daily_ks$p.value)
    }else{label <- sprintf('Difference: Not significant \n p-value: %.2f', daily_ks$p.value)}
    
    plot_list[[i]] <- ggplot(data=analysis_df) + 
      stat_ecdf(aes(x=total, color=factor(rcp)), linewidth=1) + 
      labs(title=paste(scenario, 'Daily ECDF'), x='Streamflow (mm)', y='Cumulative Frequency', color='') +
      annotate("text", x = max(analysis_df$total), y = .2, label = label, color = "black", hjust = 1, vjust = 1) + 
      scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i]), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
      scale_x_log10() + nps_theme() + theme(legend.position = "bottom")
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_ECDF_Daily.jpg"), width=200*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}


# monthly
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_monthly_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
    
    daily_ks <- ks.test((analysis_df %>% filter(rcp=='Historical'))$total, (analysis_df %>% filter(rcp!='Historical'))$total)
    if(daily_ks$p.value <= 0.05){label <- sprintf('Difference: Significant \n p-value: %.2f', daily_ks$p.value)
    }else{label <- sprintf('Difference: Not significant \n p-value: %.2f', daily_ks$p.value)}
    
    plot_list[[i]] <- ggplot(data=analysis_df) + stat_ecdf(aes(x=total, color=factor(rcp)), linewidth=1) + 
      labs(title=paste(scenario, 'Monthly ECDF'), x='Streamflow (mm)', y='Cumulative Frequency', color='') +
      annotate("text", x = max(analysis_df$total), y = .2, label = label, color = "black", hjust = 1, vjust = 1) + 
      scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i]), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
      scale_x_log10() + nps_theme() + theme(legend.position = "bottom")
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_ECDF_Monthly.jpg"), width=200*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}

# annual
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
    
    daily_ks <- ks.test((analysis_df %>% filter(rcp=='Historical'))$total, (analysis_df %>% filter(rcp!='Historical'))$total)
    if(daily_ks$p.value <= 0.05){label <- sprintf('Difference: Significant \n p-value: %.2f', daily_ks$p.value)
    }else{label <- sprintf('Difference: Not significant \n p-value: %.2f', daily_ks$p.value)}
    
    plot_list[[i]] <- ggplot(data=analysis_df) + stat_ecdf(aes(x=total, color=factor(rcp)), linewidth=1) + 
      labs(title=paste(scenario, 'Annual ECDF'), x='Streamflow (mm)', y='Cumulative Frequency', color='') +
      annotate("text", x = max(analysis_df$total), y = .2, label = label, color = "black", hjust = 1, vjust = 1) + 
      scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i]), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
      scale_x_log10() + nps_theme() + theme(legend.position = "bottom")
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_ECDF_Annual.jpg"), width=200*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}


#######################################################################
### BOXPLOTS ###

### Change in annual streamflow magnitude ###
delta_plot <- (model_annual_df %>% filter(rcp!="Historical")) %>% 
  left_join((model_annual_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(mean_total=mean(total, na.rm=TRUE)), by = "gcm") %>%
  mutate(delta_annual_mm=total-mean_total)

# All models by RCP
plot <- ggplot() + geom_boxplot(data=delta_plot, aes(x=projection, y=delta_annual_mm, fill=rcp)) + 
  labs(title='Change in annual flow magnitude', y='Change in annual flow [mm]', x='Projection', fill='RCP') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c('45' = 'orange', '85' = 'red')) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Flow_Change_RCP.jpg"), width=1200, height=500)
print(plot); dev.off()

# All models by climate future
plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta_annual_mm, fill=projection)) + 
  geom_boxplot(data=(delta_plot %>% filter(!projection %in% model_names)), aes(x=projection, y=delta_annual_mm, fill='Other')) + 
  labs(title='Change in annual flow magnitude', y='Change in annual flow [mm]', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Flow_Change.jpg"), width=1200, height=500)
print(plot); dev.off()

# Just climate futures
plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta_annual_mm, fill=projection)) + 
  labs(title='Change in annual flow magnitude', y='Change in annual flow [mm]', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Flow_Change_Futures.jpg"), width=1200, height=500)
print(plot); dev.off()


### Change in daily streamflow standard deviation ###
delta_plot_sd <- (model_df %>% filter(rcp!="Historical")) %>% group_by(projection, yr) %>% summarize(gcm=first(gcm), rcp=first(rcp), Period=first(Period), fut_sd=sd(total, na.rm=TRUE)) %>%
  left_join((model_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(hist_sd=sd(total, na.rm=TRUE)), by = "gcm") %>% 
  mutate(delta_daily_sd=fut_sd-hist_sd)

# All models by RCP
plot <- ggplot() + geom_boxplot(data=delta_plot_sd, aes(x=projection, y=delta_daily_sd, fill=rcp)) + 
  labs(title='Change in daily flow standard deviation', y='Change in standard deviation [mm]', x='Projection', fill='RCP') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c('45' = 'orange', '85' = 'red')) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "SD_Change_RCP.jpg"), width=1200, height=500)
print(plot); dev.off()

# All models by climate future
plot <- ggplot() + geom_boxplot(data=(delta_plot_sd %>% filter(projection %in% model_names)), aes(x=projection, y=delta_daily_sd, fill=projection)) + 
  geom_boxplot(data=(delta_plot_sd %>% filter(!projection %in% model_names)), aes(x=projection, y=delta_daily_sd, fill='Other')) + 
  labs(title='Change in daily flow standard deviation', y='Change in standard deviation [mm]', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "SD_Change.jpg"), width=1200, height=500)
print(plot); dev.off()

# Just climate futures
plot <- ggplot() + geom_boxplot(data=(delta_plot_sd %>% filter(projection %in% model_names)), aes(x=projection, y=delta_daily_sd, fill=projection)) + 
  labs(title='Change in daily flow standard deviation', y='Change in standard deviation [mm]', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "SD_Change_Futures.jpg"), width=1200, height=500)
print(plot); dev.off()


### Change in high flows (95th percentile) ###
high_flow_q = 0.95
high_threshold <- model_df %>% filter(rcp == "Historical") %>% group_by(gcm) %>% summarize(threshold = quantile(total, high_flow_q, na.rm = TRUE))
delta_plot <- model_df %>% filter(rcp == "Historical") %>% inner_join(high_threshold, by = "gcm") %>% 
  group_by(gcm, yr) %>% summarize(threshold=first(threshold), hist_days = sum(total > threshold, na.rm = TRUE), .groups = "drop") %>% group_by(gcm) %>% summarize(threshold=first(threshold), hist_days = mean(hist_days)) %>% 
  right_join(model_df %>% filter(rcp!="Historical"), by='gcm') %>% distinct(projection, yr, gcm, rcp, Period, threshold, total, hist_days) %>%  group_by(projection, yr) %>% summarize(gcm=first(gcm), rcp=first(rcp), Period=first(Period), hist_days=first(hist_days), fut_days = sum(total > threshold, na.rm = TRUE), delta_days=fut_days-hist_days)

# All models by RCP
plot <- ggplot() + geom_boxplot(data=delta_plot, aes(x=projection, y=delta_days, fill=rcp)) + 
  labs(title='Change in days above historic 95th percentile', y='Change in days', x='Projection', fill='RCP') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c('45' = 'orange', '85' = 'red')) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "High_Flow_Change_RCP.jpg"), width=1200, height=500)
print(plot); dev.off()

# All models by climate future
plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta_days, fill=projection)) + 
  geom_boxplot(data=(delta_plot %>% filter(!projection %in% model_names)), aes(x=projection, y=delta_days, fill='Other')) + 
  labs(title='Change in days above historic 95th percentile', y='Change in days', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "High_Flow_Change.jpg"), width=1200, height=500)
print(plot); dev.off()

# Just climate futures
plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta_days, fill=projection)) + 
  labs(title='Change in days above historic 95th percentile', y='Change in days', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "High_Flow_Change_Futures.jpg"), width=1200, height=500)
print(plot); dev.off()


### Change in low flows (5th percentile) ###
low_flow_q = 0.05
low_threshold <- model_df %>% filter(rcp == "Historical") %>% group_by(gcm) %>% summarize(threshold = quantile(total, high_flow_q, na.rm = TRUE))
delta_plot <- model_df %>% filter(rcp == "Historical") %>% inner_join(low_threshold, by = "gcm") %>% 
  group_by(gcm, yr) %>% summarize(threshold=first(threshold), hist_days = sum(total < threshold, na.rm = TRUE), .groups = "drop") %>% group_by(gcm) %>% summarize(threshold=first(threshold), hist_days = mean(hist_days)) %>% 
  right_join(model_df %>% filter(rcp!="Historical"), by='gcm') %>% distinct(projection, yr, gcm, rcp, Period, threshold, total, hist_days) %>%  group_by(projection, yr) %>% summarize(gcm=first(gcm), rcp=first(rcp), Period=first(Period), hist_days=first(hist_days), fut_days = sum(total < threshold, na.rm = TRUE), delta_days=fut_days-hist_days)

# All models by RCP
plot <- ggplot() + geom_boxplot(data=delta_plot, aes(x=projection, y=delta_days, fill=rcp)) + 
  labs(title='Change in days below historic 5th percentile', y='Change in days', x='Projection', fill='RCP') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c('45' = 'orange', '85' = 'red')) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "Low_Flow_Change_RCP.jpg"), width=1200, height=500)
print(plot); dev.off()

# All models by climate future
plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta_days, fill=projection)) + 
  geom_boxplot(data=(delta_plot %>% filter(!projection %in% model_names)), aes(x=projection, y=delta_days, fill='Other')) + 
  labs(title='Change in days below historic 5th percentile', y='Change in days', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "Low_Flow_Change.jpg"), width=1200, height=500)
print(plot); dev.off()

# Just climate futures
plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta_days, fill=projection)) + 
  labs(title='Change in days below historic 5th percentile', y='Change in days', x='Projection', fill='Future') + 
  facet_wrap(~factor(Period, levels=c("Early","Middle","Late"))) + scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
  nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
jpeg(file=paste0(outLocationPathFuture, "/", "Low_Flow_Change_Futures.jpg"), width=1200, height=500)
print(plot); dev.off()


#######################################################################
### PLOT STREAMFLOW TRENDS AND METRICS ###
# edit width of jpg figure depending on model_names

# Plot daily streamflow
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | (rcp=='Historical')))
  analysis_df <- analysis_df[order(analysis_df$date), ]
  plot_list[[i]] <- ggplot(analysis_df, aes(x = date, y = total, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Daily Streamflow (mm)", title = paste(scenario, "Daily Modeled Streamflow"), color='') +
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    nps_theme() + theme(legend.position = 'bottom') + scale_y_log10()
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Daily_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()


# Seasonal Mann-Kendall test on monthly streamflow
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- model_monthly_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | (rcp=='Historical')))
  analysis_df <- analysis_df[order(analysis_df$date), ]
  
  mod_mk <- SeasonalMannKendall(ts(analysis_df$total, start=c(startY, 1), frequency=12))
  mod_sens <- sens.slope(analysis_df$total[!is.na(analysis_df$total)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_list[[i]] <- ggplot(analysis_df, aes(x = as.yearmon(yr_mo), y = total, color=factor(rcp))) + 
    geom_line(na.rm=TRUE) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1) +
    labs(x = "Water Year", y = "Monthly Streamflow (mm)", title = paste(scenario, "Monthly Modeled Streamflow"), color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(as.yearmon(monthly_df$yr_mo)), y = max(monthly_df$total), label = label, color = "black", hjust = 1, vjust = 1)
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Monthly_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()


# Annual streamflow volume
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- model_annual_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
  mod_mk <- MannKendall(analysis_df$total)
  mod_sens <- sens.slope(analysis_df$total[!is.na(analysis_df$total)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_list[[i]] <- ggplot(analysis_df, aes(x = yr, y = total, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Annual Streamflow (mm)", title = paste(scenario, "Annual Modeled Streamflow"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$total), label = label, color = "black", hjust = 1, vjust = 1) 
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Annual_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()


# Streamflow distribution
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- model_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
    dplyr::summarize(gcm=first(gcm), rcp=first(rcp), value = mean(total, na.rm = TRUE), .groups = 'drop') %>%
    mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
  
  plot_list[[i]] <- ggplot() + 
    geom_smooth(data = (analysis_df %>% filter(!grepl(strsplit(proj, "\\.")[[1]][1], gcm) & !(grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
                method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
    geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp))), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
                method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
    geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp))), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
                method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
    labs(title=paste("Smoothed Daily Streamflow"), y=paste('Streamflow [mm]'), x=NULL, color='Model') + nps_theme() + 
    scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
    guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
}
jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_Flow_Trends.jpg"), width=300*num_models, height=150*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()



# High flows (above 95%)
high_flow_q = 0.95

plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  # Calculate days above 95th percentile
  high_flow_mm = quantile((model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp)))$total, 0.95, na.rm=TRUE)
  analysis_df <- as.data.frame(model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp))) %>% 
                                 mutate(high_flow = ifelse(total >= high_flow_mm, 1, 0)) %>%
                                 group_by(yr) %>% dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), days = sum(high_flow)))
  
  # Plot
  mod_mk <- MannKendall(analysis_df$days)
  mod_sens <- sens.slope(analysis_df$days[!is.na(analysis_df$days)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_list[[i]] <- ggplot(analysis_df, aes(x = yr, y = days, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Days per year", title = paste(scenario, "Days Above Historical 95th Percentile"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_High_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()



# Low flows (below 5%)
low_flow_q = 0.05

plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  # Calculate days below 5th percentile
  low_flow_mm = quantile((model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp)))$total, 0.05, na.rm=TRUE)
  analysis_df <- as.data.frame(model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp))) %>% 
                                 mutate(low_flow = ifelse(total <= low_flow_mm, 1, 0)) %>%
                                 group_by(yr) %>% dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), days = sum(low_flow)))
  
  # Plot
  mod_mk <- MannKendall(analysis_df$days)
  mod_sens <- sens.slope(analysis_df$days[!is.na(analysis_df$days)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_list[[i]] <- ggplot(analysis_df, aes(x = yr, y = days, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Days per year", title = paste(scenario, "Days Below Historical 5th Percentile"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Low_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()


# 50% flow date
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  # Pull projections and add water day
  analysis_df <- model_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp))) %>% 
    group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days")))) %>% filter(water_year < 2099 & water_year > 1960 & water_year!=2005 & water_year!=2006)
  
  # Calculate center of timing (CT)
  analysis_df <- as.data.frame(analysis_df) %>% mutate(tq = total*water_day) %>%
    group_by(water_year) %>% dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), ct = sum(tq)/sum(total))
  
  # Plot
  mod_mk <- MannKendall(analysis_df$ct)
  mod_sens <- sens.slope(analysis_df$ct[!is.na(analysis_df$ct)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_list[[i]] <- ggplot(analysis_df, aes(x = water_year, y = ct, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Days after October 1", title = paste(scenario, "50% Flow Date"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$ct), label = label, color = "black", hjust = 1, vjust = 1)
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_50th_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()





# Calculate Q7
hist_q7 <- model_df %>% mutate(Q7 = rollmean(total, 7, align="right", fill=NA)) %>%
  group_by(projection, water_year) %>% dplyr::summarize(gcm=first(gcm), rcp=first(rcp), min_q7 = ifelse(all(is.na(Q7)), NA, min(Q7, na.rm = TRUE)), 
                                                        max_q7 = ifelse(all(is.na(Q7)), NA, max(Q7, na.rm = TRUE)), avg_q7 = ifelse(all(is.na(Q7)), NA, mean(Q7, na.rm = TRUE)))

# Q7 min
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- hist_q7 %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
  mod_mk <- MannKendall(analysis_df$min_q7)
  mod_sens <- sens.slope(analysis_df$min_q7[!is.na(analysis_df$min_q7)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_list[[i]] <- ggplot(analysis_df, aes(x = water_year, y = min_q7, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Streamflow [mm]", title = paste(scenario, "Min. 7 Day Flow (Q7 Min)"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$min_q7), label = label, color = "black", hjust = 1, vjust = 1) 
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Q7Min_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()


# Q7 max
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- hist_q7 %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
  mod_mk <- MannKendall(analysis_df$max_q7)
  mod_sens <- sens.slope(analysis_df$max_q7[!is.na(analysis_df$max_q7)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_list[[i]] <- ggplot(analysis_df, aes(x = water_year, y = max_q7, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Streamflow [mm]", title = paste(scenario, "Max. 7 Day Flow (Q7 Max)"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$max_q7), label = label, color = "black", hjust = 1, vjust = 1) 
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Q7Max_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()




#######################################################################
### Customizable plots: metric, season, etc ###
# EDIT TO PLOT SPECIFIC FUTURE SCENARIOS

# define flow level in cfs, convert to mm
flow_level <- 3200 * 28316847*86400/(2590000000000*sqmi) 
# identify months of interest (numerical values)
mos <- c("02", "03") 
# is comparison above or below threshold?
comparison = 'above'


# calculate number of days above/below threshold
if (tolower(comparison)=='below') {
  threshold <- as.data.frame(model_df %>% filter(mo %in% mos) %>% mutate(flow = ifelse(total <= flow_level, 1, 0)) %>%
                               group_by(water_year, projection) %>% dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), days = sum(flow)))# %>%filter(!is.na(days)))
} else if (tolower(comparison)=='above'){
  threshold <- as.data.frame(model_df %>% filter(mo %in% mos) %>% mutate(flow = ifelse(total >= flow_level, 1, 0)) %>%
                               group_by(water_year, projection) %>% dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), days = sum(flow)))# %>%filter(!is.na(days)))
}

# plot
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- threshold %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))
  
  mod_mk <- MannKendall(analysis_df$days)
  mod_sens <- sens.slope(analysis_df$days[!is.na(analysis_df$days)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  
  plot_list[[i]] <- ggplot(analysis_df, aes(x = water_year, y = days, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Days", title = paste('Days', comparison, round(flow_level), 'mm, ', month.abb[as.numeric(mos)][1], '-', month.abb[as.numeric(mos)][length(mos)]), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Days_", comparison, "_", round(flow_level), "mm_", month.abb[as.numeric(mos)][1], "_", month.abb[as.numeric(mos)][length(mos)],".jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()




#######################################################################
### Non-streamflow metrics ###

# Runoff efficiency
if(point_location){maca_hist <- get_maca_hist_point(lat, lon, 1960, 2005, SiteID_FileName, model_names)
} else{maca_hist <- get_maca_hist_area(aoi, 1960, 2005, SiteID_FileName, model_names)}
if(point_location) {future_climate <- get_maca_point(lat, lon, 2006, 2100, SiteID_FileName, gcm_list)
} else future_climate <- get_maca_area(aoi, 2006, 2100, SiteID_FileName, gcm_list)

runoff_efficiency <- model_df %>% group_by(water_year, projection) %>% dplyr::summarize(gcm=first(gcm), rcp=first(rcp), ann_runoff=sum(adj_runoff))
hist_p <- maca_hist %>% mutate(water_year = sapply(date, get_water_year)) %>% group_by(water_year, projection) %>% dplyr::summarize(gcm=first(projection), ann_p=sum(pr)) %>% mutate(projection=paste0(gcm, ".Hist"), rcp="Historical")
fut_p <- future_climate %>% mutate(water_year = sapply(date, get_water_year)) %>% group_by(water_year, projection) %>% dplyr::summarize(gcm=first(strsplit(projection, "\\.")[[1]][1]), rcp=first(gsub("\\D", "", strsplit(projection, "\\.")[[1]][2])), ann_p=sum(pr))
runoff_efficiency <- runoff_efficiency %>% full_join(bind_rows(hist_p, fut_p), by = c('water_year', 'projection', 'gcm', 'rcp')) %>% filter(water_year!=2100 & water_year!=1960 & water_year!=2005 & water_year!=2006)
runoff_efficiency$efficiency <- (runoff_efficiency$ann_runoff / runoff_efficiency$ann_p) * 100

# Plot
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- runoff_efficiency %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & (grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp))) 
  mod_mk <- MannKendall(analysis_df$efficiency)
  mod_sens <- sens.slope(analysis_df$efficiency[!is.na(analysis_df$efficiency)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_mod <- ggplot(analysis_df, aes(x = water_year, y = efficiency, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Runoff Efficiency (%)", title = paste(scenario, "Runoff Efficiency"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + 
    scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c("Historical"="Historical", '45'=proj, '85'=proj)) + 
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$efficiency), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_mod)
  plot_list[[i]] <- plot_mod
}
jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_RunoffEffiency_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()







#######################################################################
#######################################################################
### PLOTS WITH GRIDMET HISTORICAL ###
# not including to provide more accurate comparisons (using MACA historical)
# 
# ### Heatmap - daily
# daily_df <- daily_df %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
# 
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   analysis_df <- daily_df %>% filter(projection=='Historical' | projection==proj)
#   plot <- ggplot(analysis_df, aes(factor(water_day), water_year, fill=total)) + geom_tile() +
#     scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), transform='log', breaks=c(min(daily_df$total), 0.01, 0.1, 10, 100, max(daily_df$total)), 
#                          labels=c(sprintf('%.3f',  min(daily_df$total)),'.01','0.1','10', '100', sprintf('%.0f', max(daily_df$total))),
#                          limits = c(min(daily_df$total, na.rm = TRUE), max(daily_df$total, na.rm = TRUE))) +
#     labs(title = paste0(scenario, " Daily Streamflow (", proj, ")"), x='Month', y='Water Year', fill = "Streamflow [mm]") +
#     scale_x_discrete(breaks=c("1", "32", "63", "94", "120", "151", "181", "212", "243", "274", "305", "335"), 
#                      labels = c("1" = "October", "32" = "November", "63" = "December", "94" = "January", "120" = "February", "151" = "March", 
#                                 "181" = "April", "212" = "May", "243" = "June", "274" = "July", "305" = "August", "335" = "September")) +
#     theme(axis.text.x = element_text(angle = 90))
#   plot_list[[i]] <- plot
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Daily_Heatmap.jpg"), width=200*num_models, height=150*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# ### Heatmap - monthly 
# # order by water year?
# monthly_df$month <- month(monthly_df$date); monthly_df$year <- year(monthly_df$date)
# 
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#
#   analysis_df <- monthly_df %>% filter(projection=='Historical' | projection==proj)
#   plot <- ggplot(analysis_df, aes(factor(month), year, fill = total)) + geom_tile() +
#     #scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
#     scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), transform='log', breaks=c(0.1, 1, 10, 100, 1000, max(monthly_df$total)), 
#                          labels=c('0.1', '1', '10', '100', '1000', sprintf('%.0f', max(monthly_df$total))),
#                          limits = c(min(monthly_df$total, na.rm = TRUE), max(monthly_df$total, na.rm = TRUE))) + 
#     labs(title = paste0(scenario, " Monthly Streamflow (", proj, ")"), x='Month', y='Water Year', fill = "Streamflow [mm]") +
#     scale_x_discrete(labels = c("1" = "January", "2" = "February", "3" = "March", 
#                                 "4" = "April", "5" = "May", "6" = "June", 
#                                 "7" = "July", "8" = "August", "9" = "September", 
#                                 "10" = "October", "11" = "November", "12" = "December")) +
#     nps_theme() + theme(axis.text.x = element_text(angle = 90))
#   plot_list[[i]] <- plot
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Monthly_Heatmap.jpg"), width=200*num_models, height=150*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# #######################################################################
# ### PLOT ECDF HISTORICAL VS FUTURE ###
# ### KS test to evaluate distributions ###
# # daily
# if(make_plots){
#   jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_ECDF_Daily.jpg"), width=500, height=300)
#   plot <- ggplot(data=daily_df %>% filter(projection %in% c(model_names,"Historical"))) + 
#     stat_ecdf(aes(x=total, color=factor(projection)), linewidth=1) + 
#     labs(title='Daily ECDF', x='Streamflow (mm)', y='Cumulative Frequency', color='') +
#     scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names)), labels = c("Historical"="Historical", setNames(scenario_names, model_names))) + 
#     scale_x_log10() + nps_theme()
#   print(plot); dev.off()
# }
# 
# # monthly
# if(make_plots){
#   jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_ECDF_Monthly.jpg"), width=500, height=300)
#   plot <- ggplot(data=monthly_df %>% filter(projection %in% c(model_names,"Historical"))) + 
#     stat_ecdf(aes(x=total, color=factor(projection)), linewidth=1) + 
#     labs(title='Monthly ECDF', x='Streamflow (mm)', y='Cumulative Frequency', color='') +
#     scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names)), labels = c("Historical"="Historical", setNames(scenario_names, model_names))) + 
#     scale_x_log10() + nps_theme()
#   print(plot); dev.off()
# }
# 
# # annual
# if(make_plots){
#   jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_ECDF_Annual.jpg"), width=500, height=300)
#   plot <- ggplot(data=annual_df %>% filter(projection %in% c(model_names,"Historical"))) + 
#     stat_ecdf(aes(x=total, color=factor(projection)), linewidth=1) + 
#     labs(title='Annual ECDF', x='Streamflow (mm)', y='Cumulative Frequency', color='') +
#     scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names)), labels = c("Historical"="Historical", setNames(scenario_names, model_names))) + 
#     scale_x_log10() + nps_theme()
#   print(plot); dev.off()
# }
# 
# 
# 
# #######################################################################
# ### PLOT STREAMFLOW TRENDS AND METRICS ###
# # edit width of jpg figure depending on model_names
# 
# # Plot daily streamflow
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
# 
#   analysis_df <- daily_df %>% filter(projection=='Historical' | projection==proj)
#   analysis_df <- analysis_df[order(analysis_df$date), ]
#   plot_mod <- ggplot(analysis_df, aes(x = date, y = total, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Daily Streamflow (mm)", title = paste(scenario, "Daily Modeled Streamflow"), color='') +
#     scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) + 
#     nps_theme() + theme(legend.position = 'bottom') + scale_y_log10()
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Daily_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# # Seasonal Mann-Kendall test on monthly streamflow
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   analysis_df <- monthly_df %>% filter(projection=='Historical' | projection==proj)
#   analysis_df <- analysis_df[order(analysis_df$date), ]
#   
#   mod_mk <- SeasonalMannKendall(ts(analysis_df$total, start=c(startY, 1), frequency=12))
#   mod_sens <- sens.slope(analysis_df$total[!is.na(analysis_df$total)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = as.yearmon(yr_mo), y = total, color=factor(projection))) + 
#     geom_line(na.rm=TRUE) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1) +
#     labs(x = "Water Year", y = "Monthly Streamflow (mm)", title = paste(scenario, "Monthly Modeled Streamflow"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') +
#     scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(as.yearmon(monthly_df$yr_mo)), y = max(monthly_df$total), label = label, color = "black", hjust = 1, vjust = 1)
#   print(plot_mod); plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Monthly_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# # Annual streamflow volume
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   analysis_df <- annual_df %>% filter(projection=='Historical' | projection==proj)
#   mod_mk <- MannKendall(analysis_df$total)
#   mod_sens <- sens.slope(analysis_df$total[!is.na(analysis_df$total)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = yr, y = total, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Annual Streamflow (mm)", title = paste(scenario, "Annual Modeled Streamflow"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$yr), y = max(analysis_df$total), label = label, color = "black", hjust = 1, vjust = 1) 
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Annual_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# # Streamflow distribution
# daily_mean_mod <- daily_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
#   dplyr::summarize(value = mean(total, na.rm = TRUE), .groups = 'drop') %>%
#   mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
# plot <- ggplot() + 
#   geom_smooth(data = (daily_mean_mod %>% filter(!projection %in% model_names)), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
#                              method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
#   geom_smooth(data = (daily_mean_mod %>% filter(projection %in% model_names)), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
#                 method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
#   geom_smooth(data = (daily_mean_mod %>% filter(projection == 'Historical')), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
#               method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
#   labs(title=paste("Smoothed Future Average Daily Streamflow"), y=paste('Streamflow [mm]'), x=NULL, color='Model') + nps_theme() + 
#   scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
#   guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
# jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_Flow_Trends.jpg"), width=600, height=300)
# print(plot)
# dev.off()
# 
# 
# # High flows (above 95%)
# high_flow_q = 0.95
# high_flow_mm = quantile((daily_df %>% filter(projection=='Historical'))$total, 0.95, na.rm=TRUE)
# 
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   # Calculate days above 95th percentile
#   analysis_df <- as.data.frame(daily_df %>% filter(projection=='Historical' | projection==proj) %>% 
#                                  mutate(high_flow = ifelse(total >= high_flow_mm, 1, 0)) %>%
#                                  group_by(yr) %>% dplyr::summarize(projection=first(projection), days = sum(high_flow)))
#   
#   # Plot
#   mod_mk <- MannKendall(analysis_df$days)
#   mod_sens <- sens.slope(analysis_df$days[!is.na(analysis_df$days)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = yr, y = days, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Days per year", title = paste(scenario, "Days Above Historical 95th Percentile"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$yr), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_High_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# 
# # Low flows (below 5%)
# low_flow_q = 0.05
# low_flow_mm = quantile(daily_df$total, 0.05, na.rm=TRUE)
# 
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   # Calculate days below 5th percentile
#   analysis_df <- as.data.frame(daily_df %>% filter(projection=='Historical' | projection==proj) %>% 
#                                  mutate(low_flow = ifelse(total <= low_flow_mm, 1, 0)) %>%
#                                  group_by(yr) %>% dplyr::summarize(projection=first(projection), days = sum(low_flow)))
#   
#   # Plot
#   mod_mk <- MannKendall(analysis_df$days)
#   mod_sens <- sens.slope(analysis_df$days[!is.na(analysis_df$days)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = yr, y = days, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Days per year", title = paste(scenario, "Days Below Historical 5th Percentile"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$yr), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Low_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# # 50% flow date
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   # Pull projections and add water day
#   analysis_df <- daily_df %>% filter(projection=='Historical' | projection==proj) %>% group_by(water_year) %>% 
#     mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days")))) %>% filter(water_year < 2099 & water_year > 1979)
#   
#   # Calculate center of timing (CT)
#   analysis_df <- as.data.frame(analysis_df) %>% mutate(tq = total*water_day) %>%
#                                  group_by(water_year) %>% dplyr::summarize(projection=first(projection), ct = sum(tq)/sum(total))
# 
#   # Plot
#   mod_mk <- MannKendall(analysis_df$ct)
#   mod_sens <- sens.slope(analysis_df$ct[!is.na(analysis_df$ct)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = water_year, y = ct, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Days after October 1", title = paste(scenario, "50% Flow Date"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$ct), label = label, color = "black", hjust = 1, vjust = 1)
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_50th_Flow_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# 
# 
# 
# # Calculate Q7
# hist_q7 <- daily_df %>% mutate(Q7 = rollmean(total, 7, align="right", fill=NA)) %>%
#   group_by(projection, water_year) %>% dplyr::summarize(min_q7 = ifelse(all(is.na(Q7)), NA, min(Q7, na.rm = TRUE)), 
#                                            max_q7 = ifelse(all(is.na(Q7)), NA, max(Q7, na.rm = TRUE)), 
#                                            avg_q7 = ifelse(all(is.na(Q7)), NA, mean(Q7, na.rm = TRUE)))
# 
# # Q7 min
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   analysis_df <- hist_q7 %>% filter(projection=='Historical' | projection==proj)
#   mod_mk <- MannKendall(analysis_df$min_q7)
#   mod_sens <- sens.slope(analysis_df$min_q7[!is.na(analysis_df$min_q7)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = water_year, y = min_q7, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Streamflow [mm]", title = paste(scenario, "Min. 7 Day Flow (Q7 Min)"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$min_q7), label = label, color = "black", hjust = 1, vjust = 1) 
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Q7Min_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# # Q7 max
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   analysis_df <- hist_q7 %>% filter(projection=='Historical' | projection==proj)
#   mod_mk <- MannKendall(analysis_df$max_q7)
#   mod_sens <- sens.slope(analysis_df$max_q7[!is.na(analysis_df$max_q7)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = water_year, y = max_q7, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Streamflow [mm]", title = paste(scenario, "Max. 7 Day Flow (Q7 Max)"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$max_q7), label = label, color = "black", hjust = 1, vjust = 1) 
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Q7Max_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# 
# 
# #######################################################################
# ### Customizable plots: metric, season, etc ###
# # EDIT TO PLOT SPECIFIC FUTURE SCENARIOS
# 
# # define flow level in cfs, convert to mm
# flow_level <- 1800 * 28316847*86400/(2590000000000*sqmi) 
# # identify months of interest (numerical values)
# mos <- c("02", "03") 
# # is comparison above or below threshold?
# comparison = 'above'
# 
# 
# # calculate number of days above/below threshold
# if (tolower(comparison)=='below') {
#   threshold <- as.data.frame(daily_df %>% filter(mo %in% mos) %>% mutate(flow = ifelse(total <= flow_level, 1, 0)) %>%
#                                     group_by(water_year, projection) %>% dplyr::summarize(projection=first(projection), days = sum(flow)))# %>%filter(!is.na(days)))
# } else if (tolower(comparison)=='above'){
#   threshold <- as.data.frame(daily_df %>% filter(mo %in% mos) %>% mutate(flow = ifelse(total >= flow_level, 1, 0)) %>%
#                                     group_by(water_year, projection) %>% dplyr::summarize(projection=first(projection), days = sum(flow)))# %>%filter(!is.na(days)))
# }
# 
# # plot
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   analysis_df <- threshold %>% filter(projection=='Historical' | projection==proj) 
#   
#   mod_mk <- MannKendall(analysis_df$days)
#   mod_sens <- sens.slope(analysis_df$days[!is.na(analysis_df$days)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   
#   plot_mod <- ggplot(analysis_df, aes(x = water_year, y = days, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Days", title = paste('Days', comparison, round(flow_level), 'mm, ', month.abb[as.numeric(mos)][1], '-', month.abb[as.numeric(mos)][length(mos)]), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_Days_", comparison, "_", round(flow_level), "mm_", month.abb[as.numeric(mos)][1], "_", month.abb[as.numeric(mos)][length(mos)],".jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# 
# 
# #######################################################################
# ### Non-streamflow metrics ###
# 
# # Calculate runoff efficiency
# runoff_efficiency <- daily_df %>% group_by(water_year, projection) %>% dplyr::summarize(projection=first(projection), ann_runoff=sum(adj_runoff))
# hist_p <- DailyClimData %>% mutate(water_year = sapply(date, get_water_year)) %>% group_by(water_year) %>% dplyr::summarize(ann_p=sum(pr)) %>% mutate(projection='Historical') %>% filter(water_year != 2023)
# fut_p <- future_climate %>% mutate(water_year = sapply(date, get_water_year)) %>% group_by(water_year, projection) %>% dplyr::summarize(projection=first(projection), ann_p=sum(pr))
# runoff_efficiency <- runoff_efficiency %>% full_join(bind_rows(hist_p, fut_p), by = c('water_year', 'projection'))
# runoff_efficiency$efficiency <- runoff_efficiency$ann_runoff / runoff_efficiency$ann_p
# 
# # Plot
# plot_list <- list()
# for (i in 1:length(model_names)){
#   proj = model_names[i]
#   scenario <- scenario_names[i]
#   
#   analysis_df <- runoff_efficiency %>% filter(projection=='Historical' | projection==proj) %>% filter(water_year != 2100)
#   mod_mk <- MannKendall(analysis_df$efficiency)
#   mod_sens <- sens.slope(analysis_df$efficiency[!is.na(analysis_df$efficiency)])
#   if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
#   }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
#   plot_mod <- ggplot(analysis_df, aes(x = water_year, y = efficiency, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
#     geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
#     labs(x = "Water Year", y = "Runoff Efficiency (%)", title = paste(scenario, "Runoff Efficiency"), color='') +
#     nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
#     annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$efficiency), label = label, color = "black", hjust = 1, vjust = 1) 
#   print(plot_mod)
#   plot_list[[i]] <- plot_mod
# }
# jpeg(file=paste0(outLocationPathFuture, "/", "Modeled_RunoffEffiency_Trends.jpg"), width=300*num_models, height=200*num_models)
# grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
