# ---------------------------------------------------------------------
# This script includes code to generate and analyze future water balance projections,
# which are either calculated using pre-calibrated coefficients or pulled from a pre-generated
# gridded water balance model produced/maintained by Mike Tercek. This code also provides 
# visualizations and preliminary analysis of these future water balance projections. 
# 
# ---------------------------------------------------------------------

#######################################################################
### GENERATE FUTURE WATER BALANCE MODELS ###

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
if(exists("future_wb_conus") & exists("future_wb_calc")){
  model_run = 'HadGEM2-CC365.rcp45'; yr = 2060
  if(make_plots){
    plot_aet <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=AET), col='black')+
      geom_line(data=future_wb_calc %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=AET), col='red')+
      labs(x='Date',y='AET [mm]', title='Actual Evapotranspiration') +
      theme(legend.position = "none") + nps_theme()
    plot_d <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(date)==yr), aes(x=date, y=deficit), col='black')+
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
### AGGREGATE FUTURE AND HISTORICAL PROJECTIONS ###

# Define which future water balance projection to use
if(calcFutureWB){
  future_wb <- future_wb_calc %>% select(date, projection, D, AET, SOIL, RUNOFF, SNOW, PET, adj_runoff)
  colnames(future_wb) <- c('date','projection','deficit','AET','soil_water','runoff','accumswe','PET','adj_runoff')
  if(!dir.exists(file.path(outLocationPath, 'WB_Calc'))) {dir.create(file.path(outLocationPath, 'WB_Calc'))}; outLocationPathFuture = file.path(outLocationPath, 'WB_Calc')
} else{
  future_wb <- future_wb_conus %>% select(-rain)
  if(!dir.exists(file.path(outLocationPath, 'WB_CONUS'))) {dir.create(file.path(outLocationPath, 'WB_CONUS'))}; outLocationPathFuture = file.path(outLocationPath, 'WB_CONUS')
}

hist_wb <- DailyWB %>% mutate(projection='Historical') %>% relocate(projection, .after = 1) %>% select(date, projection, D, AET, SOIL, RUNOFF, SNOW, PET, adj_runoff)
colnames(hist_wb) <- c('date','projection','deficit','AET','soil_water','runoff','accumswe','PET','adj_runoff')
full_wb <- bind_rows(hist_wb, future_wb); full_wb <- full_wb %>% select(where(~ !any(is.na(.)))) 

# Add temporal variables
full_wb$date <- as.Date(full_wb$date); full_wb$yr<-as.numeric(format(full_wb$date,"%Y")); full_wb$yr_mo<-format(full_wb$date,"%Y-%m")
full_wb$yday <- factor(yday(full_wb$date), levels = unique(full_wb$yday))
full_wb$water_year <- sapply(full_wb$date, get_water_year); 
full_wb <- full_wb %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
full_wb$Period<-ifelse(full_wb$yr<=2005,"Historical",ifelse (full_wb$yr>=2006 & full_wb$yr<=2050,"Early",
                                                               ifelse (full_wb$yr>=2051 & full_wb$yr<=2070,"Middle", ifelse (full_wb$yr>=2071, "Late","NA"))))

# Aggregate data to annual
annual_df <- as.data.frame(full_wb %>% group_by(projection, yr) %>%
                             dplyr::summarize(projection=first(projection), runoff=sum(runoff, na.rm=TRUE), adj_runoff = sum(adj_runoff, na.rm = TRUE),
                                              PET = sum(PET, na.rm = TRUE), AET = sum(AET, na.rm = TRUE), deficit = sum(deficit, na.rm = TRUE), 
                                              accumswe = max(accumswe, na.rm = TRUE), soil_water=mean(soil_water, na.rm=TRUE)))
annual_df$date<-as.Date(paste(annual_df$yr,"-01", "-01",sep=""))

# Aggregate data to monthly 
monthly_df <- as.data.frame(full_wb %>% group_by(projection, yr_mo) %>%
                              dplyr::summarize(projection=first(projection), runoff=sum(runoff, na.rm=TRUE), adj_runoff = sum(adj_runoff, na.rm = TRUE),
                                               PET = sum(PET, na.rm = TRUE), AET = sum(AET, na.rm = TRUE), deficit = sum(deficit, na.rm = TRUE), 
                                               accumswe = max(accumswe, na.rm = TRUE), soil_water=mean(soil_water, na.rm=TRUE)))
monthly_df$date<-as.Date(paste(monthly_df$yr_mo,"-01",sep=""))



#######################################################################
### SELECT DIVERGENT CLIMATE FUTURES ###
#scenario_names <- c("Warm Wet", "Hot Wet", "Warm Dry", "Hot Dry")
color_names <- c("#4E5BA6","#E10720","#B4B8F0","#E79C9C")
future_means <- select_climate_futures(color_names)


### Identify models in format for plotting ###
model_names <- (future_means %>% filter(!is.na(pca)))$projection; scenario_names <- (future_means %>% filter(!is.na(pca)))$pca
num_models <- length(model_names)

# Reorder for prettier plotting
order_index <- match(c("Warm Wet", "Hot Wet", "Warm Dry", "Hot Dry"), scenario_names)
model_names <- model_names[order_index]
scenario_names <- scenario_names[order_index]

# identify warm wet/hot dry [1:2] or warm dry/hot wet [3:4] or all (comment out both lines)
#model_names <- model_names[1:2]; scenario_names <- scenario_names[1:2]; color_names <- color_names[1:2]
#model_names <- model_names[3:4]; scenario_names <- scenario_names[3:4]; color_names <- color_names[3:4]



#######################################################################
### GET MACA DATA FOR CLIMATE FUTURES ###
# edit so only 4 climate futures are selected
if(point_location){
  if(!file.exists(here('Data',SiteID_FileName,paste('Historical_WB_Projections',SiteID_FileName, "1960_2005_point.csv", sep='_')))){
    maca_hist <- get_maca_hist_point(lat, lon, 1960, 2005, SiteID_FileName, model_names)
    maca_hist$gcm <- maca_hist$projection; maca_hist$projection <- paste0(maca_hist$projection, ".Hist")
    
    
    # run WB model
    hist <- NULL
    for(mod in model_names){
      hist_clim <- maca_hist %>% filter(grepl(strsplit(mod, "\\.")[[1]][1], gcm))
      
      hist_wb <- WB(hist_clim, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect, slope,
                    shade.coeff, jtemp,SWC.Max, 1, 1, 0, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon, "")
      hist_wb <- hist_wb %>% select(projection, date, RUNOFF, adj_runoff, PET, AET, D, SNOW, SOIL)
      hist_wb <- hist_wb %>% rename(runoff=RUNOFF, deficit=D, accumswe=SNOW, soil_water=SOIL)
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
      hist_wb <- hist_wb %>% select(projection, date, RUNOFF, adj_runoff, PET, AET, D, SNOW, SOIL)
      hist_wb <- hist_wb %>% rename(runoff=RUNOFF, deficit=D, accumswe=SNOW, soil_water=SOIL)
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
future_wb$rcp <- sapply(strsplit(as.character(future_wb$projection), "\\."), function(x) gsub("\\D", "", x[2]))
model_df <- rbind(hist, future_wb) 

# Add useful data
model_df$gcm <- sapply(strsplit(as.character(model_df$projection), "\\."), `[`, 1)
model_df$yr<-as.numeric(format(model_df$date,"%Y")); model_df$mo<-format(model_df$date,"%m"); model_df$yr_mo<-format(model_df$date,"%Y-%m"); model_df$day <- factor(yday(model_df$date), levels = unique(yday(model_df$date)))
model_df$water_year <- sapply(model_df$date, get_water_year); 
model_df <- model_df %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
model_df<-model_df[,c("date", "projection", "gcm", "rcp", "yr", "mo", "yr_mo", "water_year", "runoff", "adj_runoff", "PET", "AET", "deficit", "accumswe", "soil_water")]
model_df$Period<-ifelse(model_df$yr<=2005,"Historical",ifelse (model_df$yr>=2006 & model_df$yr<=2050,"Early",
                                                             ifelse (model_df$yr>=2051 & model_df$yr<=2070,"Middle", ifelse (model_df$yr>=2071, "Late","NA"))))

# Aggregate data to annual
model_annual_df <- as.data.frame(model_df %>% group_by(projection, yr) %>%
                                   dplyr::summarize(gcm=first(gcm), rcp=first(rcp), water_year=first(water_year),
                                                    runoff=sum(runoff, na.rm=TRUE), adj_runoff = sum(adj_runoff, na.rm = TRUE), PET=sum(PET, na.rm = TRUE),
                                                    AET=sum(AET, na.rm=TRUE), deficit=sum(deficit, na.rm = TRUE), accumswe=max(accumswe, na.rm = TRUE), 
                                                    soil_water=mean(soil_water, na.rm=TRUE), Period=first(Period)))
model_annual_df$date<-as.Date(paste(model_annual_df$yr,"-01", "-01",sep=""))

# Aggregate data to monthly 
model_monthly_df <- as.data.frame(model_df %>% group_by(projection, yr_mo) %>%
                                    dplyr::summarize(gcm=first(gcm), rcp=first(rcp), water_year=first(water_year),
                                                     runoff=sum(runoff, na.rm=TRUE), adj_runoff = sum(adj_runoff, na.rm = TRUE), PET=sum(PET, na.rm = TRUE),
                                                     AET=sum(AET, na.rm=TRUE), deficit=sum(deficit, na.rm = TRUE), accumswe=max(accumswe, na.rm = TRUE), 
                                                     soil_water=mean(soil_water, na.rm=TRUE), Period=first(Period)))
model_monthly_df$date<-as.Date(paste(model_monthly_df$yr_mo,"-01",sep=""))




#######################################################################
#######################################################################
### TIME SERIES PLOTS ###
# MACA historical data

# Annual AET
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(gcm==strsplit(proj, "\\.")[[1]][1] & (rcp==gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]) | grepl('Hist', rcp)))
    mod_mk <- MannKendall(analysis_df$AET)
    mod_sens <- sens.slope(analysis_df$AET[!is.na(analysis_df$AET)])
    if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
    }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
    plot_mod <- ggplot(analysis_df, aes(x = yr, y = AET)) + geom_line(aes(color=factor(rcp)), na.rm=TRUE, linewidth=1, alpha=0.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
      labs(x = "Water Year", y = "Annual AET (mm)", title = paste(scenario, "Annual Modeled AET"), color='') +
      nps_theme() + theme(legend.position = 'bottom') + 
      scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c('45'=proj, '85'=proj)) + 
      annotate("text", x = max(analysis_df$yr), y = max(analysis_df$AET), label = label, color = "black", hjust = 1, vjust = 1) 
    print(plot_mod)
    plot_list[[i]] <- plot_mod
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_AET_Trends.jpg"), width=300*num_models, height=200*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
  dev.off()
}

# Annual PET
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(gcm==strsplit(proj, "\\.")[[1]][1] & (rcp==gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]) | grepl('Hist', rcp)))
    mod_mk <- MannKendall(analysis_df$PET)
    mod_sens <- sens.slope(analysis_df$PET[!is.na(analysis_df$PET)])
    if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
    }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
    plot_mod <- ggplot(analysis_df, aes(x = yr, y = PET, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
      labs(x = "Water Year", y = "Annual PET (mm)", title = paste(scenario, "Annual Modeled PET"), color='') +
      nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c('45'=proj, '85'=proj)) + 
      annotate("text", x = max(analysis_df$yr), y = max(analysis_df$PET), label = label, color = "black", hjust = 1, vjust = 1) 
    print(plot_mod)
    plot_list[[i]] <- plot_mod
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_PET_Trends.jpg"), width=300*num_models, height=200*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
  dev.off()
}

# Annual deficit
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(gcm==strsplit(proj, "\\.")[[1]][1] & (rcp==gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]) | grepl('Hist', rcp)))
    mod_mk <- MannKendall(analysis_df$deficit)
    mod_sens <- sens.slope(analysis_df$deficit[!is.na(analysis_df$deficit)])
    if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
    }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
    plot_mod <- ggplot(analysis_df, aes(x = yr, y = deficit, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
      labs(x = "Water Year", y = "Annual Deficit (mm)", title = paste(scenario, "Annual Modeled Deficit"), color='') +
      nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c('45'=proj, '85'=proj)) + 
      annotate("text", x = max(analysis_df$yr), y = max(analysis_df$deficit), label = label, color = "black", hjust = 1, vjust = 1) 
    print(plot_mod)
    plot_list[[i]] <- plot_mod
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Deficit_Trends.jpg"), width=300*num_models, height=200*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
  dev.off()
}

# Annual SWE
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(gcm==strsplit(proj, "\\.")[[1]][1] & (rcp==gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]) | grepl('Hist', rcp)))
    mod_mk <- MannKendall(analysis_df$accumswe)
    mod_sens <- sens.slope(analysis_df$accumswe[!is.na(analysis_df$accumswe)])
    if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
    }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
    plot_mod <- ggplot(analysis_df, aes(x = yr, y = accumswe, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
      labs(x = "Water Year", y = "Max Annual SWE (mm)", title = paste(scenario, "Annual Modeled Snow Water Equivalent"), color='') +
      nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c('45'=proj, '85'=proj)) + 
      annotate("text", x = max(analysis_df$yr), y = max(analysis_df$accumswe), label = label, color = "black", hjust = 1, vjust = 1) 
    print(plot_mod)
    plot_list[[i]] <- plot_mod
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_SWE_Trends.jpg"), width=300*num_models, height=200*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
  dev.off()
}

# Annual runoff
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(gcm==strsplit(proj, "\\.")[[1]][1] & (rcp==gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]) | grepl('Hist', rcp)))
    mod_mk <- MannKendall(analysis_df$runoff)
    mod_sens <- sens.slope(analysis_df$runoff[!is.na(analysis_df$runoff)])
    if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
    }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
    plot_mod <- ggplot(analysis_df, aes(x = yr, y = runoff, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
      labs(x = "Water Year", y = "Annual Runoff (mm)", title = paste(scenario, "Annual Modeled Runoff"), color='') +
      nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c('45'=proj, '85'=proj)) + 
      annotate("text", x = max(analysis_df$yr), y = max(analysis_df$runoff), label = label, color = "black", hjust = 1, vjust = 1) 
    print(plot_mod)
    plot_list[[i]] <- plot_mod
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Runoff_Trends.jpg"), width=300*num_models, height=200*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
  dev.off()
}

# Soil water content
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(gcm==strsplit(proj, "\\.")[[1]][1] & (rcp==gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]) | grepl('Hist', rcp)))
    mod_mk <- MannKendall(analysis_df$soil_water)
    mod_sens <- sens.slope(analysis_df$soil_water[!is.na(analysis_df$soil_water)])
    if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
    }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
    plot_mod <- ggplot(analysis_df, aes(x = yr, y = soil_water, color=factor(rcp))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
      geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
      labs(x = "Water Year", y = "Mean SWC (mm)", title = paste(scenario, "Annual Modeled Soil Water Content"), color='') +
      nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', '45'=color_names[i], '85'=color_names[i], "Trend"="black"), labels = c('45'=proj, '85'=proj)) + 
      annotate("text", x = max(analysis_df$yr), y = max(analysis_df$soil_water), label = label, color = "black", hjust = 1, vjust = 1) 
    print(plot_mod)
    plot_list[[i]] <- plot_mod
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_SWC_Trends.jpg"), width=300*num_models, height=200*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
  dev.off() 
}


#######################################################################
#######################################################################
### BOX PLOTS ###

### Change in annual AET ###
delta_plot <- (model_annual_df %>% filter(rcp!="Historical")) %>% 
  left_join((model_annual_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(mean_AET=mean(AET, na.rm=TRUE)), by = "gcm") %>%
  mutate(delta = AET-mean_AET)

if(make_plots){
  # All models by RCP
  plot <- ggplot() + geom_boxplot(data=delta_plot, aes(x=projection, y=delta, fill=rcp)) + 
    labs(title='Change in annual AET', y='Change in AET [mm]', x='Projection', fill='RCP') + 
    scale_fill_manual(values = c('45' = 'orange', '85' = 'red')) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_AET_Change_RCP.jpg"), width=600, height=500)
  print(plot); dev.off()
  
  # All models by climate future
  plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta, fill=projection)) + 
    geom_boxplot(data=(delta_plot %>% filter(!projection %in% model_names)), aes(x=projection, y=delta, fill='Other')) + 
    labs(title='Change in annual AET', y='Change in AET [mm]', x='Projection', fill='Future') + 
    scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_AET_Change.jpg"), width=600, height=500)
  print(plot); dev.off()
  
  # Just climate futures
  plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta, fill=projection)) + 
    labs(title='Change in annual AET', y='Change in AET [mm]', x='Projection', fill='Future') + 
    scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_AET_Change_Futures.jpg"), width=600, height=500)
  print(plot); dev.off()
}


### Change in annual PET ###
delta_plot <- (model_annual_df %>% filter(rcp!="Historical")) %>% 
  left_join((model_annual_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(mean_PET=mean(PET, na.rm=TRUE)), by = "gcm") %>%
  mutate(delta = PET-mean_PET)

if(make_plots){
  # All models by RCP
  plot <- ggplot() + geom_boxplot(data=delta_plot, aes(x=projection, y=delta, fill=rcp)) + 
    labs(title='Change in annual PET', y='Change in PET [mm]', x='Projection', fill='RCP') + 
    scale_fill_manual(values = c('45' = 'orange', '85' = 'red')) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_PET_Change_RCP.jpg"), width=600, height=500)
  print(plot); dev.off()
  
  # All models by climate future
  plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta, fill=projection)) + 
    geom_boxplot(data=(delta_plot %>% filter(!projection %in% model_names)), aes(x=projection, y=delta, fill='Other')) + 
    labs(title='Change in annual PET', y='Change in PET [mm]', x='Projection', fill='Future') + 
    scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_PET_Change.jpg"), width=600, height=500)
  print(plot); dev.off()
  
  # Just climate futures
  plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta, fill=projection)) + 
    labs(title='Change in annual PET', y='Change in PET [mm]', x='Projection', fill='Future') + 
    scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_PET_Change_Futures.jpg"), width=600, height=500)
  print(plot); dev.off()
}


### Change in annual deficit ###
delta_plot <- (model_annual_df %>% filter(rcp!="Historical")) %>% 
  left_join((model_annual_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(mean_deficit=mean(deficit, na.rm=TRUE)), by = "gcm") %>%
  mutate(delta = deficit-mean_deficit)

if(make_plots){
  # All models by RCP
  plot <- ggplot() + geom_boxplot(data=delta_plot, aes(x=projection, y=delta, fill=rcp)) + 
    labs(title='Change in annual deficit', y='Change in deficit [mm]', x='Projection', fill='RCP') + 
    scale_fill_manual(values = c('45' = 'orange', '85' = 'red')) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Deficit_Change_RCP.jpg"), width=600, height=500)
  print(plot); dev.off()
  
  # All models by climate future
  plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta, fill=projection)) + 
    geom_boxplot(data=(delta_plot %>% filter(!projection %in% model_names)), aes(x=projection, y=delta, fill='Other')) + 
    labs(title='Change in annual deficit', y='Change in deficit [mm]', x='Projection', fill='Future') + 
    scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Deficit_Change.jpg"), width=600, height=500)
  print(plot); dev.off()
  
  # Just climate futures
  plot <- ggplot() + geom_boxplot(data=(delta_plot %>% filter(projection %in% model_names)), aes(x=projection, y=delta, fill=projection)) + 
    labs(title='Change in annual deficit', y='Change in deficit [mm]', x='Projection', fill='Future') + 
    scale_fill_manual(values = c("Other" = "gray", setNames(color_names, model_names)), labels = c(setNames(scenario_names, model_names))) + 
    nps_theme() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Deficit_Change_Futures.jpg"), width=600, height=500)
  print(plot); dev.off()
}



#######################################################################
#######################################################################
### CLIMATE FUTURE SCATTERPLOTS ### 

# AET vs D
delta_plot <- (model_annual_df %>% filter(rcp!="Historical")) %>% 
  left_join((model_annual_df %>% filter(rcp=="Historical")) %>% group_by(gcm) %>% summarize(mean_AET=mean(AET, na.rm=TRUE), mean_deficit=mean(deficit, na.rm=TRUE)), by = "gcm") %>%
  mutate(delta_AET=AET-mean_AET, delta_deficit=deficit-mean_deficit)
delta_plot <- delta_plot %>% filter(gcm!='Historical') %>% group_by(projection, Period) %>% dplyr::summarise(gcm=first(gcm), rcp=first(rcp), delta_AET=mean(delta_AET), delta_deficit=mean(delta_deficit))

annual_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                              xintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_deficit, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_deficit, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_deficit, 0.5)))
sd_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                          yintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_AET, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_AET, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_AET, 0.5)))
annual_zero <- data.frame(Period = c("Early", "Middle", "Late"), xintercept=c(0,0,0)); sd_zero <- data.frame(Period = c("Early", "Middle", "Late"), yintercept = c(0,0,0))


if(make_plots){
  # plot with RCPs
  jpeg(file=paste0(outLocationPathFuture, "/WB_Climate_Future_Scatter_RCP.jpg"), width=1000, height=400)
  plot <- ggplot(data=delta_plot, aes(x=delta_deficit,y=delta_AET,color=rcp)) + geom_point() +
    geom_text_repel(aes(label = gcm), color = 'black', max.overlaps=Inf) +
    geom_hline(data = sd_quantile, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_quantile, aes(xintercept = xintercept), color = "black") +
    #geom_hline(data = sd_zero, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_zero, aes(xintercept = xintercept), color = "black") +
    facet_wrap(~factor(Period, levels = c('Early', 'Middle', 'Late'))) + labs(title=paste('Changes in climatic water balance at',SiteID), x='Change in deficit [mm]', y='Change in AET [mm]', color='RCP') + 
    scale_color_manual(values = c("45" = "orange", "85" = "red")) + nps_theme()
  print(plot)
  dev.off()
  
  # plot with selected climate futures
  jpeg(file=paste0(outLocationPathFuture, "/WB_Climate_Future_Scatter.jpg"), width=1000, height=400)
  plot <- ggplot(data=delta_plot) + 
    geom_point(data=delta_plot %>% filter(!(projection %in% model_names)), aes(x=delta_deficit,y=delta_AET,color='Other')) +
    geom_point(data=delta_plot %>% filter(projection %in% model_names), aes(x=delta_deficit,y=delta_AET,color=projection)) + 
    geom_text_repel(data=delta_plot %>% filter(!(projection %in% model_names)), aes(x=delta_deficit,y=delta_AET, label=projection), color = 'black', size=2, max.overlaps=Inf) +
    geom_text_repel(data=delta_plot %>% filter(projection %in% model_names), aes(x=delta_deficit,y=delta_AET, label=projection, color=projection), size=4,  fontface = "bold", max.overlaps=Inf) + 
    geom_hline(data = sd_quantile, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_quantile, aes(xintercept = xintercept), color = "black") +
    #geom_hline(data = sd_zero, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_zero, aes(xintercept = xintercept), color = "black") +
    facet_wrap(~factor(Period, levels = c('Early', 'Middle', 'Late'))) + labs(title=paste('Changes in climatic water balance at',SiteID), x='Change in deficit [mm]', y='Change in AET [mm]', color='RCP') + 
    scale_color_manual(values = c("Other" = "black", setNames(color_names, model_names)), labels=c("Other" = "black", setNames(scenario_names, model_names))) + nps_theme()
  print(plot)
  dev.off()
}



#######################################################################
#######################################################################
### CLIMATE FUTURE BLOBS ###
# AET vs D
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_annual_df %>% filter(gcm==strsplit(proj, "\\.")[[1]][1] & (rcp==gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]) | grepl('Hist', rcp)))
    
    plot_list[[i]] <- ggplot(analysis_df) + geom_point(aes(x = deficit, y = AET, color = Period), size = 4, alpha = 0.6) +
      geom_encircle(aes(x = deficit, y = AET, color = Period), show.legend = FALSE, size = 2) +
      labs(title = paste(scenario, "Deficit vs. AET"), x = "Annual Deficit [mm]", y = "Annual AET [mm]") +
      #scale_color_manual(values = c("Historical"="#A3B86C", "Early"="#EBC944", "Middle"="#1496BB", "Late"='Black')) +
      theme(axis.title.x = element_text(size = 15)) + nps_theme()
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Annual_AET_D_Trends.jpg"), width=300*num_models, height=200*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
  dev.off() 
}


#######################################################################
#######################################################################
### SEASONALITY ###

# AET distribution
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
      dplyr::summarize(gcm=first(gcm), rcp=first(rcp), value = mean(AET, na.rm = TRUE), .groups = 'drop') %>%
      mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
    
    plot_list[[i]] <- ggplot() + 
      geom_smooth(data = (analysis_df %>% filter(!grepl(strsplit(proj, "\\.")[[1]][1], gcm) & !(grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp))), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp))), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      labs(title=paste("Smoothed Daily AET"), y=paste('AET [mm]'), x=NULL, color='Model') + nps_theme() + 
      scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
      guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_AET_Trends.jpg"), width=300*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}

# PET distribution
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
      dplyr::summarize(gcm=first(gcm), rcp=first(rcp), value = mean(PET, na.rm = TRUE), .groups = 'drop') %>%
      mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
    
    plot_list[[i]] <- ggplot() + 
      geom_smooth(data = (analysis_df %>% filter(!grepl(strsplit(proj, "\\.")[[1]][1], gcm) & !(grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp))), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp))), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      labs(title=paste("Smoothed Daily PET"), y=paste('PET [mm]'), x=NULL, color='Model') + nps_theme() + 
      scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
      guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_PET_Trends.jpg"), width=300*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}

# Deficit distribution
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
      dplyr::summarize(gcm=first(gcm), rcp=first(rcp), value = mean(deficit, na.rm = TRUE), .groups = 'drop') %>%
      mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
    
    plot_list[[i]] <- ggplot() + 
      geom_smooth(data = (analysis_df %>% filter(!grepl(strsplit(proj, "\\.")[[1]][1], gcm) & !(grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp))), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp))), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      labs(title=paste("Smoothed Daily Deficit"), y=paste('Deficit [mm]'), x=NULL, color='Model') + nps_theme() + 
      scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
      guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_Deficit_Trends.jpg"), width=300*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}

# Soil water
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
      dplyr::summarize(gcm=first(gcm), rcp=first(rcp), value = mean(soil_water, na.rm = TRUE), .groups = 'drop') %>%
      mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
    
    plot_list[[i]] <- ggplot() + 
      geom_smooth(data = (analysis_df %>% filter(!grepl(strsplit(proj, "\\.")[[1]][1], gcm) & !(grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp))), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp))), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      labs(title=paste("Smoothed Daily Soil Water Content"), y=paste('SWC [mm]'), x=NULL, color='Model') + nps_theme() + 
      scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
      guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_SWC_Trends.jpg"), width=300*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}

# SWE
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
      dplyr::summarize(gcm=first(gcm), rcp=first(rcp), value = mean(accumswe, na.rm = TRUE), .groups = 'drop') %>%
      mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
    
    plot_list[[i]] <- ggplot() + 
      geom_smooth(data = (analysis_df %>% filter(!grepl(strsplit(proj, "\\.")[[1]][1], gcm) & !(grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp))), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp))), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      labs(title=paste("Smoothed Daily Snow Water Equivalent"), y=paste('SWE [mm]'), x=NULL, color='Model') + nps_theme() + 
      scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
      guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_SWE_Trends.jpg"), width=300*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}

# Runoff
if(make_plots){
  plot_list <- list()
  for (i in 1:length(model_names)){
    proj = model_names[i]
    scenario <- scenario_names[i]
    
    analysis_df <- model_df %>% mutate(date_group = yday(date)) %>% group_by(projection, date_group) %>%
      dplyr::summarize(gcm=first(gcm), rcp=first(rcp), value = mean(runoff, na.rm = TRUE), .groups = 'drop') %>%
      mutate(date_group = factor(date_group, levels = c(274:366, 1:273))) %>%  arrange(date_group)
    
    plot_list[[i]] <- ggplot() + 
      geom_smooth(data = (analysis_df %>% filter(!grepl(strsplit(proj, "\\.")[[1]][1], gcm) & !(grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp) | grepl('Hist', rcp)))), aes(x = date_group, y = value, group = projection, color = "Other"), alpha = 0.1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl(gsub("\\D", "", strsplit(proj, "\\.")[[1]][2]), rcp))), aes(x = date_group, y = value, group = projection, color = projection), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      geom_smooth(data = (analysis_df %>% filter(grepl(strsplit(proj, "\\.")[[1]][1], gcm) & grepl('Hist', rcp))), aes(x = date_group, y = value, group = projection, color = "Historical"), alpha = 1,  
                  method = "gam", formula = y ~ s(x, bs = "cs"), se = FALSE) + 
      labs(title=paste("Smoothed Daily Runoff"), y=paste('Runoff [mm]'), x=NULL, color='Model') + nps_theme() + 
      scale_color_manual(values = c("Other" = "gray", 'Historical'='black',setNames(color_names, model_names)), labels = c("Other" = 'Other Models', setNames(scenario_names, model_names))) + 
      guides(alpha = "none") + scale_x_discrete(breaks = c(274,305,335,1,32,60,91,121,152,182,213,244),labels = c("Oct", "Nov","Dec","Jan","Feb","Mar","Apr","May","Jun",'Jul','Aug','Sep'))
  }
  jpeg(file=paste0(outLocationPathFuture, "/", "Smoothed_Daily_Runoff_Trends.jpg"), width=300*num_models, height=150*num_models)
  grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2); dev.off()
}
