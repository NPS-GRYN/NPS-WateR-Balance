# ---------------------------------------------------------------------
# This script includes code to generate and analyze future water balance projections,
# which are either calculated using pre-calibrated coefficients or pulled from a pre-generated
# gridded water balance model produced/maintained by Mike Tercek. This code also provides 
# visualizations and preliminary analysis of these future water balance projections. 
# 
# EDITS IN PROGRESS
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
  # Access provided file
  if(file.exists(filename_future_wb)) {
    future_wb_conus <- get_conus_wb_direct(SiteID_FileName, dataPath, filename_future_wb)
  } else{
    # Pull directly from website
    future_wb_conus <- get_conus_wb(SiteID_FileName, lat, lon, endY, 2099)
  }
  
  # If neither version can be accessed, calculate water balance
  if(anyNA(future_wb_conus)){
    calcFutureWB <- TRUE
  }
}



### Re-run water balance model to generate future projections ###
# has the name projection been changed?
if(calcFutureWB){
  if(!file.exists(here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, endY, "2100.csv", sep='_')))){
    # Get future climate data
    if(point_location) {future_climate <- get_maca_point(lat, lon, SiteID_FileName, gcm_list)
    } else future_climate <- get_maca_area(aoi, SiteID_FileName, gcm_list)
    
    # Run water balance code for each future projection
    wb_list <- vector("list", length(unique(future_climate$projection)))
    i <- 1
    for(proj in unique(future_climate$projection)){
      print(proj)
      ClimData <- future_climate %>% filter(projection==proj) %>% select(-projection)
      DailyWB_future <- WB(ClimData, gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope,
                           shade.coeff, jtemp, SWC.Max, aet_slope, aet_bias, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
      wb_list[[i]] <- cbind(proj, DailyWB_future)
      i <- i + 1
    }
    future_wb_calc <- do.call(rbind, wb_list); future_wb_calc <- future_wb_calc %>% rename(projection=proj)
    
    
    # Save calculated WB
    write.csv(future_wb_calc, here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, endY, "2100.csv", sep='_')), row.names=FALSE)
    } else{
   # Read in calculated WB
   future_wb_calc <- read.csv(here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, endY, "2100.csv", sep='_')))
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
    while (dev.cur() > 1) dev.off()
  }
}


#######################################################################
### AGGREGATE FUTURE AND HISTORICAL PROJECTIONS ###

# Define which future water balance projection to use
if(calcFutureWB){
  future_wb <- future_wb_calc %>% select(date, projection, D, AET, SOIL, RUNOFF, RAIN, SNOW, PET, adj_runoff)
  colnames(future_wb) <- c('date','projection','deficit','AET','soil_water','runoff','rain','accumswe','PET','adj_runoff')
  if(!dir.exists(file.path(outLocationPath, 'WB_Calc'))) {dir.create(file.path(outLocationPath, 'WB_Calc'))}; outLocationPathFuture = file.path(outLocationPath, 'WB_Calc')
} else{
  future_wb <- future_wb_conus 
  if(!dir.exists(file.path(outLocationPath, 'WB_CONUS'))) {dir.create(file.path(outLocationPath, 'WB_CONUS'))}; outLocationPathFuture = file.path(outLocationPath, 'WB_CONUS')
}

hist_wb <- DailyWB %>% mutate(projection='Historical') %>% relocate(projection, .after = 1) %>% select(date, projection, D, AET, SOIL, RUNOFF, RAIN, SNOW, PET, adj_runoff)
colnames(hist_wb) <- c('date','projection','deficit','AET','soil_water','runoff','rain','accumswe','PET','adj_runoff')
full_wb <- bind_rows(hist_wb, future_wb); full_wb <- full_wb %>% select(where(~ !any(is.na(.)))) 

# Add temporal variables
full_wb$date <- as.Date(full_wb$date); full_wb$yr<-as.numeric(format(full_wb$date,"%Y"))
full_wb$yday <- factor(yday(full_wb$date), levels = unique(full_wb$yday))
full_wb$water_year <- sapply(full_wb$date, get_water_year); 
full_wb <- full_wb %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
full_wb$Period<-ifelse(full_wb$yr<=2022,"Historical",ifelse (full_wb$yr>=2023 & full_wb$yr<=2050,"Early",
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
future_means <- select_climate_futures()


### Identify models in format for plotting ###
model_names <- (future_means %>% filter(!is.na(pca)))$projection; scenario_names <- (future_means %>% filter(!is.na(pca)))$pca
num_models <- length(model_names)

# identify warm wet/hot dry [1:2] or warm dry/hot wet [3:4] or all (comment out both lines)
#model_names <- model_names[1:2]; scenario_names <- scenario_names[1:2]; color_names <- color_names[1:2]
#model_names <- model_names[3:4]; scenario_names <- scenario_names[3:4]; color_names <- color_names[3:4]



#######################################################################
### PLOT FUTURE WB FOR CLIMATE FUTURES ###

# Annual AET
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- annual_df %>% filter(projection=='Historical' | projection==proj)
  mod_mk <- MannKendall(analysis_df$AET)
  mod_sens <- sens.slope(analysis_df$AET[!is.na(analysis_df$AET)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_mod <- ggplot(analysis_df, aes(x = yr, y = AET, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Annual Streamflow (mm)", title = paste(scenario, "Annual Modeled Streamflow"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$AET), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_mod)
  plot_list[[i]] <- plot_mod
}
jpeg(file=paste0(outLocationPathFuture, "/", "Annual_AET_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()

# Annual PET
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- annual_df %>% filter(projection=='Historical' | projection==proj)
  mod_mk <- MannKendall(analysis_df$PET)
  mod_sens <- sens.slope(analysis_df$PET[!is.na(analysis_df$PET)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_mod <- ggplot(analysis_df, aes(x = yr, y = PET, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Annual Streamflow (mm)", title = paste(scenario, "Annual Modeled Streamflow"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$PET), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_mod)
  plot_list[[i]] <- plot_mod
}
jpeg(file=paste0(outLocationPathFuture, "/", "Annual_PET_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()

# Annual deficit
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- annual_df %>% filter(projection=='Historical' | projection==proj)
  mod_mk <- MannKendall(analysis_df$deficit)
  mod_sens <- sens.slope(analysis_df$deficit[!is.na(analysis_df$deficit)])
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)
  }else{label <- sprintf('Trend: Not significant\n p-value: %.2f\n Estimated slope: %.2f', mod_mk$sl, mod_sens$estimates)}
  plot_mod <- ggplot(analysis_df, aes(x = yr, y = deficit, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Annual Streamflow (mm)", title = paste(scenario, "Annual Modeled Streamflow"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$deficit), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_mod)
  plot_list[[i]] <- plot_mod
}
jpeg(file=paste0(outLocationPathFuture, "/", "Annual_Deficit_Trends.jpg"), width=300*num_models, height=200*num_models)
grid.arrange(grobs = plot_list, ncol=num_models/2, nrow=num_models/2)
dev.off()
