# ---------------------------------------------------------------------
# This script includes code to generate future projections of streamflow using the IHACRES 
# rainfall-streamflow methodology. This requires future projections of adjusted runoff,
# which are either calculated using pre-calibrated coefficients or pulled from a pre-generated
# gridded water balance model produced/maintained by Mike Tercek. This code also provides 
# preliminary visualizations of these future streamflow projections. 
# 
# EDITS IN PROGRESS
# implement code to pull and use Mike's data (see DataFunctions.R script)
# ---------------------------------------------------------------------



#######################################################################
### GENERATE FUTURE WATER BALANCE MODELS ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')

### Use Mike Tercek's pre-generated gridded CONUS water balance model for future projections ###
# EDIT
future_wb_conus <- get_conus_wb(SiteID_FileName, lat, lon, endY, 2099)
future_wb_conus$adj_runoff <-get_adj_runoff(future_wb_conus$runoff, gw_add = gw_add, vfm = vfm)

# not sure what workflow should be if they already have this file (i.e. directly from Mike)
# if(file.exists(file.path(dataPath, paste("WB",SiteID_FileName,"2023_2100.csv", sep = "_")))){
#   future_wb <- read.csv(file.path(dataPath, paste("WB",SiteID_FileName,"2023_2100.csv", sep = "_"))) 
#   future_wb$Date <- as.Date(future_wb$Date, '%m/%d/%Y')
#   future_wb<-subset(future_wb, projection !="MIROC-ESM-CHEM.rcp85")   # drop "MIROC-ESM-CHEM.rcp85" because it doesn't have an associated RCP 4.5
#   
#   #convert to mm and fix column names
#   future_wb<-cbind(future_wb[,c("Date","projection")], 25.4*(future_wb[,c(which(colnames(future_wb)=="Deficit.in"):ncol(future_wb))]))
#   colnames(future_wb)<- c("Date", "projection", "Deficit", "AET", "soil_water", "runoff", "rain", "accumswe", "PET")
#   
#   #adjust for ground water addition and volume forcing multiplier
#   future_wb$adj_runoff<- get_adj_runoff(future_wb$runoff, gw_add = gw_add, vfm = vfm)



### Re-run water balance model to generate future projections ###
# has the name projection been changed?
if(!file.exists(here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, endY, "2100.csv", sep='_')))){
  # Get future climate data
  future_climate <- get_maca_point(lat, lon, SiteID_FileName)
  future_climate$date <- as.Date(future_climate$date)

  # Run water balance code for each future projection
  future_wb_calc <- NULL
  for(projection in unique(future_climate$projection)){
    ClimData <- future_climate %>% filter(projection==projection) %>% select(-projection)
    DailyWB_future <- WB(ClimData, gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope,
                         shade.coeff, jtemp, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
    DailyWB_future <- cbind(projection, DailyWB_future)  
    future_wb_calc <- rbind(future_wb_calc, DailyWB_future)
  }
  # Save calculated WB
  future_wb_calc <- future_wb_calc %>% rename(projection = run)
  write.csv(future_wb_calc, here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, endY, "2100.csv", sep='_')), row.names=FALSE)
} else{
  # Read in calculated WB 
  future_wb_calc <- read.csv(here('Data',SiteID_FileName,paste('WB_calc',SiteID_FileName, endY, "2100.csv", sep='_')))
  future_wb_calc$date <- as.Date(future_wb_calc$date, '%m/%d/%Y')
}



### Compare the two future water balance projections, just for fun ###
# plot AET, deficit, adj runoff for a sample year and a sample model
model_run = 'HadGEM2-CC365.rcp45'; yr = 2060
if(make_plots){
  plot_aet <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(Date)==yr), aes(x=Date, y=AET), col='black')+
    geom_line(data=future_wb_calc %>% filter(projection==model_run & year(Date)==yr), aes(x=Date, y=AET), col='red')+
    labs(x='Date',y='AET [mm]', title='Actual Evapotranspiration') +
    theme(legend.position = "none") + nps_theme()
  plot_d <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(Date)==yr), aes(x=Date, y=Deficit), col='black')+
    geom_line(data=future_wb_calc %>% filter(projection==model_run & year(Date)==yr), aes(x=Date, y=D), col='red')+
    labs(x='Date',y='Deficit [mm]', title='Deficit') +
    theme(legend.position = "none") + nps_theme()
  plot_run <- ggplot() + geom_line(data=future_wb_conus %>% filter(projection==model_run & year(Date)==yr), aes(x=Date, y=adj_runoff, col='Gridded WB'))+
    geom_line(data=future_wb_calc %>% filter(projection==model_run& year(Date)==yr), aes(x=Date, y=adj_runoff, col='Calculated WB'))+
    labs(x='Date',y='Adjusted Runoff [mm]', title='Adjusted Runoff') +  
    scale_color_manual(values = c("Gridded WB"="black", "Calculated WB"="red"), name="WB Projections") + nps_theme()
  
  nameReduce = gsub(pattern = " ",replacement = "_", x = paste(SiteID, "Future WB Projection Comparison"))
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=2000, height=600)
  grid.arrange(plot_aet, plot_d, plot_run, ncol = 3, widths=c(1,1,1.3), top = textGrob(paste('WB Projection Comparisons for', model_run, ':', yr),gp=gpar(fontsize=30)))
  dev.off() 
}



#######################################################################
### GENERATE FUTURE STREAMFLOW PROJECTIONS ###

### Define which future water balance to use ###
if(calcFutureWB){
  future_wb <- future_wb_calc
} else{future_wb <- future_wb_conus}

### Run IHACRES model for each future projection ###
gcms<-unique(future_wb$projection)
futures <- NULL
for (j in 1:length(gcms)){
  # Subset one model
  fut_ro <- subset(future_wb, projection == gcms[j]); print(gcms[j])
  data <- data.frame(fut_ro$adj_runoff)
  colnames(data) <- c("adj_runoff")
  
  # Run IHACRES model
  DailyDrainFuture <- Drain(data, q0, qa, qb, s0, sa, sb, v0, va, vb)
  
  # Save streamflow projection to futures dataframe
  drainage_qsvt <- cbind(gcms[j], fut_ro$Date, DailyDrainFuture)
  colnames(drainage_qsvt)[] <- c("projection", "date", "adj_runoff", "quick", "slow", "veryslow", "total")
  futures <-rbind(futures, drainage_qsvt)
}

# Filter out futures that overlap with the date of the historic flow, extract GCM and RCP
futures<- futures[futures$date>endDate,]
futures$gcm <- sapply(X = strsplit(futures$projection, split=".rcp"), FUN = "[", 1) 
futures$rcp <- as.numeric(x = sapply(strsplit(futures$projection, split=".rcp"), FUN = "[", 2)) 



#######################################################################
### COMPILE FUTURE AND HISTORICAL PROJECTIONS ###

# Get historical model projections
hist_flow_mod <- data.frame(date=DailyDrain$date, adj_runoff=DailyDrain$adj_runoff, quick=DailyDrain$Quick,
                            slow=DailyDrain$Slow, veryslow=DailyDrain$Very_Slow, total=DailyDrain$total, 
                            projection = rep('Historical', length(DailyDrain$date)), gcm=rep('Historical', length(DailyDrain$date)), rcp=rep('Hist', length(DailyDrain$date)))


# Combine historical and future model projections
daily_df<-rbind(futures, hist_flow_mod)
daily_df$date <- as.Date(daily_df$date); daily_df$yr<-as.numeric(format(daily_df$date,"%Y")); daily_df$mo<-format(daily_df$date,"%m"); daily_df$yr_mo<-format(daily_df$date,"%Y-%m")
daily_df<-daily_df[,c("date", "projection", "gcm", "rcp", "yr", "mo", "yr_mo", "adj_runoff", "quick", "slow", "veryslow", "total")]
daily_df$Period<-ifelse(daily_df$yr<=2022,"Historical",ifelse (daily_df$yr>=2023 & daily_df$yr<=2050,"Early",
                                           ifelse (daily_df$yr>=2051 & daily_df$yr<=2070,"Middle", ifelse (daily_df$yr>=2071, "Late","NA"))))

# Aggregate data to annual
annual_df <- as.data.frame(daily_df %>% group_by(gcm, rcp, yr) %>%
  dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                   slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE), total = sum(total, na.rm = TRUE), 
                   Period=first(Period)))
annual_df$date<-as.Date(paste(annual_df$yr,"-01", "-01",sep=""))

# Aggregate data to monthly 
monthly_df <- as.data.frame(daily_df %>% group_by(gcm,rcp, yr_mo) %>%
  dplyr::summarize(projection=first(projection), gcm=first(gcm), rcp=first(rcp), adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                   slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE) ,total = sum(total, na.rm = TRUE),
                   Period=first(Period)))
monthly_df$date<-as.Date(paste(monthly_df$yr_mo,"-01",sep=""))

# Mean of future daily streamflow projections
mean_daily_df <- daily_df %>% filter(projection!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total))



#######################################################################
### SELECT DIVERGENT CLIMATE FUTURES ###
cf_names <- c('Warm Wet', 'Hot Dry', 'Warm Dry', 'Hot Wet')

# Pull meteorological data and group by year
# historical
hist_climate_ann <- as.data.frame(DailyClimData %>% group_by(year(date)) %>%
                                    dplyr::summarize(year=first(year(date)), pr = sum(pr, na.rm = TRUE),
                                                     tmmn = mean(tmmn, na.rm = TRUE), tmmx = mean(tmmx, na.rm = TRUE)))
hist_climate_ann$t_avg <- (hist_climate_ann$tmmn + hist_climate_ann$tmmx) / 2

# future
future_climate <- get_maca_point(lat, lon, SiteID_FileName); future_climate$date <- as.Date(future_climate$date, '%d/%m/%Y')
future_climate_ann <- as.data.frame(future_climate %>% group_by(projection, year(date)) %>%
                                                  dplyr::summarize(projection=first(projection), year=first(year(date)), pr = sum(pr, na.rm = TRUE),
                                                                   tmmn = mean(tmmn, na.rm = TRUE), tmmx = mean(tmmx, na.rm = TRUE)))
future_climate_ann$t_avg <- (future_climate_ann$tmmn + future_climate_ann$tmmx) / 2

# calculate averages over historical period
hist_avg_precip <- mean(hist_climate_ann$pr)
hist_avg_tavg <- mean(hist_climate_ann$t_avg)

# calculate averages over future period (30 yr average centered around 2050: 2035-2065)
# something is wrong: mismatch in magnitude between historical + future precip
future_climate_plot <- future_climate_ann %>% filter(year >= 2035 & year <= 2065) %>% group_by(projection) %>% 
  dplyr::summarize(pr=mean(pr, na.rm=TRUE), t_avg=mean(t_avg, na.rm=TRUE), pr_delta=mean(pr, na.rm=TRUE)-hist_avg_precip,
                   tavg_delta=mean(t_avg, na.rm=TRUE)-hist_avg_tavg)
future_centroid_pr <- mean(future_climate_plot$pr_delta)
future_centroid_tavg <- mean(future_climate_plot$tavg_delta)

# Create individual columns for gcm/rcp
future_climate_plot$gcm <- sapply(strsplit(future_climate_plot$projection, split = "\\."), `[`, 1)
future_climate_plot$rcp <- sapply(strsplit(future_climate_plot$projection, split = "\\."), `[`, 2)

# Identify the furthest futures in each quadrant
### EDIT THIS ###
future_climate_plot$distance_from_pr <- future_climate_plot$pr_delta-future_centroid_pr 
future_climate_plot$distance_from_tavg <- future_climate_plot$tavg_delta-future_centroid_tavg

# plot for visualization
plot <- ggplot(data=future_climate_plot, aes(x=tavg_delta, y=pr_delta, color=rcp)) + geom_point() +
  geom_text_repel(aes(label = gcm), color = 'black', max.overlaps=Inf) +
  geom_hline(aes(yintercept=quantile(pr_delta, 0.5)), color = "black", linetype='dashed') + geom_vline(aes(xintercept=quantile(tavg_delta, 0.5)), color = "black", linetype='dashed') +
  geom_rect(aes(xmin = quantile(tavg_delta, 0.25), xmax = quantile(tavg_delta, 0.75), ymin = quantile(pr_delta, 0.25), ymax = quantile(pr_delta, 0.75)), color = "black", size=1, alpha=0) +
  labs(title=paste('Changes in climate means by 2050 at',SiteID), x='Change in annual average temperature [C]', y='Change in annual average precipitation [mm]', color='RCP') + 
  scale_color_manual(values = c("rcp45" = "orange", "rcp85" = "red")) + nps_theme()
print(plot)
dev.off()



#######################################################################
### PLOT FUTURE STREAMFLOW ### 

if(make_plots){
  # Daily streamflow projections for all models
  name = paste(SiteID, "Daily Streamflow")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=1300, height=800)
  plot<- ggplot() + geom_line(data=daily_df %>% filter(rcp!="Hist"), aes(x = date, y = total, color = 'Future')) + 
    facet_wrap(~projection, ncol=6)+ ylab("Streamflow (mm)") + xlab("Year")+ ggtitle(name) + nps_theme() +
    geom_line(data=daily_df %>% filter(rcp=="Hist")%>% select(-projection), aes(x = date, y = total, color=Period)) + 
    scale_color_manual(values = c( "Historical" = "blue", "Future" = "black"), limits=c('Historical','Future'), name = "Period")
  print(plot)
  dev.off()
  
  # Monthly streamflow projections for all models
  name = paste(SiteID, "Monthly Streamflow")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=1300, height=800)
  plot<- ggplot() + geom_line(data=monthly_df %>% filter(rcp!="Hist"), aes(x = date, y = total, color = 'Future')) + 
    facet_wrap(~projection, ncol=6)+ ylab("Streamflow (mm)") + xlab("Year")+ ggtitle(name) + nps_theme() +
    geom_line(data=monthly_df %>% filter(rcp=="Hist")%>% select(-projection), aes(x = date, y = total, color=Period)) + 
    scale_color_manual(values = c( "Historical" = "blue", "Future" = "black"), limits=c('Historical','Future'), name = "Period")
  print(plot)
  dev.off()
  
  # Annual streamflow projections for all models
  name = paste(SiteID, "Annual Streamflow")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=1300, height=800)
  plot<- ggplot() + geom_line(data=annual_df %>% filter(rcp!="Hist"), aes(x = date, y = total, color = 'Future')) + 
    facet_wrap(~projection, ncol=6)+ ylab("Streamflow (mm)") + xlab("Year")+ ggtitle(name) + nps_theme() +
    geom_line(data=annual_df %>% filter(rcp=="Hist")%>% select(-projection), aes(x = date, y = total, color=Period)) + 
    scale_color_manual(values = c( "Historical" = "blue", "Future" = "black"), limits=c('Historical','Future'), name = "Period")
  print(plot)
  dev.off()
  

  # Time series of daily, monthly, annual streamflow
  plot_daily <- ggplot() + labs(title='Daily Modeled Streamflow', x='Year', y='Streamflow [mm]', color='Model') + nps_theme() + 
    geom_line(data=(daily_df %>% filter(projection!='Historical')), aes(x=date, y=total), col='gray', alpha = 0.7) + 
    geom_line(data=(daily_df %>% filter(projection=='Historical')), aes(x=date, y=total), col='blue',alpha=1, linewidth=1.5) +
    geom_line(data=daily_df %>% filter(projection!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total), col='black', alpha=1, linewidth=1.5) + 
    guides(alpha = "none") + theme(legend.position="none")
  plot_mon <-  ggplot() + labs(title='Monthly Modeled Streamflow', x='Year', y='Streamflow [mm]', color='Model') + nps_theme() + 
    geom_line(data=(monthly_df %>% filter(projection!='Historical')), aes(x=date, y=total), col='gray', alpha = 0.7) + 
    geom_line(data=(monthly_df %>% filter(projection=='Historical')), aes(x=date, y=total), col='blue', alpha=1, linewidth=1.5) +
    geom_line(data=monthly_df %>% filter(projection!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total), col='black',alpha=1, linewidth=1.5) + 
    guides(alpha = "none") + theme(legend.position="none")
  plot_ann <-  ggplot() + labs(title='Annual Modeled Streamflow', x='Year', y='Streamflow [mm]', color='Model') + nps_theme() + 
    geom_line(data=(annual_df %>% filter(projection!='Historical')), aes(x=date, y=total, col='Individual Model Projections'), alpha = 0.7) + 
    geom_line(data=(annual_df %>% filter(projection=='Historical')), aes(x=date, y=total, col='Historical'), alpha=1, linewidth=1.5) +
    geom_line(data=annual_df %>% filter(projection!='Historical') %>% group_by(date) %>% dplyr::summarize(mean_total=mean(total)), 
              aes(x=date, y=mean_total, col='Mean Model Projections'), alpha=1, linewidth=1.5) + 
    scale_color_manual(values = c('Historical'='blue', "Individual Model Projections" = "gray", 'Mean Model Projections'="black"),
                       labels= c(expression('Historical', "Individual\nModel Projections"), expression('Mean Model\nProjections'))) + 
    guides(alpha = "none") 
  nameReduce = gsub(pattern = " ",replacement = "_", x = paste(SiteID, "Streamflow Projections Time Series"))
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=2000, height=600)
  grid.arrange(plot_daily, plot_mon, plot_ann, ncol = 3, widths=c(1,1,1.3))
  dev.off()
  
  # Annual streamflow projections for all models with trends
  name = paste(SiteID, "Annual Streamflow Trends")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=2600, height=1600)
  plot<- ggplot() + geom_line(data=annual_df %>% filter(rcp!="Hist"), aes(x = date, y = total, color = 'Future')) +
    geom_smooth(data=annual_df %>% filter(rcp!="Hist"), aes(x = date, y = total, color = 'Future'), method = "lm", se = FALSE) + 
    facet_wrap(~projection)+ ylab("Streamflow (mm)") + xlab("Year")+ ggtitle(name) + nps_theme() +
    geom_line(data=annual_df %>% filter(rcp=="Hist")%>% select(-projection), aes(x = date, y = total, color=Period)) + 
    geom_smooth(data=annual_df %>% filter(rcp=="Hist")%>% select(-projection), aes(x = date, y = total, color=Period), method = "lm", se = FALSE) + 
    scale_color_manual(values = c( "Historical" = "blue", "Future" = "black"), limits=c('Historical','Future'), name = "Period")
  print(plot)
  dev.off()
  
  
  # Time series of all the models together, with user-specified individual models in bold
  name = paste(SiteID, "Annual Streamflow Projections Time Series")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=1300, height=800)
  plot <- ggplot() + labs(title=paste(SiteID, 'Annual Streamflow Projections Time Series'), y='Streamflow [mm]', x='Years', color='Model') +  
    geom_line(data=(annual_df %>% filter(projection=='Historical')), aes(x=yr, y=total, group=projection, color=ifelse(projection=='Historical', 'Historical', projection)), alpha=1, linewidth=1.5) + 
    geom_line(data=(annual_df %>% filter(!projection %in% individual_models)), aes(x=yr, y=total, group = projection, color=ifelse(projection %in% individual_models, projection, "Other")), alpha = 0.7) + 
    geom_line(data=(annual_df %>% filter(projection %in% individual_models)), aes(x=yr, y=total, group=projection, color=ifelse(projection %in% individual_models, projection, "Other")), alpha=1, linewidth=1.5) +
    scale_color_manual(values = c("Other" = "gray", 'Historical'='blue', setNames(sample(colors(), length(individual_models)), individual_models))) +
    guides(alpha = "none") + nps_theme()
  print(plot)
  dev.off()
  
  
  # Climate future scatterplot: annual magnitude vs daily standard deviation
  name = paste(SiteID, "Streamflow Climate Future Scatterplot")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  annual_df$delta_annual_mm <- annual_df$total - mean((annual_df %>% filter(gcm=='Historical'))$total)
  delta_plot <- daily_df %>% group_by(gcm, rcp, yr) %>% 
    dplyr::summarize(Period=first(Period), delta_daily_sd = sd(total)-sd((daily_df %>% filter(gcm=='Historical'))$total))
  delta_plot <- delta_plot %>% left_join(annual_df %>% select(gcm, rcp, yr, delta_annual_mm, Period), 
                                         by = c("gcm", "rcp", "yr","Period"))
  delta_plot <- delta_plot %>% filter(gcm!='Historical') %>% group_by(gcm,rcp,Period) %>% dplyr::summarise(delta_daily_sd=mean(delta_daily_sd), delta_annual_mm=mean(delta_annual_mm))
  annual_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                                xintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_annual_mm, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_annual_mm, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_annual_mm, 0.5)))
  sd_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                            yintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_daily_sd, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_daily_sd, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_daily_sd, 0.5)))
  annual_zero <- data.frame(Period = c("Early", "Middle", "Late"), xintercept=c(0,0,0)); sd_zero <- data.frame(Period = c("Early", "Middle", "Late"), yintercept = c(0,0,0))
  
  jpeg(file=paste0(outLocationPath, "/", nameReduce, ".jpg"), width=1000, height=400)
  plot <- ggplot(data=delta_plot, aes(x=delta_annual_mm,y=delta_daily_sd,color=rcp)) + geom_point() +
    geom_text_repel(aes(label = gcm), color = 'black', max.overlaps=Inf) +
    geom_hline(data = sd_quantile, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_quantile, aes(xintercept = xintercept), color = "black") +
    #geom_hline(data = sd_zero, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_zero, aes(xintercept = xintercept), color = "black") +
    facet_wrap(~factor(Period, levels = c('Early', 'Middle', 'Late'))) + labs(title=paste('Changes in streamflow at',SiteID), x='Change in annual magnitude [mm]', y='Change in daily standard deviation [mm]', color='RCP') + 
    scale_color_manual(values = c("45" = "orange", "85" = "red")) + nps_theme()
  print(plot)
  dev.off()
}


# trend plots (like historical)
model1 <- ''; model2 <- ''
model_names <- c('Warm Wet', 'Hot Dry')


# heatmap
daily_df$day <- factor(yday(daily_df$date), levels = unique(yday(daily_df$date)))
jpeg(file=paste0(outLocationPathHist, "/", "Modeled_Daily_Heatmap.jpg"), width=600, height=400)
ggplot(daily_df, aes(day, yr, fill = as.numeric(total))) + geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), trans='log', breaks=c(min(daily_df$total), 0.01, 0.1, 10, 100, max(daily_df$total)), 
                       labels=c(sprintf('%.3f',  min(daily_df$total)),'.01','0.1','10', '100', sprintf('%.0f',  max(daily_df$total)))) +
  labs(title = "Modeled daily streamflow", x='Month', fill = "Streamflow [mm]") +
  scale_x_discrete(breaks=c("1", "32", "60", "91", "121", "152", "182", "213", "244", "274", "305", "335"), 
                   labels = c("1" = "January", "32" = "February", "60" = "March", "91" = "April", "121" = "May", "152" = "June", 
                              "182" = "July", "213" = "August", "244" = "September", "274" = "October", "305" = "November", "335" = "December")) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()



