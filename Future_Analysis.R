# ---------------------------------------------------------------------
# This script includes code to generate future projections of streamflow using the IHACRES 
# rainfall-streamflow methodology. This requires future projections of adjusted runoff,
# which are either calculated using pre-calibrated coefficients or pulled from a pre-generated
# gridded water balance model produced/maintained by Mike Tercek. This code also provides 
# preliminary visualizations of these future streamflow projections. 
# 
# EDITS IN PROGRESS
# implement code to pull and use Mike's data (see DataFunctions.R script)
# add future scenarios to plots of the future
# remember to upload the model performance file to GitHub
# make sure everything runs smoothly
# figure out what's the appropriate line of reasoning for the future wb models
# ---------------------------------------------------------------------



#######################################################################
### GENERATE FUTURE WATER BALANCE MODELS ###
gcm_list <- c('BNU-ESM', 'CCSM4', 'CNRM-CM5', 'CSIRO-Mk3-6-0', 'CanESM2','GFDL-ESM2G', 'HadGEM2-CC365', 
              'IPSL-CM5A-LR', 'MIROC5', 'MIROC-ESM-CHEM','MRI-CGCM3', 'NorESM1-M', 'inmcm4')

# Remove low skill models using list from Rupp et al. 2016
low_skill_models = read.delim('./Data/GCM_skill_by_region.txt', header=TRUE) %>% 
  filter(Region == ifelse(region %in% Region, region, "mean")) %>% top_n(n=round(length(gcm_list)*percent_skill_cutoff), wt=Rank)
#gcm_list[!gcm_list %in% low_skill_models$GCM]

### Use Mike Tercek's pre-generated gridded CONUS water balance model for future projections ###
future_wb_conus <- get_conus_wb(SiteID_FileName, lat, lon, endY, 2099)
future_wb_conus$adj_runoff <-get_adj_runoff(future_wb_conus$runoff, gw_add = gw_add, vfm = vfm)

# This is the code to read in the file if it was provided directly by Mike Tercek
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

# Define which future water balance projection to use
if(calcFutureWB){
  future_wb <- future_wb_calc
} else{future_wb <- future_wb_conus}

# Run IHACRES model for each future projection
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
daily_df$date <- as.Date(daily_df$date); daily_df$yr<-as.numeric(format(daily_df$date,"%Y")); daily_df$mo<-format(daily_df$date,"%m"); daily_df$yr_mo<-format(daily_df$date,"%Y-%m"); daily_df$day <- factor(yday(daily_df$date), levels = unique(yday(daily_df$date)))
daily_df$water_year <- sapply(daily_df$date, get_water_year); 
daily_df <- daily_df %>% group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
daily_df<-daily_df[,c("date", "projection", "gcm", "rcp", "yr", "mo", "yr_mo", "water_year", "adj_runoff", "quick", "slow", "veryslow", "total")]
daily_df$Period<-ifelse(daily_df$yr<=2022,"Historical",ifelse (daily_df$yr>=2023 & daily_df$yr<=2050,"Early",
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
cf_names <- c("Warm Wet", "Hot Wet", "Central", "Warm Dry", "Hot Dry")

### Pull and aggregate meteorological data ###
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
# probably an error in get_maca_point :(
future_means <- future_climate_ann %>% filter(year >= 2035 & year <= 2065) %>% group_by(projection) %>% 
  dplyr::summarize(pr=mean(pr, na.rm=TRUE), t_avg=mean(t_avg, na.rm=TRUE), pr_delta=mean(pr, na.rm=TRUE)-hist_avg_precip,
                   tavg_delta=mean(t_avg, na.rm=TRUE)-hist_avg_tavg)
future_centroid_pr <- mean(future_means$pr_delta)
future_centroid_tavg <- mean(future_means$tavg_delta)

# Create individual columns for gcm/rcp
future_means$gcm <- sapply(strsplit(future_means$projection, split = "\\."), `[`, 1)
future_means$rcp <- sapply(strsplit(future_means$projection, split = "\\."), `[`, 2)


### Identify the furthest futures in each quadrant using principal component analysis (PCA) ###
# Adapted from Amber's climate futures code: https://github.com/nationalparkservice/CCRP_automated_climate_futures/blob/master/scripts/Plot_Table_Creation.R

### Label each climate future with quadrants
Pr0 = as.numeric(quantile(future_means$pr_delta, 0)); Pr25 = as.numeric(quantile(future_means$pr_delta, 0.25)); PrAvg = as.numeric(mean(future_means$pr_delta)); Pr75 = as.numeric(quantile(future_means$pr_delta, 0.75)); Pr100 = as.numeric(quantile(future_means$pr_delta, 1))
Tavg0 = as.numeric(quantile(future_means$tavg_delta, 0)); Tavg25 = as.numeric(quantile(future_means$tavg_delta, 0.25)) ; Tavg = as.numeric(mean(future_means$tavg_delta)); Tavg75 = as.numeric(quantile(future_means$tavg_delta, 0.75)); Tavg100 = as.numeric(quantile(future_means$tavg_delta, 1))

future_means$CF1 = as.numeric((future_means$tavg_delta<Tavg & future_means$pr_delta>Pr75) | future_means$tavg_delta<Tavg25 & future_means$pr_delta>PrAvg)
future_means$CF2 = as.numeric((future_means$tavg_delta>Tavg & future_means$pr_delta>Pr75) | future_means$tavg_delta>Tavg75 & future_means$pr_delta>PrAvg)
future_means$CF3 = as.numeric((future_means$tavg_delta>Tavg25 & future_means$tavg_delta<Tavg75) & (future_means$pr_delta>Pr25 & future_means$pr_delta<Pr75))
future_means$CF4 = as.numeric((future_means$tavg_delta<Tavg & future_means$pr_delta<Pr25) | future_means$tavg_delta<Tavg25 & future_means$pr_delta<PrAvg)
future_means$CF5 = as.numeric((future_means$tavg_delta>Tavg & future_means$pr_delta<Pr25) | future_means$tavg_delta>Tavg75 & future_means$pr_delta<PrAvg)


#Assign full name of climate future to new variable CF
future_means$CF <- NULL
for(i in 1:5) {
  future_means$CF[future_means[[paste0("CF", i)]] == 1] <- cf_names[i]
}
future_means <- future_means %>% select(-CF1, -CF2, -CF3, -CF4, -CF5)

### Selection with corners method
# Identify corners
lx = min(future_means$tavg_delta); ux = max(future_means$tavg_delta); ly = min(future_means$pr_delta); uy = max(future_means$pr_delta)
ww = c(lx,uy); wd = c(lx,ly); hw = c(ux,uy); hd = c(ux,ly)

# Calculate Euclidean distance of each point from corners
pts <- future_means
pts$WW.distance <- sqrt((pts$tavg_delta - ww[1])^2 + (pts$pr_delta - ww[2])^2); pts$WD.distance <- sqrt((pts$tavg_delta - wd[1])^2 + (pts$pr_delta - wd[2])^2)
pts$HW.distance <- sqrt((pts$tavg_delta - hw[1])^2 + (pts$pr_delta - hw[2])^2); pts$HD.distance <- sqrt((pts$tavg_delta - hd[1])^2 + (pts$pr_delta - hd[2])^2)

# Select scenarios based on shortest distance
pts %>% filter(CF == "Warm Wet") %>% slice(which.min(WW.distance)) %>% .$projection -> ww
pts %>% filter(CF == "Warm Dry") %>% slice(which.min(WD.distance)) %>% .$projection -> wd
pts %>% filter(CF == "Hot Wet") %>% slice(which.min(HW.distance)) %>% .$projection -> hw
pts %>% filter(CF == "Hot Dry") %>% slice(which.min(HD.distance)) %>% .$projection -> hd

future_means %>% mutate(corners = ifelse(projection == ww,"Warm Wet",
                                         ifelse(projection == wd, "Warm Dry",
                                                ifelse(projection == hw, "Hot Wet",
                                                       ifelse(projection == hd, "Hot Dry",NA))))) -> future_means

### PCA
FM <- future_means %>% select("projection","pr_delta","tavg_delta") %>%  
  remove_rownames %>% column_to_rownames(var="projection") 
CF_GCM = data.frame(projection = future_means$projection, CF = future_means$CF)
pca.df <- as.data.frame(prcomp(FM, center = TRUE,scale. = TRUE)$x)

# Save results of PCA
if(make_plots){
  ggsave("PCA-loadings.jpg", plot=autoplot(pca, data = FM, loadings = TRUE,label=TRUE), width=8, height=5, path = outLocationPath) 
}
write.csv(pca.df, paste0(outLocationPath, "/PCA-loadings.csv"))

#Take the min/max of each of the PCs
PCs <-rbind(data.frame(projection = c(rownames(pca.df)[which.min(pca.df$PC1)],rownames(pca.df)[which.max(pca.df$PC1)]),PC="PC1"),
            data.frame(projection = c(rownames(pca.df)[which.min(pca.df$PC2)],rownames(pca.df)[which.max(pca.df$PC2)]),PC="PC2"))

# Assigns CFs to diagonals
diagonals <- rbind(data.frame(CF = cf_names[c(1,5)],diagonals=factor("diagonal1")),data.frame(CF = cf_names[c(4,2)],diagonals=factor("diagonal2")))

# Aggregate and add PCA selections to dataframes
PCA <- CF_GCM %>% filter(projection %in% PCs$projection) %>% left_join(diagonals,by="CF") %>% right_join(PCs,by="projection")
future_means %>% mutate(pca = ifelse(projection %in% PCs$projection[which(PCs$PC=="PC1")], as.character(CF), #assign PCs to quadrants and select those projections
                                     ifelse(projection %in% PCs$projection[which(PCs$PC=="PC2")], as.character(CF),NA))) -> future_means #future_means

# Check whether duplicate models were selected and fix if so
if(length(setdiff(cf_names[cf_names != "Central"],future_means$pca)) > 0){ #if a quadrant is missing 
  future_means$pca[which(future_means$corners == setdiff(cf_names[cf_names != "Central"],future_means$pca))] = setdiff(cf_names[cf_names != "Central"],future_means$pca) #assign corners selection to that CF
  if(nrow(PCA[duplicated(PCA$projection),]) > 0) { #If there is a redundant projection
    future_means$pca = future_means$pca #Do nothing - otherwise end up with empty quadrant. This line could be removed and make the previous statement inverse but it makes it more confusing what's going on that way
  } else{
    future_means$pca[which(future_means$projection == ID.redundant.gcm(PCA))] = NA #Removes the projection that is in redundant diagonal
  }
}


### Plot for visualization ###
plot <- ggplot(data=future_means, aes(x=tavg_delta, y=pr_delta, color=rcp)) + geom_point() +
  geom_text_repel(aes(label = gcm), color = 'black', max.overlaps=Inf) +
  geom_hline(aes(yintercept=PrAvg), color = "black", linetype='dashed') + geom_vline(aes(xintercept=Tavg), color = "black", linetype='dashed') +
  geom_rect(aes(xmin = Tavg25, xmax = Tavg75, ymin = Pr25, ymax = Pr75), color = "black", linewidth=1, alpha=0) +
  labs(title=paste('Changes in climate means by 2050 at',SiteID), x='Change in annual average temperature [C]', y='Change in annual average precipitation [mm]', color='RCP') + 
  scale_color_manual(values = c("rcp45" = "orange", "rcp85" = "red")) + nps_theme()
print(plot)
dev.off()


# Identify models in format for plotting
model_names <- (future_means %>% filter(!is.na(pca)))$projection; scenario_names <- (future_means %>% filter(!is.na(pca)))$pca
color_names <- c("#12045C","#E10720","#9A9EE5","#F3D3CB")

# identify warm wet/hot dry [1:2] or warm dry/hot wet [3:4] or all (comment out both lines)
model_names <- model_names[1:2]; scenario_names <- scenario_names[1:2]; color_names <- color_names[1:2]
#model_names <- model_names[3:4]; scenario_names <- scenario_names[3:4]; color_names <- color_names[3:4]


#######################################################################
### PLOT FUTURE STREAMFLOW FOR ALL MODELS ### 

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



#######################################################################
### PLOT FUTURE STREAMFLOW FOR IDENTIFIED MODELS ###

### Heatmap - daily
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- daily_df %>% filter(projection=='Historical' | projection==proj)
  plot <- ggplot(analysis_df, aes(factor(water_day), water_year, fill = as.numeric(total))) + geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), trans='log', breaks=c(min(daily_df$total), 0.01, 0.1, 10, 100, max(daily_df$total)), 
                         labels=c(sprintf('%.3f',  min(daily_df$total)),'.01','0.1','10', '100', sprintf('%.0f', max(daily_df$total))),
                         limits = c(0.00000001, max(daily_df$total, na.rm = TRUE))) +
    labs(title = paste0(scenario, " Daily Streamflow (", proj, ")"), x='Month', y='Water Year', fill = "Streamflow [mm]") +
    scale_x_discrete(breaks=c("1", "32", "63", "94", "120", "151", "181", "212", "243", "274", "305", "335"), 
                     labels = c("1" = "October", "32" = "November", "63" = "December", "94" = "January", "120" = "February", "151" = "March", 
                                "181" = "April", "212" = "May", "243" = "June", "274" = "July", "305" = "August", "335" = "September")) +
    theme(axis.text.x = element_text(angle = 90))
  plot_list[[i]] <- plot
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Daily_Heatmap.jpg"), width=600, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()


### Heatmap - monthly 
# order by water year?
monthly_df$month <- month(monthly_df$date); monthly_df$year <- year(monthly_df$date)

plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  plot <- ggplot(monthly_df, aes(factor(month), year, fill = total)) + geom_tile() +
    #scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
    scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), trans='log', breaks=c(0.1, 1, 10, 100, 1000, max(monthly_df$total)), 
                         labels=c('0.1', '1', '10', '100', '1000', sprintf('%.0f', max(monthly_df$total))),
                         limits = c(min(monthly_df$total, na.rm = TRUE), max(monthly_df$total, na.rm = TRUE))) + 
    labs(title = paste0(scenario, " Monthly Streamflow (", proj, ")"), x='Month', y='Water Year', fill = "Streamflow [mm]") +
    scale_x_discrete(labels = c("1" = "January", "2" = "February", "3" = "March", 
                                "4" = "April", "5" = "May", "6" = "June", 
                                "7" = "July", "8" = "August", "9" = "September", 
                                "10" = "October", "11" = "November", "12" = "December")) +
    nps_theme() + theme(axis.text.x = element_text(angle = 90))
  plot_list[[i]] <- plot
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Monthly_Heatmap.jpg"), width=600, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()






#######################################################################
### PLOT STREAMFLOW TRENDS AND METRICS ###
# edit width of jpg figure depending on model_names

# General Mann-Kendall test on daily streamflow ? should probably be seasonal
# this is too much data to provide a useful visualization - delete?
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]

  analysis_df <- daily_df %>% filter(projection=='Historical' | projection==proj)
  meas_mk <- SeasonalMannKendall(ts(analysis_df$total, start=c(year(startDate), 1), frequency=365))
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  plot_meas <- ggplot(analysis_df, aes(x = date, y = total, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Daily Streamflow (mm)", title = paste(scenario, "Daily Modeled Streamflow"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$date), y = max(analysis_df$total), label = label, color = "black", hjust = 1, vjust = 1) + scale_y_log10()
  print(plot_meas)
  plot_list[[i]] <- plot_meas
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Daily_Flow_Trends.jpg"), width=1200, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()

# Seasonal Mann-Kendall test on monthly streamflow
meas_mk <- SeasonalMannKendall(ts(monthly_df$total, start=c(year(startDate), 1), frequency=12))
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Monthly_Flow_Trends.jpg"), width=600, height=400)
plot_meas <- ggplot(monthly_df, aes(x = as.yearmon(yr_mo), y = total)) + geom_line(aes(color = 'Modeled'), na.rm=TRUE, linewidth=1) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
  labs(x = "Water Year", y = "Monthly Streamflow (mm)", title = "Monthly Modeled Streamflow", color='') +
  nps_theme() + theme(legend.position = 'bottom') +
  scale_color_manual(values = c("Modeled" = "black", "Trend" = "red")) +
  annotate("text", x = max(as.yearmon(monthly_df$yr_mo)), y = max(monthly_df$total), label = label, color = "black", hjust = 1, vjust = 1)
plot_meas
dev.off()


# Annual streamflow volume
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- annual_df %>% filter(projection=='Historical' | projection==proj)
  meas_mk <- MannKendall(analysis_df$total)
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  plot_meas <- ggplot(analysis_df, aes(x = yr, y = total, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Annual Streamflow (mm)", title = paste(scenario, "Annual Modeled Streamflow"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$total), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_meas)
  plot_list[[i]] <- plot_meas
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Annual_Flow_Trends.jpg"), width=1200, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()


# High flows (above 95%)
high_flow_q = 0.95
high_flow_mm = quantile((daily_df %>% filter(projection=='Historical'))$total, 0.95, na.rm=TRUE)

plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  # Calculate days above 95th percentile
  analysis_df <- as.data.frame(daily_df %>% filter(projection=='Historical' | projection==proj) %>% 
                                 mutate(high_flow = ifelse(total >= high_flow_mm, 1, 0)) %>%
                                 group_by(yr) %>% dplyr::summarize(projection=first(projection), days = sum(high_flow)))
  
  # Plot
  meas_mk <- MannKendall(analysis_df$days)
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  plot_meas <- ggplot(analysis_df, aes(x = yr, y = days, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Days per year", title = paste(scenario, "Days Above Historical 95th Percentile"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_meas)
  plot_list[[i]] <- plot_meas
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_High_Flow_Trends.jpg"), width=1200, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()



# Low flows (below 5%)
low_flow_q = 0.05
low_flow_mm = quantile(daily_df$total, 0.05, na.rm=TRUE)

plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  # Calculate days below 5th percentile
  analysis_df <- as.data.frame(daily_df %>% filter(projection=='Historical' | projection==proj) %>% 
                                 mutate(low_flow = ifelse(total <= low_flow_mm, 1, 0)) %>%
                                 group_by(yr) %>% dplyr::summarize(projection=first(projection), days = sum(low_flow)))
  
  # Plot
  meas_mk <- MannKendall(analysis_df$days)
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  plot_meas <- ggplot(analysis_df, aes(x = yr, y = days, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Days per year", title = paste(scenario, "Days Below Historical 5th Percentile"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$yr), y = max(analysis_df$days), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_meas)
  plot_list[[i]] <- plot_meas
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Low_Flow_Trends.jpg"), width=1200, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()


# 50% flow date
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  # Pull projections and add water day
  analysis_df <- daily_df %>% filter(projection=='Historical' | projection==proj) %>% 
    group_by(water_year) %>% mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))
  
  # Calculate center of timing (CT)
  analysis_df <- as.data.frame(analysis_df) %>% mutate(tq = total*water_day) %>%
                                 group_by(water_year) %>% dplyr::summarize(projection=first(projection), ct = sum(tq)/sum(total))

  # Plot
  meas_mk <- MannKendall(analysis_df$ct)
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  plot_meas <- ggplot(analysis_df, aes(x = water_year, y = ct, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Days after October 1", title = paste(scenario, "50% Flow Date"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$ct), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_meas)
  plot_list[[i]] <- plot_meas
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_50th_Flow_Trends.jpg"), width=1200, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()





# Calculate Q7
hist_q7 <- daily_df %>% mutate(Q7 = rollmean(total, 7, align="right", fill=NA)) %>%
  group_by(projection, water_year) %>% dplyr::summarize(min_q7 = ifelse(all(is.na(Q7)), NA, min(Q7, na.rm = TRUE)), 
                                           max_q7 = ifelse(all(is.na(Q7)), NA, max(Q7, na.rm = TRUE)), 
                                           avg_q7 = ifelse(all(is.na(Q7)), NA, mean(Q7, na.rm = TRUE)))

# Q7 min
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- hist_q7 %>% filter(projection=='Historical' | projection==proj)
  meas_mk <- MannKendall(analysis_df$min_q7)
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  plot_meas <- ggplot(analysis_df, aes(x = water_year, y = min_q7, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Streamflow [mm]", title = paste(scenario, "Min. 7 Day Flow (Q7 Min)"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$min_q7), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_meas)
  plot_list[[i]] <- plot_meas
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Q7Min_Trends.jpg"), width=1200, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()


# Q7 max
plot_list <- list()
for (i in 1:length(model_names)){
  proj = model_names[i]
  scenario <- scenario_names[i]
  
  analysis_df <- hist_q7 %>% filter(projection=='Historical' | projection==proj)
  meas_mk <- MannKendall(analysis_df$max_q7)
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  plot_meas <- ggplot(analysis_df, aes(x = water_year, y = max_q7, color=factor(projection))) + geom_line(na.rm=TRUE, linewidth=1, alpha=0.7) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend'), linetype='dashed', linewidth=1.5) +
    labs(x = "Water Year", y = "Streamflow [mm]", title = paste(scenario, "Max. 7 Day Flow (Q7 Max)"), color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c('Historical'='black', setNames(color_names, model_names), "Trend"="black")) +
    annotate("text", x = max(analysis_df$water_year), y = max(analysis_df$max_q7), label = label, color = "black", hjust = 1, vjust = 1) 
  print(plot_meas)
  plot_list[[i]] <- plot_meas
}
jpeg(file=paste0(outLocationPath, "/", "Modeled_Q7Max_Trends.jpg"), width=1200, height=400)
grid.arrange(grobs = plot_list, ncol=length(model_names))
dev.off()




#######################################################################
### Customizable plots: metric, season, etc ###
# EDIT TO PLOT SPECIFIC FUTURE SCENARIOS

# define flow level in cfs and convert to mm
flow_level <- 1800 #* 28316847*86400/(2590000000000*watershed_area) 
# identify months of interest (numerical values)
mos <- c(2, 3) 
# is comparison above or below threshold?
comparison = 'above'


# calculate number of days above/below threshold
if (tolower(comparison)=='below') {
  hist_threshold <- as.data.frame(DailyStream %>% mutate(flow = ifelse(CFS <= flow_level, 1, 0)) %>%
                                    group_by(waterYear) %>% dplyr::summarize(days = sum(flow)))# %>%filter(!is.na(days)))
} else if (tolower(comparison)=='above'){
  hist_threshold <- as.data.frame(DailyStream %>% mutate(flow = ifelse(CFS >= flow_level, 1, 0)) %>%
                                    group_by(waterYear) %>% dplyr::summarize(days = sum(flow)))# %>%filter(!is.na(days)))
}

# plot
meas_mk <- MannKendall(hist_threshold$days)
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
jpeg(file=paste0(outLocationPath, "/", "Days_", comparison, "_", flow_level, "_", month.abb[mos][1], "_", month.abb[mos][length(mos)],".jpg"), width=600, height=400)
plot_meas <- ggplot(hist_threshold, aes(x = waterYear, y = days)) + geom_line(aes(color = 'Modeled'), na.rm=TRUE, linewidth=1) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
  labs(x = "Water Year", y = "Days", title = paste('Days', comparison, flow_level, 'cfs, ', month.abb[mos][1], '-', month.abb[mos][length(mos)]), color='') +
  nps_theme() + theme(legend.position = 'bottom') +
  scale_color_manual(values = c("Modeled" = "black", "Trend" = "red")) +
  annotate("text", x = max(hist_threshold$waterYear), y = max(hist_threshold$days, na.rm=TRUE), label = label, color = "black", hjust = 1, vjust = 1)
plot_meas
dev.off()




