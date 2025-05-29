# ---------------------------------------------------------------------
# This script includes code to analyze historical measurements of streamflow. In addition to
# assessing historical trends, this script compares historical measurements to modeled
# streamflow, estimated using the IHACRES  rainfall-streamflow methodology.
# This code also provides preliminary analyses and visualizations of historical streamflow. 
# All trend analyses assume p < 0.05 is significant; this can be changed in the code
# 
# EDITS IN PROGRESS
# doesn't work when run in conjunction with main script
# make script so it can be run independently or in conjunction with main script
# figure out how to export 3D plot
# ---------------------------------------------------------------------

#######################################################################
# Source in function files
library('here')
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
# redefine gage site ID so this script can be run independently
GageSiteID <- GageSiteID
SiteID <- SiteID
#GageSiteID <- '03460000'
#SiteID <- 'Cataloochee'
#make_plots <- TRUE

SiteID_FileName <- gsub(pattern = " ", x = SiteID, replacement = "")

# create folder to store results of historical analysis
if(!dir.exists(here('Data', SiteID_FileName))) {dir.create(here('Data', SiteID_FileName))}; dataPath <- here('Data', SiteID_FileName)
if(!dir.exists(here('Output', SiteID_FileName, 'Historical'))) {dir.create(here('Output', SiteID_FileName, 'Historical'))}; outLocationPathHist <- here('Output', SiteID_FileName, 'Historical')



#######################################################################
### Retrieve and format streamflow data ###

# Get streamflow data
gage_data <- get_gage_data(GageSiteID, FALSE, FALSE, dataPath)
DailyStream <- gage_data$DailyStream
meas_flow_daily <- as.data.frame(gage_data$meas_flow_daily); meas_flow_mon <- gage_data$meas_flow_mon
meas_flow_daily$date <- as.Date(rownames(meas_flow_daily)); rownames(meas_flow_daily) <- NULL; colnames(meas_flow_daily) <- c('MeasMM', 'date')

# Get watershed data
INFO <- readNWISInfo(GageSiteID, "00060", interactive = FALSE)
watershed_area <- INFO$drain_area_va  # in square miles

# information for EGRET package
#printqUnitCheatSheet()
#printFluxUnitCheatSheet()

# Remove incomplete water years of data 
# Start and end dates correspond to beginning/end of water year
startDate_meas <- DailyStream$date[1]
if((month(startDate_meas) >= 10) & (day(startDate_meas) > 1)){
  startDate_meas <- as.Date(paste(year(startDate_meas)+1, 10, 1, sep='-'))
} else{
  startDate_meas <- as.Date(paste(year(startDate_meas), 10, 1, sep='-'))
}

endDate_meas <- DailyStream$date[nrow(DailyStream)]
if(month(endDate_meas) >= 10){
  endDate_meas <- as.Date(paste(year(endDate_meas), 9, 30, sep='-'))
} else{
  endDate_meas <- as.Date(paste(year(endDate_meas)-1, 9, 30, sep='-'))
}

# Add date information to daily dataframe
meas_flow_daily$water_year <- sapply(meas_flow_daily$date, get_water_year)
meas_flow_daily$year <- year(meas_flow_daily$date); meas_flow_daily$day <- yday(meas_flow_daily$date)
meas_flow_daily <- meas_flow_daily[meas_flow_daily$day != 366, ]

# Add water day
meas_flow_daily <- meas_flow_daily %>% group_by(water_year) %>%
  mutate(water_day = (as.integer(difftime(date,ymd(paste0(water_year - 1 ,'-09-30')), units = "days"))))

# Add date information to monthly dataframe
meas_flow_mon$Month <- month(as.Date(paste(meas_flow_mon$YrMon, "-01", sep=""), format="%Y-%m-%d"))
meas_flow_mon$Year <- year(as.Date(paste(meas_flow_mon$YrMon, "-01", sep=""), format="%Y-%m-%d"))

# Select only complete water years for all dataframes
# this can be commented out if all data wants to be considered; limited to full water years for statistical analyses
meas_flow_daily <- meas_flow_daily[meas_flow_daily$date >= startDate_meas & meas_flow_daily$date <= endDate_meas, ]
DailyStream <- DailyStream[DailyStream$date >= startDate_meas & DailyStream$date <= endDate_meas, ]
meas_flow_mon <- meas_flow_mon[as.Date(paste(meas_flow_mon$YrMon, "-01", sep=""), format="%Y-%m-%d") >= startDate_meas & as.Date(paste(meas_flow_mon$YrMon, "-01", sep=""), format="%Y-%m-%d") <= endDate_meas, ]


#######################################################################
### EGRET plots ###

# Aggregate in EGRET package format for plotting
eList <- as.egret(INFO, DailyStream)

# Adjust window in EGRET formatting based on available years of data
half_window <- (year(endDate_meas) - year(startDate_meas)) / 2
eList <- setPA(eList, window = half_window)

# plot summary statistics 
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Summary_Statistics.jpg"), width=1000, height=600); par(mfrow = c(2, 4), oma = c(0, 0, 4, 0))
  for(i in 1:8){
    plotFlowSingle(eList, istat=i, qUnit=1, printStaName=FALSE, lwd=2)
  }
  title(main=paste("Summary Statistics for", SiteID), outer=TRUE, cex.main = 1.5)
  dev.off()
}

# changes in variability - adjust window
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Variability_Change.jpg"), width=800, height=300); par(mfrow = c(1, 3), oma = c(0, 0, 4, 0))
  plotSDLogQ(eList, window=half_window, printStaName=FALSE)
  eList <- setPA(eList, paStart = 6, paLong = 3); plotSDLogQ(eList, window=4, printStaName=FALSE)
  eList <- setPA(eList, paStart = 12, paLong = 3); plotSDLogQ(eList, window=4, printStaName=FALSE)
  title(main=paste("Variability for", SiteID), outer=TRUE, cex.main = 1.5)
  dev.off()
}

# reset the window to the entire water year
eList <- setPA(eList, paStart=10, paLong=12)

# days above certain thresholds
plotQTimeDaily(eList, qLower=100)


#######################################################################
### Historical measured streamflow summary plots ###

### Heatmap - monthly 
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Historical_Monthly_Heatmap.jpg"), width=600, height=400)
  plot <- ggplot(meas_flow_mon, aes(factor(Month), Year, fill = MeasMM)) + geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
    labs(title = "Measured monthly streamflow", x='Month', fill = "Streamflow [mm]") +
    scale_x_discrete(labels = c("1" = "January", "2" = "February", "3" = "March", 
                                "4" = "April", "5" = "May", "6" = "June", 
                                "7" = "July", "8" = "August", "9" = "September", 
                                "10" = "October", "11" = "November", "12" = "December")) +
    nps_theme() + theme(axis.text.x = element_text(angle = 90))
  print(plot); dev.off()
}


### Heatmap - daily
meas_flow_daily$day <- factor(meas_flow_daily$day, levels = unique(meas_flow_daily$day))
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Historical_Daily_Heatmap.jpg"), width=600, height=400)
  plot <- ggplot(meas_flow_daily, aes(day, year, fill = as.numeric(MeasMM))) + geom_tile() +
    scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), trans='log', breaks=c(min(meas_flow_daily$MeasMM), 0.01, 0.1, 10, 100, max(meas_flow_daily$MeasMM)), 
                         labels=c(sprintf('%.3f',  min(meas_flow_daily$MeasMM)),'.01','0.1','10', '100', sprintf('%.0f',  max(meas_flow_daily$MeasMM)))) +
    labs(title = "Measured daily streamflow", x='Month', fill = "Streamflow [mm]") +
    scale_x_discrete(breaks=c("1", "32", "60", "91", "121", "152", "182", "213", "244", "274", "305", "335"), 
                     labels = c("1" = "January", "32" = "February", "60" = "March", "91" = "April", "121" = "May", "152" = "June", 
                                "182" = "July", "213" = "August", "244" = "September", "274" = "October", "305" = "November", "335" = "December")) +
    theme(axis.text.x = element_text(angle = 90))
  print(plot); dev.off()
}


### 3D raster hydrograph
# must save using viewer panel in R Studio
meas_flow_daily$day <- as.numeric(meas_flow_daily$day)
z_matrix <- reshape((meas_flow_daily %>% select(c('day','year','MeasMM')) %>% arrange(as.numeric(day))), idvar = "year", timevar = "day", direction = "wide")
rownames(z_matrix) <- z_matrix$year; colnames(z_matrix) <- sub("MeasMM.", "", colnames(z_matrix))
z_matrix <- as.matrix(z_matrix %>% select(-'year'))
fig <- plot_ly() %>% layout(title = paste(SiteID, "3D Hydrograph"),
                            scene = list(yaxis = list(title = 'Year'), xaxis = list(title = 'Day of year'), zaxis = list(title = 'Streamflow [mm]')), 
                            coloraxis = list(colorbar = list(title = 'Streamflow [mm]')))
for(year in unique(meas_flow_daily$year)) {
  year_data <- meas_flow_daily[meas_flow_daily$year == year, ]
  fig <- fig %>% add_trace(x = year_data$day, y = rep(year, nrow(year_data)), z = year_data$MeasMM, color = year_data$MeasMM, #colors = c('YlGnBu'), 
                           type = 'scatter3d', mode = 'lines', line = list(width = 4), 
                           name = as.character(year), showlegend=FALSE)
  # simulate fill - not sure how to do this
  # fig <- fig %>% add_trace(x = year_data$day, y = rep(year, nrow(year_data)), z = rep(min(year_data$Meas), nrow(year_data)), 
  #   fill = 'tozeroy', color = year_data$Meas, colors = 'YlGnBu', type = 'scatter3d', mode = 'lines',
  #   fill= 'tonexty', fillcolor = 'rgba(0, 100, 255, 0.3)', line = list(width = 0), showlegend = FALSE)
}
fig



#######################################################################
### Historical measured streamflow trend plots for specific metrics ###

# General Mann-Kendall test on daily streamflow ? should probably be seasonal
meas_mk <- SeasonalMannKendall(ts(meas_flow_daily$MeasMM, start=c(year(startDate_meas), 1), frequency=365))
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Daily_Flow_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(meas_flow_daily, aes(x = date, y = MeasMM)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    #geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Streamflow (mm)", title = "Daily Measured Streamflow", color='') +
    nps_theme() + theme(legend.position = 'bottom') + scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    #annotate("text", x = max(meas_flow_daily$date), y = max(meas_flow_daily$MeasMM), label = label, color = "black", hjust = 1, vjust = 1) + 
    scale_y_log10()
  print(plot_meas); dev.off()
}

# Seasonal Mann-Kendall test on monthly streamflow
meas_mk <- SeasonalMannKendall(ts(meas_flow_mon$MeasMM, start=c(year(startDate_meas), 1), frequency=12))
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Monthly_Flow_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(meas_flow_mon, aes(x = as.yearmon(YrMon), y = MeasMM)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Monthly Streamflow (mm)", title = "Monthly Measured Streamflow", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(as.yearmon(meas_flow_mon$YrMon)), y = max(meas_flow_mon$MeasMM), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}



# Annual streamflow volume
meas_flow_ann <- meas_flow_daily %>% group_by(water_year) %>% summarize(MeasMM = sum(MeasMM))
meas_mk <- MannKendall(meas_flow_ann$MeasMM)
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Annual_Volume_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(meas_flow_ann, aes(x = water_year, y = MeasMM)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Annual Streamflow (mm)", title = "Annual Measured Streamflow", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(meas_flow_ann$water_year), y = max(meas_flow_ann$MeasMM), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}


# High flows (above 95%)
high_flow_q = 0.95
high_flow_mm = quantile(meas_flow_daily$MeasMM, 0.95, na.rm=TRUE)

hist_high <- as.data.frame(meas_flow_daily %>% mutate(high_flow = ifelse(MeasMM >= high_flow_mm, 1, 0)) %>%
                                 group_by(water_year) %>% dplyr::summarize(days = sum(high_flow)))
meas_mk <- MannKendall(hist_high$days)
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Annual_High_Flow_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(hist_high, aes(x = water_year, y = days)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Number of days per year", title = "Days Above Historical 95th Percentile", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(hist_high$water_year), y = max(hist_high$days, na.rm=TRUE), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}


# Low flows (below 5%)
low_flow_q = 0.05
low_flow_mm = quantile(meas_flow_daily$MeasMM, 0.05, na.rm=TRUE)

hist_low <- as.data.frame(meas_flow_daily %>% mutate(low_flow = ifelse(MeasMM <= low_flow_mm, 1, 0)) %>%
                             group_by(water_year) %>% dplyr::summarize(days = sum(low_flow)))
meas_mk <- MannKendall(hist_low$days)
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Annual_Low_Flow_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(hist_low, aes(x = water_year, y = days)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Number of days per year", title = "Days Below Historical 5th Percentile", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(hist_low$water_year), y = max(hist_low$days, na.rm=TRUE), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}


# 50% flow date
hist_ct <- as.data.frame(meas_flow_daily %>% mutate(tq = MeasMM*water_day) %>%
                          group_by(water_year) %>% dplyr::summarize(ct = sum(tq)/sum(MeasMM)))
meas_mk <- MannKendall(hist_ct$ct)
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Annual_50th_Flow_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(hist_ct, aes(x = water_year, y = ct)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Days after Oct 1", title = "Historical 50% Flow Date", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(hist_ct$water_year), y = max(hist_ct$ct, na.rm=TRUE), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}


# Q7 min %>% filter(!is.na(waterYear))
hist_q7 <- DailyStream %>% group_by(waterYear) %>% dplyr::summarize(min_q7 = ifelse(all(is.na(Q7)), NA, min(Q7, na.rm = TRUE)), 
                                                                    max_q7 = ifelse(all(is.na(Q7)), NA, max(Q7, na.rm = TRUE)), 
                                                                    avg_q7 = ifelse(all(is.na(Q7)), NA, mean(Q7, na.rm = TRUE)))
meas_mk <- MannKendall(hist_q7$min_q7)
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Annual_Q7Min_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(hist_q7, aes(x = waterYear, y = min_q7)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Streamflow [m3/s]", title = "Minimum 7 Day Flow (Q7 Min)", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(hist_q7$waterYear), y = max(hist_q7$min_q7, na.rm=TRUE), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}

# Q7 max
meas_mk <- MannKendall(hist_q7$max_q7)
if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
}else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Annual_Q7Max_Trends.jpg"), width=600, height=400)
  plot_meas <- ggplot(hist_q7, aes(x = waterYear, y = max_q7)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Streamflow [m3/s]", title = "Maximum 7 Day Flow (Q7 Max)", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(hist_q7$waterYear), y = max(hist_q7$max_q7, na.rm=TRUE), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}



#######################################################################
### Customizable plots: metric, season, etc ###

# define flow level in cfs 
flow_level <- 1250
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
if(make_plots){
  jpeg(file=paste0(outLocationPathHist, "/", "Days_", comparison, "_", flow_level, "_", month.abb[mos][1], "_", month.abb[mos][length(mos)],".jpg"), width=600, height=400)
  plot_meas <- ggplot(hist_threshold, aes(x = waterYear, y = days)) + geom_line(aes(color = 'Measured'), na.rm=TRUE, linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Water Year", y = "Days", title = paste('Days', comparison, flow_level, 'cfs, ', month.abb[mos][1], '-', month.abb[mos][length(mos)]), color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(hist_threshold$waterYear), y = max(hist_threshold$days, na.rm=TRUE), label = label, color = "black", hjust = 1, vjust = 1)
  print(plot_meas); dev.off()
}



