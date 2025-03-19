# ---------------------------------------------------------------------
# This script includes code to analyze historical measurements of streamflow. In addition to
# assessing historical trends, this script compares historical measurements to modeled
# streamflow, estimated using the IHACRES  rainfall-streamflow methodology.
# This code also provides preliminary analyses and visualizations of historical streamflow. 
# 
# EDITS IN PROGRESS
# moving code from main script into here
# make script so it can be run independently or in conjunction with main script
# figure out how to export 3D plot
# ---------------------------------------------------------------------


#######################################################################
# redefine gage site ID so this script can be run independently
GageSiteID <- GageSiteID



#######################################################################
### Get streamflow data in EGRET package format ###
gage_data <- get_gage_data(GageSiteID, FALSE, dataPath)
DailyStream <- gage_data$DailyStream
meas_flow_daily <- as.data.frame(gage_data$meas_flow_daily); meas_flow_mon <- gage_data$meas_flow_mon
meas_flow_daily$Date <- rownames(meas_flow_daily); rownames(meas_flow_daily) <- NULL; colnames(meas_flow_daily) <- c('MeasMM', 'Date')

INFO <- readNWISInfo(GageSiteID, "00060", interactive = FALSE)
eList <- as.egret(INFO, DailyStream)

# information for EGRET package
#printqUnitCheatSheet()
#printFluxUnitCheatSheet()

startY_meas <- year(DailyStream$Date[1]); endY_meas <-year(DailyStream$Date[nrow(DailyStream)])
half_window <- min(endY, endY_meas) - max(startY_meas, startY)
eList <- setPA(eList, window = half_window)

# plot summary statistics 
jpeg(file=paste0(outLocationPath, "/", "Summary_Statistics.jpg"), width=1000, height=600); par(mfrow = c(2, 4), oma = c(0, 0, 4, 0))
for(i in 1:8){
  plotFlowSingle(eList, istat=i, qUnit=1, printStaName=FALSE)
}
title(main=paste("Summary Statistics for", SiteID), outer=TRUE, cex.main = 1.5)
dev.off()

# changes in variability - adjust window
jpeg(file=paste0(outLocationPath, "/", "Variability_Change.jpg"), width=800, height=300); par(mfrow = c(1, 3), oma = c(0, 0, 4, 0))
plotSDLogQ(eList, window=4, printStaName=FALSE)
eList <- setPA(eList, paStart = 6, paLong = 3); plotSDLogQ(eList, window=4, printStaName=FALSE)
eList <- setPA(eList, paStart = 12, paLong = 3); plotSDLogQ(eList, window=4, printStaName=FALSE)
title(main=paste("Variability for", SiteID), outer=TRUE, cex.main = 1.5)
dev.off()

# reset the window to the entire water year
eList <- setPA(eList, paStart=10, paLong=12)

# days above certain thresholds
plotQTimeDaily(eList, qLower=100)


#######################################################################
### Plots of historical measured streamflow - general ###

### Heatmap - monthly 
meas_flow_mon$Month <- month(as.Date(paste(meas_flow_mon$YrMon, "-01", sep=""), format="%Y-%m-%d"))
meas_flow_mon$Year <- year(as.Date(paste(meas_flow_mon$YrMon, "-01", sep=""), format="%Y-%m-%d"))
jpeg(file=paste0(outLocationPath, "/", "Historical_Monthly_Heatmap.jpg"), width=600, height=400)
ggplot(meas_flow_mon, aes(factor(Month), Year, fill = MeasMM)) + geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu")) +
  labs(title = "Measured monthly streamflow", x='Month', fill = "Streamflow [mm]") +
  scale_x_discrete(labels = c("1" = "January", "2" = "February", "3" = "March", 
                              "4" = "April", "5" = "May", "6" = "June", 
                              "7" = "July", "8" = "August", "9" = "September", 
                              "10" = "October", "11" = "November", "12" = "December")) +
  nps_theme() + theme(axis.text.x = element_text(angle = 90))
dev.off()


### Heatmap - daily
meas_flow_daily$year <- year(meas_flow_daily$Date); meas_flow_daily$day <- yday(meas_flow_daily$Date)
meas_flow_daily <- meas_flow_daily[meas_flow_daily$day != 366, ]
meas_flow_daily$day <- factor(meas_flow_daily$day, levels = unique(meas_flow_daily$day))
jpeg(file=paste0(outLocationPath, "/", "Historical_Daily_Heatmap.jpg"), width=600, height=400)
ggplot(meas_flow_daily, aes(day, year, fill = as.numeric(MeasMM))) + geom_tile() +
  scale_fill_gradientn(colors = brewer.pal(9, "YlGnBu"), trans='log', breaks=c(min(hist_flow_daily_df$Meas), 0.01, 0.1, 10, 100, max(hist_flow_daily_df$Meas)), 
                       labels=c(sprintf('%.3f',  min(hist_flow_daily_df$Meas)),'.01','0.1','10', '100', sprintf('%.0f',  max(hist_flow_daily_df$Meas)))) +
  labs(title = "Measured daily streamflow", x='Month', fill = "Streamflow [mm]") +
  scale_x_discrete(breaks=c("1", "32", "60", "91", "121", "152", "182", "213", "244", "274", "305", "335"), 
                   labels = c("1" = "January", "32" = "February", "60" = "March", "91" = "April", "121" = "May", "152" = "June", 
                              "182" = "July", "213" = "August", "244" = "September", "274" = "October", "305" = "November", "335" = "December")) +
  theme(axis.text.x = element_text(angle = 90))
dev.off()


### 3D raster hydrograph
meas_flow_daily$day <- as.numeric(meas_flow_daily$day)
z_matrix <- reshape((meas_flow_daily %>% select(c('day','year','MeasMM')) %>% arrange(as.numeric(day))), idvar = "year", timevar = "day", direction = "wide")
rownames(z_matrix) <- z_matrix$year; colnames(z_matrix) <- sub("MeasMM.", "", colnames(z_matrix))
z_matrix <- as.matrix(z_matrix %>% select(-'year'))
fig <- plot_ly() %>% layout(title = paste(SiteID, "3D Hydrograph"),
                            scene = list(yaxis = list(title = 'Year'), xaxis = list(title = 'Day of year'), zaxis = list(title = 'Streamflow [mm]')), 
                            coloraxis = list(colorbar = list(title = 'Streamflow [mm]')))
for(year in unique(hist_flow_daily_df$year)) {
  year_data <- hist_flow_daily_df[hist_flow_daily_df$year == year, ]
  fig <- fig %>% add_trace(x = year_data$day, y = rep(year, nrow(year_data)), z = year_data$Meas, color = year_data$Meas, #colors = c('YlGnBu'), 
                           type = 'scatter3d', mode = 'lines', line = list(width = 4), 
                           name = as.character(year), showlegend=FALSE)
  # simulate fill - still not sure how to do this
  # fig <- fig %>% add_trace(x = year_data$day, y = rep(year, nrow(year_data)), z = rep(min(year_data$Meas), nrow(year_data)), 
  #   fill = 'tozeroy', color = year_data$Meas, colors = 'YlGnBu', type = 'scatter3d', mode = 'lines',
  #   fill= 'tonexty', fillcolor = 'rgba(0, 100, 255, 0.3)', line = list(width = 0), showlegend = FALSE)
}
fig


#######################################################################
### Plots of historical measured streamflow - analysis of metrics ###

meas_flow_ann <- meas_flow_daily %>% group_by(year) %>% summarize(MeasMM = sum(MeasMM))
