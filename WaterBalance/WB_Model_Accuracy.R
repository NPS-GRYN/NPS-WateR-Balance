# ---------------------------------------------------------------------
# This script contains code for evaluating the accuracy of the water balance model
# by comparing its estimations against observations for the historical time 
# period. 
# This script is not intended to be run independently; it must be called from 
# the main Run_WB_Model.R script.
#
# ---------------------------------------------------------------------

#######################################################################
### HISTORICAL AET #####

# Daily, monthly, annual aggregation with measured and modeled ET
hist_aet_daily <- merge(xts(with(DailyWB, cbind(AET)), order.by = as.Date(DailyWB$date)), DailyET); hist_aet_daily <- hist_aet_daily[complete.cases(hist_aet_daily),]
colnames(hist_aet_daily) <- c("Mod", "Meas")
hist_aet_monthly <- apply.monthly(xts(with(DailyWB, cbind(AET)), order.by = as.Date(DailyWB$date)), function(x) {colSums(x, na.rm = TRUE)})
index(hist_aet_monthly) <- as.Date(format(index(hist_aet_monthly), "%Y-%m-01"))
hist_aet_monthly <- merge(hist_aet_monthly, MonthlyET); hist_aet_monthly <- hist_aet_monthly[complete.cases(hist_aet_monthly),]; colnames(hist_aet_monthly) <- c("Mod", "Meas")
hist_aet_ann <- apply.yearly(hist_aet_monthly, function(x) {colSums(x, na.rm = TRUE)})


#######################################################################
### Create summary plots ###

# scatterplot of Historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
# there are two trend lines in the scatter plot because the intercept is set to 0 in one and allowed to vary in the other
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_AET_Scatter.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  meas <- coredata(hist_aet_daily$Meas); mod <- coredata(hist_aet_daily$Mod)
  plot(mod, meas, main =paste("Daily Average AET\nNSE:",round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # monthly
  mod <- coredata(hist_aet_monthly$Mod); meas <- coredata(hist_aet_monthly$Meas)
  plot(mod, meas, main =paste("Monthly Total AET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # annual
  mod <- coredata(hist_aet_ann$Mod); meas <- coredata(hist_aet_ann$Meas)
  plot(mod,meas, main =paste("Annual Total AET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled AET (mm)", ylab = "Measured AET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  dev.off()
}

# time series plot of historical Measured vs Modeled AET for daily, monthly, annual aggregation
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_AET_TimeSeries.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  plot(hist_aet_daily[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Daily AET (mm)", main = "Daily", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # monthly
  plot(hist_aet_monthly[,c("Mod", "Meas")], type = "l", lwd = 2, xlab = "Date", ylab = "Monthly Sum AET (mm)", main = "Monthly", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # annual
  plot(hist_aet_ann[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Annual Sum AET (mm)", main = "Annual", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}

#######################################################################
### HISTORICAL PET #####
# ETo from OpenET is for a grass reference crop, so isn't accurate in areas 
# with more complex vegetation cover. These results should not be considered
# unless the model calibration has been PET to ETo.

# Daily, monthly, annual aggregation with measured and modeled ET
hist_pet_daily <- merge(xts(with(DailyWB, cbind(PET)), order.by = as.Date(DailyWB$date)), DailyETo); hist_pet_daily <- hist_pet_daily[complete.cases(hist_pet_daily),]; colnames(hist_pet_daily) <- c("Mod", "Meas")
hist_pet_daily$Mod <- hist_pet_daily$Mod * k_c
hist_pet_monthly <- apply.monthly(xts(with(DailyWB, cbind(PET)), order.by = as.Date(DailyWB$date)), function(x) {colSums(x, na.rm = TRUE)})
index(hist_pet_monthly) <- as.Date(format(index(hist_pet_monthly), "%Y-%m-01"))
hist_pet_monthly <- merge(hist_pet_monthly, MonthlyETo); hist_pet_monthly <- hist_pet_monthly[complete.cases(hist_pet_monthly),]; colnames(hist_pet_monthly) <- c("Mod", "Meas")
hist_pet_monthly$Mod <- hist_pet_monthly$Mod * k_c
hist_pet_ann <- apply.yearly(hist_pet_monthly, function(x) {colSums(x, na.rm = TRUE)})


#######################################################################
### Create summary plots ###

# scatterplot of Historical Measured vs Modeled PET for daily, monthly, annual aggregation
# there are two trend lines in the scatter plot because the intercept is set to 0 in one and allowed to vary in the other
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_PET_Scatter.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  meas <- coredata(hist_pet_daily$Meas); mod <- coredata(hist_pet_daily$Mod)
  plot(mod, meas, main =paste("Daily Average PET\nNSE:",round(NSE(mod, meas),digits=2)), xlab = "Modeled PET (mm)", ylab = "Measured PET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # monthly
  mod <- coredata(hist_pet_monthly$Mod); meas <- coredata(hist_pet_monthly$Meas)
  plot(mod, meas, main =paste("Monthly Total PET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled PET (mm)", ylab = "Measured PET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # annual
  mod <- coredata(hist_pet_ann$Mod); meas <- coredata(hist_pet_ann$Meas)
  plot(mod,meas, main =paste("Annual Total PET\nNSE:", round(NSE(mod, meas),digits=2)), xlab = "Modeled PET (mm)", ylab = "Measured PET (mm)")
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  dev.off()
}

# time series plot of historical Measured vs Modeled PET for daily, monthly, annual aggregation
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_PET_TimeSeries.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  plot(hist_pet_daily[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Daily PET (mm)", main = "Daily", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # monthly
  plot(hist_pet_monthly[,c("Mod", "Meas")], type = "l", lwd = 2, xlab = "Date", ylab = "Monthly Sum PET (mm)", main = "Monthly", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # annual
  plot(hist_pet_ann[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Annual Sum PET (mm)", main = "Annual", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}

