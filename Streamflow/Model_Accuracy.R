# ---------------------------------------------------------------------
# This script contains code for evaluating the accuracy of the streamflow model
# by comparing its estimations against observations for the historical time 
# period. 
# This script is not intended to be run independently; it must be called from 
# the main Run_Streamflow_Model.R script.
#
# EDITS IN PROGRESS:
# make sure it all works
# ---------------------------------------------------------------------


#######################################################################
### Create dataframes for aggregations ###

# Daily aggregation with measured and modeled streamflow
hist_flow_daily <- merge(xts(with(DailyDrain, cbind(adj_runoff, Quick, Slow, Very_Slow, total)), order.by = as.Date(DailyDrain$date)), meas_flow_daily)
hist_flow_daily <- hist_flow_daily[complete.cases(hist_flow_daily),]
colnames(hist_flow_daily) <- c("adj_runoff", "quick", "slow", "very slow", "Mod", "Meas")

# Monthly aggregation with measured and modeled streamflow
MeasMod <- MeasModWB(DailyDrain = DailyDrain, meas_flow_mon = meas_flow_mon, cutoffYear = cutoffYear)

# Annual aggregation with measured and modeled streamflow
hist_flow_ann <- as.data.frame(matrix(NA, nrow = nrow(apply.yearly(hist_flow_daily[,"adj_runoff"], sum)),
                                      ncol = ncol(hist_flow_daily), dimnames = list(c(), colnames(hist_flow_daily))))
for(i in 1:ncol(hist_flow_daily)){
  hist_flow_ann[,i] <- apply.yearly(hist_flow_daily[,i], sum)
}


#######################################################################
### Create summary plots ### - DOUBLE CHECK ALL OF THESE FOR ACCURACY/CONSISTENTCY

# scatterplot of Historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
# there are two trend lines in the scatter plot because the intercept is set to 0 in one and allowed to vary in the other
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_Scatter.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  meas <- coredata(hist_flow_daily$Meas); mod <- coredata(hist_flow_daily$Mod)
  plot(mod, meas, main ="Daily Average Streamflow", xlab = "Modeled Streamflow (mm)", ylab = "Measured Streamflow (mm)")
  text(0.25*par("usr")[2], 0.85*par("usr")[4], paste('NSE:',round(NSE(mod, meas),digits=2)), cex = 3)
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # monthly
  mod <- MeasMod$Mod; meas <- MeasMod$Meas
  plot(mod, meas, main ="Monthly Total Streamflow", xlab = "Modeled Streamflow (mm)", ylab = "Measured Streamflow (mm)")
  text(0.35*par("usr")[2], 0.85*par("usr")[4], paste('NSE:',round(NSE(mod, meas),digits=2)), cex = 3)
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  
  # annual
  mod <- coredata(hist_flow_ann$Mod); meas <- coredata(hist_flow_ann$Meas)
  plot(mod,meas, main ="Annual Total Streamflow", xlab = "Modeled Streamflow (mm)", ylab = "Measured Streamflow (mm)")
  text(0.6*par("usr")[2], 0.85*par("usr")[4], paste('NSE:',round(NSE(mod, meas),digits=2)), cex = 3)
  abline(lm(meas ~ 0 + mod), col= "red")
  abline(lm(meas ~ mod), col= "red")
  dev.off()
}

# time series plot of historical Measured vs Modeled Streamflow for daily, monthly, annual aggregation
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_TimeSeries.jpg"), width=1000, height=400); par(mfrow=c(1,3))
  
  # daily
  plot(hist_flow_daily[,c('Mod','Meas')], type = "l", lwd = 2, xlab = "Date", ylab = "Daily Streamflow (mm)", main = "Daily", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # monthly
  plot(xts(MeasMod[,c("Mod", "Meas")], order.by = ym(MeasMod$YrMon)), 
       type = "l", lwd = 2, xlab = "Date", ylab = "Monthly Sum Streamflow (mm)", main = "Monthly", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  
  # annual
  plot(xts(hist_flow_ann[,c('Mod','Meas')], order.by=as.Date(index(apply.yearly(hist_flow_daily, sum)))), 
       type = "l", lwd = 2, xlab = "Date", ylab = "Annual Sum Streamflow (mm)", main = "Annual", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}

# time series plot of historical annual streamflow trends (measured and modeled)
# trend analysis assumes p value of < 0.05 is significant
if(make_plots){
  meas_mk <- MannKendall(hist_flow_ann$Meas)
  if(meas_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', meas_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', meas_mk$sl)}
  
  plot_meas <- ggplot(hist_flow_ann, aes(x = index(Meas), y = Meas)) + geom_line(aes(color = 'Measured'), linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Date", y = "Annual Streamflow (mm)", title = "Annual Measured Streamflow", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Measured" = "black", "Trend" = "red")) +
    annotate("text", x = max(index(hist_flow_ann$Meas)), y = max(hist_flow_ann$Meas), label = label, color = "black", hjust = 1, vjust = 1)
  plot_meas
  
  mod_mk <- MannKendall(hist_flow_ann$Mod)
  if(mod_mk$sl <= 0.05){label <- sprintf('Trend: Significant \n p-value: %.2f', mod_mk$sl)
  }else{label <- sprintf('Trend: Not significant \n p-value: %.2f', mod_mk$sl)}
  
  plot_mod <- ggplot(hist_flow_ann, aes(x = index(Mod), y = Mod)) + geom_line(aes(color = 'Modeled'), linewidth=1) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, aes(color = 'Trend')) +
    labs(x = "Date", y = "Annual Streamflow (mm)", title = "Annual Modeled Streamflow", color='') +
    nps_theme() + theme(legend.position = 'bottom') +
    scale_color_manual(values = c("Modeled" = "black", "Trend" = "red")) +
    annotate("text", x = max(index(hist_flow_ann$Mod)), y = max(hist_flow_ann$Mod), label = label, color = "black", hjust = 1, vjust = 1)
  plot_mod
  
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_Trends.jpg"), width=1000, height=400)
  grid.arrange(plot_meas, plot_mod, ncol = 2) 
  dev.off()
}

# log plot of historical measured vs modeled daily streamflow 
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_Measured_Modeled_Daily_Log.jpg"), width=700, height=400)
  plot(hist_flow_daily[,c('Mod','Meas')], type = "l", log=TRUE, lwd = 2, xlab = "Date", ylab = "Streamflow (mm)", main = "Log Plot of Daily Historical Streamflow", col=c('red','black'))
  print(xts::addLegend("topleft", legend.names = c("Modeled", "Measured"), lty=1, col= c("red", "black")))
  dev.off()
}



#######################################################################
### Assess model accuracy ###

### Model accuracy: high and low quantiles ###
high_flow_q = 0.99
high_flow_meas = quantile(hist_flow_daily$Meas, high_flow_q)
high_flow_mod = quantile(hist_flow_daily$Mod, high_flow_q) 

low_flow_q = 0.1
low_flow_meas = quantile(hist_flow_daily$Meas, low_flow_q)
low_flow_mod = quantile(hist_flow_daily$Mod, low_flow_q)


### Simple quantile regression ###
historic_75 = quantile(hist_flow_daily$Meas, 0.75); historic_25 = quantile(hist_flow_daily$Meas, 0.25)
high_flow = hist_flow_daily[hist_flow_daily$Meas > high_flow_meas, ]

if(make_plots){
  #jpeg(file=paste0(figPath, "/", paste0(gsub(" ", "_", title), "_Measured_Modeled_Scatter_75th_Percentile.jpg")))
  plot(coredata(high_flow$Mod), coredata(high_flow$Meas), main='High Flow (75th Percentile) Daily Historical Streamflow', xlab = "Modeled Streamflow", ylab = "Measured Streamflow", 
       xlim=c(pmin(min(high_flow$Mod), min(high_flow$Meas)), pmax(max(high_flow$Mod), max(high_flow$Meas))), 
       ylim=c(pmin(min(high_flow$Meas), min(high_flow$Mod)), pmax(max(high_flow$Meas), max(high_flow$Mod))))
  #abline(lm(coredata(high_flow$Meas) ~ 0 + coredata(high_flow$Mod)), col= "red")
  abline(lm(coredata(high_flow$Meas) ~ coredata(high_flow$Mod)), col= "red")
  
  # calculate statistics to locate on the plot
  nse_plot = NSE(coredata(high_flow$Mod), coredata(high_flow$Meas))
  r2_plot = R2(coredata(high_flow$Mod), coredata(high_flow$Meas))
  
  y_txt <- max(coredata(high_flow$Meas)) * (1/4)
  x_txt <- max(coredata(high_flow$Mod)) * (8/9)
  
  text(x = x_txt, y = y_txt, labels = sprintf("NSE: %.2f \nR2: %.2f",nse_plot, r2_plot), col = "red", cex = 1.2) 
  #dev.off()
}


### Quantile regression plot ###
model <- rq(Meas ~ Mod, data = hist_flow_daily, tau = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99))
r2_values <- sapply(model$tau, function(t) calculate_pseudo_r2(model, hist_flow_daily, t)) # calculate psuedo r2 for each quantile
r2_results <- data.frame(tau = model$tau, pseudo_r2 = r2_values)

colors = brewer.pal(n = 7, name = "RdYlBu") 
line_data <- data.frame(intercept = coef(model)[seq(1, length(coef(model)), by = 2)], slope = coef(model)[seq(2, length(coef(model)), by = 2)],
                        color = colors, tau = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99))
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "Historical_QuantileRegression_Scatter.jpg"), width=800, height=500)
  plot <- ggplot(hist_flow_daily, aes(Mod,Meas)) + geom_point() +   geom_abline(data = line_data, aes(intercept = intercept, slope = slope, color = factor(tau)), linewidth=1.1) +
    scale_color_manual(values = colors, labels = c(0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 0.99), name = "Quantiles") +
    theme(legend.position = "right")+ scale_x_continuous(limits = c(1, NA)) + labs(title="Quantile Regression of Daily Historical Streamflow",y="Measured Streamflow", x="Modeled Streamflow")+
    scale_color_manual(values = colors, labels = c(bquote(bold("0.01:") * " -0.047"), bquote(bold("0.1:")*" 0.132"), bquote(bold("0.25:")*" 0.395"), bquote(bold("0.5:")*" 0.624"), bquote(bold("0.75:")*" 0.602"), bquote(bold("0.9:")*" 0.323"), bquote(bold("0.99:")*" -2.09")), name = "Pseudo R2 by Quantile") + 
    nps_theme() 
  plot
  dev.off()
}
