# ---------------------------------------------------------------------
# This script contains code for optimizing the streamflow model. This occurs in
# two steps: first the water balance model is optimized, then the IHACRES 
# streamflow model is optimized. Both optimizations are included here. 
# This script is not intended to be run independently; it must be called from 
# the main Run_Streamflow_Model.R script.
#
# EDITS IN PROGRESS:
# make the workflow make sense!!
#look into optim(): increment jrange in whole numbers ?
# ---------------------------------------------------------------------


# Create dataframe to store optimization results
results <- data.frame(SiteID = SiteID, start = startDate, end = endDate, PETMethod = PETMethod, optimization = optimization,
                      GridMET = GridMET, lon = lon, lat = lat,
                      startY = startY, startM = startM, startD = startD, endY = endY, endM = endM, endD = endD,
                      cutoffYear = cutoffYear, NonZeroDrainInitCoeff = NonZeroDrainInitCoeff, incompleteMonths = incompleteMonths)


#######################################################################
#######################################################################
### First optimization ###
# Optimize water balance variables according to the NSE of monthly summed streamflow over historical period

parms<- c(gw_add = gw_add, vfm = vfm, jrange = jrange, hock =  hock, hockros = hockros,dro = dro, mondro = mondro,
          aspect = aspect,slope= slope, shade.coeff= shade.coeff, SWC.Max = SWC.Max, jtemp = jtemp)

#run the optimization routine
strtTimeM <-Sys.time()
set.seed(123) #this ensures reproducibility each time
WBcoeffs <- tibble()

# Use genetic algorithm (GA) for optimization
optMonth_init <- ga(type = "real-valued", fitness = function(x) 
  WB_Optim(c(gw_add=x[1], vfm=x[2], jrange=x[3], hock=x[4], hockros=x[5], dro=x[6], mondro=x[7], aspect=x[8], slope=x[9], shade.coeff=x[10], SWC.Max=x[11], jtemp=x[12]), 
           meas_flow_daily = meas_flow_daily, cutoffYear = cutoffYear, q0=q0, s0=s0, v0=v0,qa=qa, qb=qb, sa=sa, sb=sb,va=va, vb=vb, Soil.Init = Soil.Init, 
           Snowpack.Init = Snowpack.Init, T.Base = T.Base, PETMethod= PETMethod, DailyClimData = DailyClimData, lat=lat,lon=lon, meas_flow_mon = meas_flow_mon), 
  lower=WB_lower, upper=WB_upper)
elpTimeM <- Sys.time() - strtTimeM

# Define the water balance variables from the best run
optValuesM <- data.frame(nseM = optMonth_init@fitnessValue, optMonth_init@solution)
gw_add=optValuesM$gw_add; vfm=optValuesM$vfm; jrange=optValuesM$jrange; hock=optValuesM$hock
hockros=optValuesM$hockros; dro=optValuesM$dro; mondro=optValuesM$mondro; aspect=optValuesM$aspect
slope=optValuesM$slope; shade.coeff=optValuesM$shade.coeff; SWC.Max=optValuesM$SWC.Max; jtemp=optValuesM$jtemp

# store and save results
results = data.frame(results, optValuesM, elpTimeM = elpTimeM)
saveRDS(WBcoeffs, file = paste0(outLocationPath, "/WBcoeffs.rds"))


### PLOTS
if(make_plots){
  if (dev.cur() != 1) dev.off()
  # Parallel coordinates plot of WB coeffs
  jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_ParallelCoords.jpg"), width=1600, height=800, res=100)
  print(ggparcoord(data=WBcoeffs, columns=1:12, groupColumn=13, scale="uniminmax") + 
          scale_color_gradient(low = "black", high = "gray90") + nps_theme()) # for red: "darkred" and "#fee5d9"
  dev.off()
  
  # Scatterplots of WB coeffs
  WBcoeffs_long <- reshape(WBcoeffs, varying = names(WBcoeffs)[1:12], v.names = "WBcoeffs", timevar = "Variable", 
                           times = names(WBcoeffs)[1:12], direction = "long")
  jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_Scatter.jpg"), width=800, height=500)
  print(ggplot(WBcoeffs_long, aes(x = WBcoeffs, y = nseM)) + geom_point() +
          facet_wrap(~ Variable, scales = 'free') + nps_theme() +
          labs(title = 'Water Balance Coefficients', x='', y = 'Monthly NSE'))
  dev.off()
  
  # Range plot of optimal WB coeffs
  wb_optim <- data.frame(var=c('Groundwater Addition', 'Volume Forcing Multiplier', 'Jennings Temperature Range','Hock','Hock Rain on Snow','Direct Runoff','Mondro','Aspect','Slope','Shade Coefficient','Max Soil Water Content','Jennings Temperature'),
                         value=c(gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope, shade.coeff, SWC.Max, jtemp),
                         lower=WB_lower, upper=WB_upper)
  jpeg(file=paste0(outLocationPath, "/", "WB_Optim_Coefficients.jpg"), width=600, height=400); par(mfrow = c(3, 4))
  for (i in 1:12){
    len <- ceiling(wb_optim$upper[i])-floor(wb_optim$lower[i])
    plot(c(wb_optim$lower[i], wb_optim$upper[i]), c(0, 0), type = "n", xlab = "", ylab = "",
         main = wb_optim$var[i], xlim = range(c(wb_optim$lower[i]-(len/5), wb_optim$upper[i]+(len/5))), ylim = c(-1, 1), xaxt = "n", yaxt = "n")
    segments(wb_optim$lower[i], 0, wb_optim$upper[i], 0, col = "black", lwd = 2)
    points(wb_optim$value[i], 0, col = "red", pch = 19, cex = 1.5)
    text(wb_optim$value[i], 0.3, sprintf("%.2f", wb_optim$value[i]), col='red')
    axis(1, at = seq(floor(wb_optim$lower[i]), ceiling(wb_optim$upper[i]), by=(len/5)))
  }
  dev.off()
}




#######################################################################
#######################################################################
### Second optimization ###
# Optimize IHACRES A and B coefficients according to the NSE of daily streamflow over historical period
# check how it works w 2 flow components - used to be bad

#Re run the Water Balance model as an input for the second optimization with the optimized water balance variables
DailyWB<- WB(DailyClimData, gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope,
             shade.coeff, jtemp,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod,lat, lon)

# Define variables for optimization
IHACRES_lower <- c(l_qa=0, l_qb=0, l_sa=0, l_sb=0, l_va=0)
IHACRES_upper <- c(u_qa=1, u_qb=1, u_sa=1, u_sb=1, u_va=1)
strtTimeD <-Sys.time()
set.seed(123)


### use GA (genetic algorithm) to explore parameter space and estimate IHACRES flow coefficients ###
IHcoeffs <- tibble()
optDaily_init <- ga(type = "real-valued", fitness = function(x) IHACRES_optim(c(qa=x[1], qb=x[2], sa=x[3], sb=x[4], va=x[5]), q0=q0, s0=s0, v0=v0, DailyWB=DailyWB, meas_flow_daily= meas_flow_daily, cutoffYear = cutoffYear), lower =IHACRES_lower, upper = IHACRES_upper)
if(flow_components==3){
  qa <- optDaily_init@solution[1]; qb <- optDaily_init@solution[2]; sa <- optDaily_init@solution[3]; sb <- optDaily_init@solution[4]; va <- optDaily_init@solution[5]; vb<-calc_vb(qa, qb, sa, sa, va)
} else if(flow_components==2){
  optSol <- IHcoeffs[which.max(IHcoeffs$nseD), ]
  qa<-optSol$qa; qb<-optSol$qb; sa<-optSol$sa; sb<-optSol$sb; va<-optSol$va; vb<-optSol$vb
}
saveRDS(IHcoeffs, file = paste0(outLocationPath, "/IHcoeffs_initial.rds"))

# plot 
if(make_plots){
  # Parallel coordinates plot of IH coeffs
  jpeg(file=paste0(outLocationPath, "/", "IHACRES_Coeffs_ParallelCoords.jpg"), width=600, height=400)
  print(ggparcoord(data=IHcoeffs, columns=1:6, groupColumn=7, scale="uniminmax") + 
          scale_color_gradient(low = "black", high = "gray90") + nps_theme()) # for red: "darkred" and "#fee5d9"
  dev.off()
  
  # Scatterplots of IH coeffs
  IHcoeffs_long <- reshape(IHcoeffs, varying = names(IHcoeffs)[1:6], v.names = "IHcoeffs", timevar = "Variable", 
                           times = names(IHcoeffs)[1:6], direction = "long")
  jpeg(file=paste0(outLocationPath, "/", "IHACRES_Coeffs_Scatter.jpg"), width=600, height=500)
  print(ggplot(IHcoeffs_long, aes(x = IHcoeffs, y = nseD)) + geom_point() +
          facet_wrap(~ Variable, scales = 'free') + nps_theme() +
          labs(title = 'IHACRES Coefficients', x='', y = 'Daily NSE'))
  dev.off()
}


### Run optim() to find optimal IHACRES flow coefficients ###
IHcoeffs <- tibble()
IHACRES_parms <- c(qa=qa, qb=qb, sa=sa, sb=sb, va=va)
optDaily <- optim(par = IHACRES_parms, fn = IHACRES_optim, method = "L-BFGS-B",
                  lower = IHACRES_lower, upper = IHACRES_upper, hessian=TRUE, control = list(fnscale = -1, factr = '1e-2'),
                  #these parameters are carried through to IHACRES_optim
                  q0=q0, s0=s0, v0=v0, DailyWB=DailyWB, meas_flow_daily=meas_flow_daily, cutoffYear=cutoffYear)

elpTimeD <- Sys.time() - strtTimeD

# redefine IHACRES variables from optimized run
# EDIT THIS - unnecessarily complicated
optValuesD <- data.frame(nseD = optDaily$value, t(as.matrix(optDaily$par)))#data.frame(t(data.frame(optDaily$par)))
optValuesD$vb <- calc_vb(optValuesD$qa, optValuesD$qb, optValuesD$sa, optValuesD$sb, optValuesD$va)
qa=optValuesD$qa; qb=optValuesD$qb; sa=optValuesD$sa; sb=optValuesD$sb; va=optValuesD$va; vb=optValuesD$vb

# store and save results
results<- data.frame(results, qa=optValuesD$qa, qb=optValuesD$qb, sa=optValuesD$sa, sb=optValuesD$sb, va=optValuesD$va, vb=optValuesD$vb, nseD=optDaily$value, elpTimeD=elpTimeD)
saveRDS(IHcoeffs, file = paste0(outLocationPath, "/IHcoeffs_final.rds"))


### Plot optimized IHACRES values
if(make_plots){
  jpeg(file=paste0(outLocationPath, "/", "IHACRES_Coeffs_Final.jpg")); par(mfrow = c(3,2))
  for (i in 2:7){
    val = optValuesD[1,i]
    plot(c(0, 1), c(0, 0), type = "n", xlab = "", ylab = "",
         main = colnames(optValuesD)[i], xlim = c(-0.2, 1.2), ylim = c(-1, 1))
    segments(0, 0, 1, 0, col = "black", lwd = 2)
    points(val, 0, col = "red", pch = 19, cex = 1.5)
    text(val, 0.3, sprintf("%.4f", val), col='red')
  }
  dev.off()
}


### Rerun entire model with optimal variables ###
#DailyWB<- WB(DailyClimData, gw_add, vfm, jrange,hock, hockros, dro, mondro, aspect,
#             slope, shade.coeff, jtemp, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
#DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)

### Save results as RDS file ###
saveRDS(results, file = paste0(outLocationPath, "/optim_results.rds"))