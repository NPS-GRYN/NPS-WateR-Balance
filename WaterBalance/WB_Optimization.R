# ---------------------------------------------------------------------
# This script contains code for optimizing the water balance model by comparing
# modeled ET to OpenET estimates of ET. This optimization can occur either for potential 
# ET (PET) or actual ET (AET); code for both are included here. 
# This script is not intended to be run independently; it must be called from 
# the main Run_WBModel.R script.
# ---------------------------------------------------------------------


### Create file for results ###
results <- data.frame(SiteID = SiteID, start = startDate, end = endDate, PETMethod = PETMethod, optimization = optimization,
                      GridMET = GridMET, lon = lon, lat = lat,
                      startY = startY, startM = startM, startD = startD, endY = endY, endM = endM, endD = endD,
                      incompleteMonths = incompleteMonths, optimization_var = optimization_var)

parms<- c(gw_add = gw_add, vfm = vfm, jrange = jrange, hock =  hock, hockros = hockros,dro = dro, mondro = mondro,
          aspect = aspect,slope= slope, shade.coeff= shade.coeff, SWC.Max = SWC.Max, jtemp = jtemp, 
          et_slope=et_slope, et_bias=et_bias)

# Run optimization routine
strtTimeM <-Sys.time()
set.seed(123) 
WBcoeffs <- tibble()

### Run optimization routine ###
# use GA on a monthly time scale
# possible parameters: popSize=200, maxiter=100, pmutation = 0.2, pcrossover = 0.8, elitism = 10, run = 50, monitor = TRUE 
optMonth_init <- ga(type = "real-valued", fitness = function(x)
  WB_Optim_ET(c(gw_add=x[1], vfm=x[2], jrange=x[3], hock=x[4], hockros=x[5], dro=x[6], mondro=x[7], aspect=x[8], slope=x[9], shade.coeff=x[10], SWC.Max=x[11], jtemp=x[12], et_slope=x[13], et_bias=x[14]),
               Soil.Init = Soil.Init, Snowpack.Init = Snowpack.Init, T.Base = T.Base, PETMethod= PETMethod,
               DailyClimData = DailyClimData, lat=lat, lon=lon, meas_et = MonthlyET, interval='monthly', optimization_var=optimization_var),
  lower=WB_lower, upper=WB_upper, popSize=200, maxiter=100, maxFitness = 0.98, pmutation = 0.2, pcrossover = 0.8, elitism = 10, run = 50)
elpTimeM <- Sys.time() - strtTimeM


# Define the water balance variables from the best run
optValuesM <- data.frame(nseM = optMonth_init@fitnessValue, optMonth_init@solution)
gw_add=optValuesM$gw_add; vfm=optValuesM$vfm; jrange=optValuesM$jrange; hock=optValuesM$hock
hockros=optValuesM$hockros; dro=optValuesM$dro; mondro=optValuesM$mondro; aspect=optValuesM$aspect
slope=optValuesM$slope; shade.coeff=optValuesM$shade.coeff; SWC.Max=optValuesM$SWC.Max; jtemp=optValuesM$jtemp
et_slope=optValuesM$et_slope; et_bias=optValuesM$et_bias

# store and save results
results = data.frame(results, optValuesM, elpTimeM = elpTimeM)
saveRDS(WBcoeffs, file = paste0(outLocationPath, "/WBcoeffs.rds"))

# print time
print(paste("Time elapsed:", elpTimeM))

### PLOTS
if(make_plots){
  if (dev.cur() != 1) dev.off()
  # Parallel coordinates plot of WB coeffs
  jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_ParallelCoords.jpg"), width=1600, height=800, res=100)
  print(ggparcoord(data=WBcoeffs, columns=1:14, groupColumn=15, scale="uniminmax") + 
          scale_color_gradient(low = "black", high = "gray90") + nps_theme()) # for red: "darkred" and "#fee5d9"
  dev.off()
  
  # Scatterplots of WB coeffs
  WBcoeffs_long <- reshape(WBcoeffs, varying = names(WBcoeffs)[1:14], v.names = "WBcoeffs", timevar = "Variable", 
                           times = names(WBcoeffs)[1:14], direction = "long")
  jpeg(file=paste0(outLocationPath, "/", "WB_Coeffs_Scatter.jpg"), width=800, height=500)
  print(ggplot(WBcoeffs_long, aes(x = WBcoeffs, y = nse)) + geom_point() +
          facet_wrap(~ Variable, scales = 'free') + nps_theme() +
          labs(title = 'Water Balance Coefficients', x='', y = 'Monthly NSE'))
  dev.off()
  
  # Range plot of optimal WB coeffs
  wb_optim <- data.frame(var=c('Groundwater Addition', 'Volume Forcing Multiplier', 'Jennings Temperature Range','Hock','Hock Rain on Snow','Direct Runoff','Mondro','Aspect','Slope','Shade Coefficient','Max Soil Water Content','ET Slope','ET Bias','Jennings Temperature'),
                         value=c(gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope, shade.coeff, SWC.Max, et_slope, et_bias, jtemp),
                         lower=WB_lower, upper=WB_upper)
  jpeg(file=paste0(outLocationPath, "/", "WB_Optim_Coefficients.jpg"), width=800, height=600); par(mfrow = c(3, 5))
  for (i in 1:14){
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