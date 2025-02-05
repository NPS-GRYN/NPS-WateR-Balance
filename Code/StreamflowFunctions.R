# ---------------------------------------------------------------------
# This script includes functions that support the calibration and running of the Water Balance
# and IHACRES flow models.
# 
# EDITS IN PROGRESS
# improve documentation + code readability
# ---------------------------------------------------------------------

# Theme for plots
windowsFonts("Frutiger LT Std 55 Roman" = windowsFont("Frutiger LT Std 55 Roman"))
nps_theme <- function(base_size = 20, base_family="Frutiger LT Std 55 Roman") {
  theme_bw(base_size = base_size, base_family = "Frutiger LT Std 55 Roman") %+replace%
    theme(axis.text.x = element_text(family="Frutiger LT Std 55 Roman", size = base_size * 0.8), complete = TRUE)}


# Check validity of date vector
# Args:
#   vec: vector of dates
# Returns:
#   Warning if the vector is not in the format year, month, day
CheckVecDates= function(vec){
  if(length(vec)!=3|names(vec)[1]!= "year" |names(vec)[2]!= "month"|names(vec)[3]!="day"){
    warning("start and end must be vectors of length 3 with names c(year, month, day)."); stop()
  }
}

# Validate the number of columns in a data frame
# Args:
#   inputData
#   n_col
# Returns:
#   Warning if the number of columns does not match n_col
CheckNumCol = function(inputData, n_col){
  if(ncol(inputData)!= n_col | !is.data.frame(inputData)){
    warning(... = paste("Input should be a data frame and should have",n_col , "columns.")); stop()
  }
}

# Calculate IHACRES vb from other IHACRES coefficients, to ensure conservation of mass
# Args:
#   IHACRES coefficients (qa, qb, sa, sb, va)
# Returns: 
#   Calculated value of vb
calc_vb = function(qa, qb, sa, sb, va){
  vb <- (1-(qb/(1-qa)+sb/(1-sa)))*(1-va)
  return(vb)
}

# Calculate adjusted runoff and intermediate variables given meteorological data as input
# Args: 
#   DailyWB: dataframe with meteorological data @@
#   Water Balance parameters(gw_add, vfm, jrange ,hock ,hockros,dro,mondro , aspect, 
#   slope, shade.coeff, jtemp ,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod): parameters for water balance model
#   lat, lon: Latitude and longitude of site, in degrees
# Returns:
#   DailyWB dataframe with @@
WB= function(DailyWB, gw_add, vfm , jrange ,hock ,hockros,dro,mondro , aspect, 
             slope, shade.coeff, jtemp ,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon){
  #add month and day of year columns, GridMet comes with a date column, but the functions below need month and day of year
  DailyWB$month<- as.numeric(format(as.Date(DailyWB$date, format="%Y-%m-%d"),"%m"))
  DailyWB$yday<- yday(DailyWB$date)
  DailyWB$tmean_C<-apply(cbind(DailyWB$tmmn,DailyWB$tmmx), MARGIN = 1, mean)
  DailyWB$RAINDRO = get_dro(DailyWB = DailyWB, mondro = mondro, dro=dro, tmean =DailyWB$tmean_C, jtemp = jtemp,jrange = jrange, precip = DailyWB$pr, month = DailyWB$month)
  DailyWB$PRECIP = get_precip(precip = DailyWB$pr, raindro = DailyWB$RAINDRO)
  DailyWB$F = get_freeze(jtemp =jtemp , tmean = DailyWB$tmean_C,jrange = jrange)
  DailyWB$RAIN = get_rain(ppt = DailyWB$PRECIP,freeze =  DailyWB$F)
  DailyWB$SNOW = get_snow(ppt = DailyWB$PRECIP, freeze = DailyWB$F)
  DailyWB$MELT = get_melt(rain = DailyWB$RAIN,hockros = hockros, hock = hock, tmean = DailyWB$tmean_C, jtemp = jtemp, snow = DailyWB$SNOW, sp.0 = Snowpack.Init, jrange = jrange)
  DailyWB$PACK = get_snowpack(snow = DailyWB$SNOW, melt = DailyWB$MELT, sp.0 = Snowpack.Init)
  DailyWB$W = get_w(DailyWB$MELT, DailyWB$RAIN)
  DailyWB$PET = switch(PETMethod, 
                    Oudin=  {get_OudinPET(doy = DailyWB$yday, lat = lat,snowpack =  DailyWB$PACK, tmean = DailyWB$tmean_C,slope =  slope, aspect = aspect, shade.coeff = shade.coeff)},
                    Penman= {warning("Penman PET Method is not supported yet"); stop()},
                    Hamon = {warning("Hamon PET Method is not supported yet"); stop()})
  DailyWB$W_PET = get_w_pet(w = DailyWB$W, pet = DailyWB$PET)
  DailyWB$SOIL = get_soil(w = DailyWB$W, swc.0 = Soil.Init, pet = DailyWB$PET, swc.max = SWC.Max,w_pet = DailyWB$W_PET)
  DailyWB$DELTA_SOIL = get_d_soil(DailyWB$SOIL, swc.0 = Soil.Init)
  DailyWB$AET = get_AET(w = DailyWB$W, pet = DailyWB$PET, swc = DailyWB$SOIL, swc.0 = Soil.Init)
  DailyWB$RUNOFF = get_runoff(w = DailyWB$W, d_soil = DailyWB$DELTA_SOIL, AET = DailyWB$AET, RainDRO = DailyWB$RAINDRO)  # this is W - ET - DELTA_SOIL
  DailyWB$D = get_deficit(pet = DailyWB$PET, AET = DailyWB$AET)
  DailyWB$GDD = get_GDD(tmean = DailyWB$tmean_C, tbase = T.Base)
  DailyWB$adj_runoff = get_adj_runoff(orig = DailyWB$RUNOFF, gw_add = gw_add, vfm = vfm)
  return(DailyWB)
}

# Calculate total flow based on water balance output
# Args:
#   DailyDrain: dataframe containing @@
#   Initial IHACRES coefficients (q0, s0, v0): initial values for IHACRES model. Must be between 0 and 1 (inclusive)
#   IHACRES coefficients (qa, qb, sa, sb, va, vb): numeric values between 0 and 1 (inclusive) for IHACRES model 
# Returns:
#   Data frame including total streamflow and intermediate calculations of quick, slow, and very slow flow
Drain = function(DailyDrain, q0, s0, v0, qa, qb, sa, sb, va, vb){
  days <- nrow(DailyDrain)
  # Calculate flows for first day only, based on initial flow conditions
  q<- c(); q[1] = DailyDrain$adj_runoff[1]*qb+q0*qa
  s<- c(); s[1] = DailyDrain$adj_runoff[1]*sb+s0*sa
  v<- c(); v[1] = DailyDrain$adj_runoff[1]*vb+v0*va
  
  # Calculate flows for the rest of the days, based on previous day's flow
  for (i in 2:days){
    q[i] = DailyDrain$adj_runoff[i] * qb + q[i-1] * qa
    s[i] = DailyDrain$adj_runoff[i] * sb + s[i-1] * sa
    v[i] = DailyDrain$adj_runoff[i] * vb + v[i-1] * va
  }     
  DailyDrain$Quick<- q
  DailyDrain$Slow<- s
  DailyDrain$Very_Slow<- v
  DailyDrain$total<-q+s+v
  return(DailyDrain)
}

# Aggretate modeled and measured streamflow at a monthly scale
# Args:
#   DailyDrain: dataframe of modeled streamflow
#   meas_flow_mon: dataframe of measured streamflow
#   cutoffYear: year at which calibration to historical data begins
# Returns:
#   Dataframe with monthly summed measured and modeled streamflow
MeasModWB = function(DailyDrain, meas_flow_mon, cutoffYear){
  Mod<- DailyDrain[,c("date","total")]
  Mod$YrMon<- format(as.Date(Mod$date, format="%Y-%m-%d"),"%Y-%m")
  ModAgg<- aggregate(Mod$total, by=list(Mod$YrMon), FUN=sum)
  colnames(ModAgg)<- c("YrMon", "ModMM")
  MeasMod<- dplyr::full_join(meas_flow_mon, ModAgg, by = join_by(YrMon))
  colnames(MeasMod)<- c("YrMon", "Meas","Mod")
  MeasMod<- MeasMod[complete.cases(MeasMod),]
  MeasMod<- subset(MeasMod, MeasMod$YrMon>cutoffYear)
  row.names(MeasMod)<- NULL
  return(MeasMod)
}



# Calculate initial drainage coefficients by averaging January average flows for 10 year calibration period
# Generally not used because calibration is not very accurate; if not used, initial coefficients will be 0
# Args:
#   DailyClimData: dataframe of climate data including @@ at a daily time step 
#   Water Balance parameters(gw_add, vfm , jtemp , jrange ,hock ,hockros,dro, 
#   mondro , aspect,slope, shade.coeff, SWC.Max,Soil.Init, Snowpack.Init, T.Base,PETMethod): parameters for water balance model
#   Initial IHACRES coefficients (q0, s0, v0): initial values for IHACRES model. Must be between 0 and 1 (inclusive)
#   IHACRES coefficients (qa, qb, sa, sb, va, vb): numeric values between 0 and 1 (inclusive) for IHACRES model 
#   lat, lon: Latitude and longitude of site, in degrees
#   cutoffYear: year at which calibration to historical data begins
# Returns:
#   Dataframe (?) of initial drainage coefficients  
get_Init_Drain_Coef = function(DailyClimData, gw_add, vfm , jtemp , jrange ,hock ,hockros,dro,
                               mondro , aspect,slope, shade.coeff, SWC.Max,Soil.Init, Snowpack.Init, T.Base, PETMethod,
                               q0, s0, v0, qa, qb, sa, sb, va, vb, lat, lon, cutoffYear){

  DailyWB<- WB(DailyClimData, gw_add, vfm, jtemp, jrange, hock, hockros, dro, mondro, aspect, slope, 
               shade.coeff, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
  DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
  #add a YrMon Column to DailyDrain
  DailyDrain$YrMon<- format(as.Date(DailyDrain$date, format="%Y-%m-%d"),"%Y-%m")
  
  #average the quick, slow, and very slow flow by month
  QuickFlowMon<- aggregate(DailyDrain$Quick, by=list(DailyDrain$YrMon), FUN=mean)
  SlowFlowMon<- aggregate(DailyDrain$Slow, by=list(DailyDrain$YrMon), FUN=mean)
  VSFlowMon<- aggregate(DailyDrain$Very_Slow, by=list(DailyDrain$YrMon), FUN=mean)
  Mon<- aggregate(DailyDrain$month, by=list(DailyDrain$YrMon), FUN=unique)
  
  #Change the Column names so the merge works
  colnames(QuickFlowMon)<- c("YrMon", "Quick")
  colnames(SlowFlowMon)<- c("YrMon", "Slow")
  colnames(VSFlowMon)<- c("YrMon", "Very_Slow")
  colnames(Mon)<- c("YrMon", "Month")

  #merge all by YrMon
  df_list <- list(QuickFlowMon, SlowFlowMon, VSFlowMon, Mon)
  FlowMon<- Reduce(function(x, y) merge(x, y, all=TRUE), df_list)

  #subset by Jan
  JanFlow <- subset(FlowMon, FlowMon$Month==1)
  JanFlow <- subset(JanFlow, JanFlow$YrMon>cutoffYear)

  #apply the mean over the Flow columns to get the average
  InitConditions<- apply(JanFlow[,c("Quick","Slow","Very_Slow")], MARGIN = 2, FUN = mean)
  return(InitConditions)
}

# Water Balance function for first optimization routine
# Args:
#   parms: list of water balance parameters for optimizing (gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope, shade.coeff, SWC.Max, jtemp)
#   meas_flow_daily_xts: xts object containing measured streamflow at daily time step
#   cutoffYear: year at which calibration to historical data begins
#   Initial IHACRES coefficients (q0, s0, v0): initial values for IHACRES model. Must be between 0 and 1 (inclusive)
#   IHACRES coefficients (qa, qb, sa, sb, va, vb): numeric values between 0 and 1 (inclusive) for IHACRES model 
#   Default Water Balance parameters (Soil.Init, Snowpack.Init, T.Base, PETMethod): water balance parameters that are not included in optimization 
#   DailyClimData: 
#   lat, lon: Latitude and longitude of site, in degrees
#   meas_flow_mon: xts ? object containing measured streamflow at a monthly timestep
# Returns:
#   Monthly NSE of modeled streamflow
WB_Optim = function(parms,
                    #these are passed into WB_optim() from the global environment
                    meas_flow_daily_xts, cutoffYear, q0, s0, v0, qa, qb, sa, sb, va, vb,  
                    Soil.Init, Snowpack.Init, T.Base, PETMethod, DailyClimData,
                    lat, lon, meas_flow_mon){
  gw_add=parms[["gw_add"]]; vfm=parms[["vfm"]]; jrange=parms[["jrange"]];hock=parms[["hock"]];hockros = parms[["hockros"]]
  dro=parms[["dro"]]; mondro=parms[["mondro"]];aspect=parms[["aspect"]]; slope=parms[["slope"]];
  shade.coeff=parms[["shade.coeff"]]; SWC.Max=parms[["SWC.Max"]]; jtemp=parms[["jtemp"]]
  
  # run WB
  DailyWB <- WB(DailyClimData, gw_add, vfm , jrange, hock, hockros, dro, mondro, aspect, 
               slope, shade.coeff, jtemp, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
  
  # run drain
  DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)

  # calculate monthly NSE
  MeasMod <- MeasModWB(DailyDrain, meas_flow_mon, cutoffYear)
  x<- MeasMod$Mod
  y<- MeasMod$Meas
  nseM = NSE(x, y)
  Coeffs = data.frame(gw_add=gw_add, vfm= vfm, jrange=jrange, hock =hock, hockros=hockros, dro =dro, mondro = mondro, 
                                  aspect = aspect, slope =slope, shade.coeff=shade.coeff, SWC.Max=SWC.Max, jtemp = jtemp, nseM=nseM)
  WBcoeffs <<- rbind(WBcoeffs, Coeffs) # (<<-) is a global assignment operator
  print(str_c("nseM ", round(nseM, 4)))
  return(nseM)
}



# Water Balance function for optimization of only water balance based on AET
# Args:
#   parms: list of water balance parameters for optimizing (gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect, slope, shade.coeff, SWC.Max, jtemp)
#   cutoffYear: year at which calibration to historical data begins
#   Default Water Balance parameters (Soil.Init, Snowpack.Init, T.Base, PETMethod): water balance parameters that are not included in optimization 
#   DailyClimData: 
#   lat, lon: Latitude and longitude of site, in degrees
#   meas_aet_mon: xts ? object containing measured AET at a monthly timestep, from OpenET
# Returns:
#   Monthly NSE of modeled streamflow
WB_Optim_AET = function(parms, Soil.Init, Snowpack.Init, T.Base, PETMethod, DailyClimData, lat, lon, meas_aet_mon){
  gw_add=parms[["gw_add"]]; vfm=parms[["vfm"]]; jrange=parms[["jrange"]];hock=parms[["hock"]];hockros = parms[["hockros"]]
  dro=parms[["dro"]]; mondro=parms[["mondro"]];aspect=parms[["aspect"]]; slope=parms[["slope"]];
  shade.coeff=parms[["shade.coeff"]]; SWC.Max=parms[["SWC.Max"]]; jtemp=parms[["jtemp"]]
  
  # run WB
  DailyWB <- WB(DailyClimData, gw_add, vfm , jrange, hock, hockros, dro, mondro, aspect, 
                slope, shade.coeff, jtemp, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)
  
  # Aggregate WB AET monthly 
  Mod <- DailyWB[,c("date","AET")]; colnames(Mod) <- c("Date", "Mod AET")
  #Mod$YrMon <- format(as.Date(Mod$date, format="%Y-%m-%d"),"%Y-%m"); meas_aet_mon$Date <- format(as.Date(meas_aet_mon$Date),"%Y-%m")
  #ModAgg <- aggregate(Mod$AET, by=list(Mod$YrMon), FUN=sum)
  #colnames(ModAgg) <- c("Date", "Mod AET")
  MeasMod <- dplyr::full_join(meas_aet_mon, Mod, by = join_by(Date))
  colnames(MeasMod)<- c("YrMon", "Meas","Mod")
  MeasMod<- MeasMod[complete.cases(MeasMod),]
  
  # calculate monthly NSE
  x<- MeasMod$Mod
  y<- MeasMod$Meas
  nseM = NSE(x, y)
  Coeffs = data.frame(gw_add=gw_add, vfm= vfm, jrange=jrange, hock =hock, hockros=hockros, dro =dro, mondro = mondro, 
                      aspect = aspect, slope =slope, shade.coeff=shade.coeff, SWC.Max=SWC.Max, jtemp = jtemp, nseM=nseM)
  WBcoeffs <<- rbind(WBcoeffs, Coeffs) # (<<-) is a global assignment operator
  print(str_c("nseM ", round(nseM, 4)))
  return(nseM)
}



# IHACRES Flow function for second optimization
# Args:
#   parms: list of numeric values between 0 and 1 (inclusive) corresponding to qa, qb, sa, sb, va
#   Initial IHACRES coefficients (q0, s0, v0): initial values for IHACRES model
#   Daily WB: 
#   meas_flow_daily_xts: xts object containing measured streamflow at daily time step
#   cutoffYear: year at which calibration to historical data begins
# Returns:
#   Daily NSE of modeled streamflow
IHACRESFlow <- function(parms, q0, s0, v0, DailyWB, meas_flow_daily_xts, cutoffYear){ 
  if(flow_components==3){
    qa=parms[['qa']]; qb=parms[['qb']]; sa=parms[['sa']]; sb=parms[['sb']]; va=parms[['va']]
  } else if(flow_components==2){
    qa=0; qb=0; sa=parms[['sa']]; sb=parms[['sb']]; va=parms[['va']]
  } else{print('invalid number of flow parameters')}
  
  # vb must be greater than 0 for conservation of mass
  # returning NA doesn't work with optim() so instead returns a very low NSE
  vb <- calc_vb(qa,qb,sa,sb,va)
  if(vb < 0 | is.na(vb)){
    print(paste('not valid', vb))
    return(-100000000)
  }
  #obtain the total/streamflow from Drain()
  DailyDrain <- Drain(DailyWB, q0, s0, v0, qa, qb, sa, sb, va, vb)
  
  # Calculate NSE of daily measured vs modeled streamflow
  Mod_xts<- xts(DailyDrain$total, ymd(DailyDrain$date))
  DrainMeasMod_xts<- merge.xts(meas_flow_daily_xts, Mod_xts, all=FALSE)
  df_d<- data.frame(date = index(DrainMeasMod_xts), total_mod = as.numeric(coredata(DrainMeasMod_xts$Mod_xts)), 
                    total_meas = as.numeric(coredata(DrainMeasMod_xts$meas_flow_daily_xts)))
  df_d<- df_d[complete.cases(df_d),]#the measured data has NA's because it is missing values
  df_d = df_d[df_d$date >= paste(cutoffYear, "01", "01", sep = "-"), ]
  nseD <- NSE(sim = df_d$total_mod,obs = df_d$total_meas)
  Coeffs = data.frame(qa=qa, qb=qb, sa=sa, sb=sb, va=va, vb=vb, nseD=nseD)
  IHcoeffs <<- rbind(IHcoeffs, Coeffs) 
  print(str_c("nseD ", round(nseD, 4)))
  return(nseD)
} 

