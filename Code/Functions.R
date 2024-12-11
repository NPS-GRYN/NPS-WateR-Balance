#this function checks the validity of date vectors inputted by the user used to make dates
#it used used inside get_elev_Daymet
CheckVecDates= function(vec){
  if(length(vec)!=3|names(vec)[1]!= "year" |names(vec)[2]!= "month"|names(vec)[3]!="day"){
    warning("start and end must be vectors of length 3 with names c(year, month, day).") ;stop()
  }
}

#this function validates the number of columns in a data frame
CheckNumCol = function(inputData, n_col){
  if(ncol(inputData)!= n_col | !is.data.frame(inputData)){
    warning(... = paste("Input should be a data frame and should have",n_col , "columns."))
    stop()
  }
}

#starting with meteorological data such as temperature, precip, and other variables, 
#returns a dataframe with all intermediate calculations up to adjusted runoff
WB= function(DailyClimData, gw_add, vfm , jrange ,hock ,hockros,dro,mondro , aspect, 
             slope, shade.coeff, jtemp ,SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon){
  DailyWB = DailyClimData
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
  DailyWB$DSOIL = get_d_soil(DailyWB$SOIL, swc.0 = Soil.Init)
  DailyWB$AET = get_AET(w = DailyWB$W, pet = DailyWB$PET, swc = DailyWB$SOIL, swc.0 = Soil.Init)
  DailyWB$W_ET_DSOIL = get_runoff(w = DailyWB$W, d_soil = DailyWB$DSOIL, AET = DailyWB$AET, RainDRO = DailyWB$RAINDRO)
  DailyWB$D = get_deficit(pet = DailyWB$PET, AET = DailyWB$AET)
  DailyWB$GDD = get_GDD(tmean = DailyWB$tmean_C, tbase = T.Base)
  DailyWB$adj_runoff = get_adj_runoff(orig = DailyWB$W_ET_DSOIL, gw_add = gw_add, vfm = vfm)
  return(DailyWB)
}

#This function takes the output from WB and returns a data frame with all intermediate calculations through total flow
#It breaks the adjusted runoff into pieces that flow at different times, quick, slow, and very slow. 
#Total flow is the sum of those on each day. Total flow is synonymous with Stream Flow
#The function does the calculations in two parts. The first part evaluates flow for the first day, the next for days 2 to the end
#It is done this way becayse there are initial flow conditions for the first day of flow that replaced by the previous days flow in subsequent calculations

Drain = function(DailyWB, q0, qa, qb, s0, sa, sb, v0, va, vb){
  DailyDrain = DailyWB
  days <- nrow(DailyDrain)
  #Calculate flows for first day only
  q<- c()
  q[1] = DailyDrain$adj_runoff[1]*qb+q0*qa
  
  s<- c()
  s[1] = DailyDrain$adj_runoff[1]*sb+s0*sa
  
  v<- c()
  v[1] = DailyDrain$adj_runoff[1]*vb+v0*va
  
  #calculate flows for the rest of the days
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

#This function is called MeasModWB because it compares the Measures vs modeled stream/ total flow 
#returns a data frame with monthly summed measured and modeled stream/total
MeasModWB = function(DailyDrain, MeasAgg, cutoffYear){
  Mod<- DailyDrain[,c("date","total")]
  Mod$YrMon<- format(as.Date(Mod$date, format="%Y-%m-%d"),"%Y-%m")
  ModAgg<- aggregate(Mod$total, by=list(Mod$YrMon), FUN=sum)
  colnames(ModAgg)<- c("YrMon", "ModMM")
  MeasMod<- dplyr::full_join(MeasAgg, ModAgg, by = join_by(YrMon))
  colnames(MeasMod)<- c("YrMon", "Meas","Mod")
  MeasMod<- MeasMod[complete.cases(MeasMod),]
  MeasMod<- subset(MeasMod, MeasMod$YrMon>cutoffYear)
  row.names(MeasMod)<- NULL
  return(MeasMod)
}

#this function returns new initial drainage coefficients averaging the January average flow 
#This is not usually used, but is done to fix the the inaccuracies in the first 10 or 20 years of modeled stream flow
#They are inaccurate because it takes 10-20 years to calibrate up to the correct values of flow
#If not used, these Initial Drain Coefficients will all be 0
get_Init_Drain_Coef = function(DailyClimData, gw_add, vfm , jtemp , jrange ,hock ,hockros,dro,
                               mondro , aspect,slope, shade.coeff, SWC.Max,Soil.Init, Snowpack.Init, T.Base, DailyWB, 
                               q0, qa, qb, s0, sa, sb, v0, va, vb, PETMethod, lat, lon, cutoffYear){

  DailyWB<- WB(DailyClimData = DailyClimData, gw_add = gw_add, vfm = vfm, jtemp = jtemp, jrange = jrange,
               hock = hock, hockros = hockros, dro = dro, mondro = mondro, aspect = aspect, slope = slope, 
               shade.coeff = shade.coeff, SWC.Max = SWC.Max, Soil.Init = Soil.Init, Snowpack.Init = Snowpack.Init, T.Base = T.Base, PETMethod = PETMethod, lat = lat, lon = lon)
  DailyDrain <- Drain(DailyWB = DailyWB, q0 = q0, qa = qa, qb = qb, s0 = s0, sa = sa, sb = sb, v0 = v0, va = va, vb = vb)
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
  JanFlow<- subset(FlowMon, FlowMon$Month==1)
  JanFlow<- subset(JanFlow, JanFlow$YrMon>cutoffYear)

  #apply the mean over the Flow columns to get the average
  InitConditions<- apply(JanFlow[,c("Quick","Slow","Very_Slow")], MARGIN = 2, FUN = mean)
  return(InitConditions)
}

#This is used in the first optimization routine and returns the NSE of stream flows on a monthly time step
#This code is in essence the same as the code in WB(), Drain(), and MeasModWB()
#If you make changes to those functions, paste the new code below in between the correct lines
WB_Optim = function(parms,
                    #these are passed into WB_optim() from the global environment
                    Meas_xts, cutoffYear, q0, s0, v0, qa, qb, sa, sb,va,
                    SWC.Max, Soil.Init, Snowpack.Init, T.Base, DailyClimData,
                    PETMethod, lat, lon, MeasAgg){
  gw_add = parms[["gw_add"]];vfm = parms[["vfm"]]; jrange = parms[["jrange"]];hock = parms[["hock"]];hockros = parms[["hockros"]]
  dro = parms[["dro"]];mondro = parms[["mondro"]];aspect = parms[["aspect"]];slope = parms[["slope"]];
  shade.coeff = parms[["shade.coeff"]];jtemp = parms[["jtemp"]]
  
  ###############   paste WB() here    ###################
  DailyWB = DailyClimData
  
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
  DailyWB$DSOIL = get_d_soil(DailyWB$SOIL, swc.0 = Soil.Init)
  DailyWB$AET = get_AET(w = DailyWB$W, pet = DailyWB$PET, swc = DailyWB$SOIL, swc.0 = Soil.Init)
  DailyWB$W_ET_DSOIL = get_runoff(w = DailyWB$W, d_soil = DailyWB$DSOIL, AET = DailyWB$AET, RainDRO = DailyWB$RAINDRO)
  DailyWB$D = get_deficit(pet = DailyWB$PET, AET = DailyWB$AET)
  DailyWB$GDD = get_GDD(tmean = DailyWB$tmean_C, tbase = T.Base)
  DailyWB$adj_runoff = get_adj_runoff(orig = DailyWB$W_ET_DSOIL, gw_add = gw_add, vfm = vfm)
  
  ###############   paste Drain() here    ###################
  DailyDrain = DailyWB
  days <- nrow(DailyDrain)
  vb = (1-(qb/(1-qa)+sb/(1-sa)))*(1-va)
  q<- c()
  q[1] = DailyDrain$adj_runoff[1]*qb+q0*qa

  s<- c()
  s[1] = DailyDrain$adj_runoff[1]*sb+s0*sa

  v<- c()
  v[1] = DailyDrain$adj_runoff[1]*vb+v0*va

  for (i in 2:days){
    q[i] = DailyDrain$adj_runoff[i] * qb + q[i-1] * qa

    s[i] = DailyDrain$adj_runoff[i] * sb + s[i-1] * sa

    v[i] = DailyDrain$adj_runoff[i] * vb + v[i-1] * va
  }
  DailyDrain$Quick<- q
  DailyDrain$Slow<- s
  DailyDrain$Very_Slow<- v
  DailyDrain$total<-q+s+v

  ###############   paste MeasModWB() here    ###################
  Mod<- DailyDrain[,c("date","total")]
  Mod$YrMon<- format(as.Date(Mod$date, format="%Y-%m-%d"),"%Y-%m")
  ModAgg<- aggregate(Mod$total, by=list(Mod$YrMon), FUN=sum)
  colnames(ModAgg)<- c("YrMon", "ModMM")
  MeasMod<- dplyr::full_join(MeasAgg, ModAgg, by = join_by(YrMon))
  colnames(MeasMod)<- c("YrMon", "Meas","Mod")
  MeasMod<- MeasMod[complete.cases(MeasMod),]
  MeasMod<- subset(MeasMod, MeasMod$YrMon>cutoffYear)
  row.names(MeasMod)<- NULL

  #################################################
  x<- MeasMod$Mod
  y<- MeasMod$Meas
  nseM = NSE(x, y)
  Coeffs = data.frame(gw_add=gw_add, vfm= vfm, jrange=jrange, hock =hock, hockros=hockros, dro =dro, mondro = mondro, 
                                  aspect = aspect, slope =slope, shade.coeff=shade.coeff, jtemp = jtemp, nseM=nseM)
  WBcoeffs<<- rbind(WBcoeffs, Coeffs) # (<<-) is a global assignment operator
  print(str_c("nseM ", round(nseM, 4)))
  return(nseM)
}

#This is the function used to perform the second optimization
#This function calculates the NSE of measured vs modeled stream flow on a daily basis
IHACRESFlow <- function(q0, s0, v0, qa, qb, sa, sb, va, DailyWB, Meas_xts, cutoffYear){      
  vb = (1-(qb/(1-qa)+sb/(1-sa)))*(1-va)
  #vb has to be greater than 0 for conservation of mass
  if(vb<0){ 
    nseD = NA
    return(nseD)
  }
  #obtain the total/stream flow from Drain()
  DailyDrain = DailyDrain <- Drain(DailyWB = DailyWB, q0 = q0, qa = qa, qb = qb, 
                                   s0 = s0, sa = sa, sb = sb, v0 = v0, 
                                   va = va,vb = vb)
  
  #merge the daily measured vs modeled stream flow and calculate the NSE of it
  Mod_xts<- xts(DailyDrain$total, ymd(DailyDrain$date))
  DrainMeasMod_xts<- merge.xts(Meas_xts,Mod_xts)
  df_d<- data.frame(date = index(DrainMeasMod_xts), total_mod = as.numeric(coredata(DrainMeasMod_xts$Mod_xts)), 
                    total_meas = as.numeric(coredata(DrainMeasMod_xts$Meas_xts)))
  df_d<- df_d[complete.cases(df_d),]#the measured data has NA's because it is missing values
  df_d = df_d[df_d$date >= paste(cutoffYear, "01", "01", sep = "-"), ]
  nseD <- NSE(sim = df_d$total_mod,obs = df_d$total_meas)
  coeffs = data.frame(qa=qa, qb= qb, sa=sa, sb =sb, va=va, vb  =vb, nseD = nseD)
  return(coeffs) 
}   

