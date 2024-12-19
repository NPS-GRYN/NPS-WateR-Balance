## GENERATE FUTURE PROJECTIONS OF STREAMFLOW ##
library(ggrepel)


fut_ro1 <- read.csv(file.path(dataPath, paste(GaugeSiteID, "watershed_avg_water_balance_future.csv", sep = "_"))) 
gcms<-unique(fut_ro1$GCM) #all models

#drop "MIROC-ESM-CHEM.rcp85" because it doesn't have an associated RCP 4.5
fut_ro1<-subset(fut_ro1, GCM !="MIROC-ESM-CHEM.rcp85")
gcms<-unique(fut_ro1$GCM)

#use the first model type to grab the dates
fut_ro_date<-subset(fut_ro1, GCM == gcms[1]);head(fut_ro_date);tail(fut_ro_date)
fut_date<-as.Date(fut_ro_date$Date);fut_date#format needed for working with xts
df_date<-as.data.frame(fut_date)#format needed for use with data frames

#define the length of the future data in days
d<-length(fut_date)
futures<-NULL #empty object to hold data after modeled runoff converted to modeled stream flow

# j iterates over the models, i is being defined in the drain function and iterates over days
for (j in 1:length(gcms)){
  #subset for one model
  fut_ro<-subset(fut_ro1, GCM == gcms[j])
  print(gcms[j])
  
  #convert to mm and fix column names
  fut_ro<-cbind(fut_ro[,c("Date","GCM")], 25.4*(fut_ro[,c(which( colnames(fut_ro)=="Deficit.in"):ncol(fut_ro))])) #grabs all columns to the from deficit to the end
  colnames(fut_ro)<- c("Date", "GCM", "Deficit", "AET", "soil_water", "runoff", "rain", "accumswe", "PET")
  
  #adjust for ground water addition and volume forcing multiplier
  fut_ro$adj_runoff<- get_adj_runoff(fut_ro$runoff, gw_add = gw_add, vfm = vfm)
  
  #create a data frame for the Drain function
  data = data.frame(fut_ro$adj_runoff)
  colnames(data) = c("adj_runoff")
  
  #run the Drain function
  DailyDrainFuture <- Drain(DailyWB = data, q0 = q0, qa = qa, qb = qb, s0 = s0, sa = sa, sb = sb, v0 = v0, va = va, vb = vb)
  
  #save this model run to futures
  drainage_qsvt <- cbind(gcms[j],df_date,DailyDrainFuture)
  colnames(drainage_qsvt)[] <- c("model","date","adj_runoff","quick","slow","veryslow","total")
  futures<-rbind(futures, drainage_qsvt)#stack results from each model into a heap for plotting
}
#filter out futures that are overlapping with the date of the historic flow
futures<- futures[futures$date>endDate,]

#add the yr_mo and yr columns to futures
futures$yr<-format(as.Date(futures$date),"%Y")
futures$mo<-format(as.Date(futures$date),"%m")
futures$yr_mo<-format(as.Date(futures$date),"%Y-%m")

# reorder columns to match order of column in projections  
futures<-futures[,c("model", "date", "yr", "mo", "yr_mo", 
                    "adj_runoff", "quick", "slow", "veryslow", "total")]

#extract gcm and rcp and create compiled
base<-strsplit(futures$model, split=".rcp")
gcm <- sapply(X = base, FUN = "[", 1) 
rcp <- as.numeric(x = sapply(base, FUN = "[", 2)) 
future<-cbind(gcm, rcp, futures[,which(colnames(futures)=="date"):ncol(futures)]) 
colnames(DailyHistFlow)[colnames(DailyHistFlow) == 'model'] <- 'gcm'
DailyHistFlow$rcp<-NA
#reorder the columns to have gcm, rcp, and then date through total
DailyHistFlow<-DailyHistFlow[,c(which(colnames(DailyHistFlow)=="gcm"), 
                                which(colnames(DailyHistFlow)=="rcp"),
                                which(colnames(DailyHistFlow)=="date"):which(colnames(DailyHistFlow)=="total"))] 
compiled<-rbind(DailyHistFlow, future)

#daily data frame
hist<-subset(compiled, gcm == "Historic")
hist$rcp<-"Hist"
fut<-subset(compiled, gcm!="Historic")

#delete then rebuild compiled with rcp values
daily_df<-as.data.frame(rbind(hist,fut))
daily_df$Period<-ifelse(daily_df$yr<=2022,"Historical",
                                   ifelse (daily_df$yr>=2023 & daily_df$yr<=2050,"Early",
                                           ifelse (daily_df$yr>=2051 & daily_df$yr<=2070,"Middle",
                                                   ifelse (daily_df$yr>=2071, "Late","NA"))))
daily_df$yr<-as.numeric(daily_df$yr)

#annual data frame
annual_df <- as.data.frame(daily_df %>% group_by(gcm, rcp, yr) %>%
  dplyr::summarize(adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                   slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE), total = sum(total, na.rm = TRUE), 
                   Period=first(Period)))
#create a time index by forcing month and day to be Jan 1st
annual_df$date<-as.Date(paste(annual_df$yr,"-01", "-01",sep=""))

#monthly data frame
monthly_df <- as.data.frame(daily_df %>%
  group_by(gcm,rcp, yr_mo) %>%
  dplyr::summarize(adj_runoff = sum(adj_runoff, na.rm = TRUE),quick = sum(quick, na.rm = TRUE),
                   slow = sum(slow, na.rm = TRUE),veryslow = sum(veryslow, na.rm = TRUE) ,total = sum(total, na.rm = TRUE),
                   Period=first(Period)))
#create a time index by forcing month and day to be Jan 1st
monthly_df$date<-as.Date(paste(monthly_df$yr_mo,"-01",sep=""))


### Make Daily, Monthly, and Annual xts plots of stream flow
if(make_plots){
  #make Daily Stream Flow xts plot and save as pdf
  name = paste(ClimateSiteID, "Daily Stream Flow")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  pdf(file=paste0(outLocationPath, "/", nameReduce, ".pdf"))
  par(mfrow = c(1,1))
  plot<-ggplot(data = daily_df) + geom_line(aes(x=date, y = total, colour= rcp))+
    facet_wrap(~gcm)+ ylab("Stream flow (mm)") + xlab("Day")+
    ggtitle(name)
  dev.off()
  print(plot)
  
  #make monthly Stream Flow xts plot and save as pdf
  name = paste(ClimateSiteID, "Monthly Stream Flow")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  pdf(file=paste0(outLocationPath, "/", nameReduce, ".pdf"))
  par(mfrow = c(1,1))
  plot<-ggplot(data = monthly_df) + geom_line(aes(x=date, y = total, colour= rcp))+
    facet_wrap(~gcm)+ ylab("Stream flow (mm)") + xlab("Month")+
    ggtitle(name)
  dev.off()
  print(plot)
  
  #make annual Stream Flow xts plot and save as pdf
  name = paste(ClimateSiteID, "Annual Stream Flow")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  pdf(file=paste0(outLocationPath, "/", nameReduce, ".pdf"))
  par(mfrow = c(1,1))
  plot<-ggplot(data = annual_df) + geom_line(aes(x=date, y = total, colour= rcp))+
    facet_wrap(~gcm)+ ylab("Stream flow (mm)") + xlab("Year")+
    ggtitle(name)
  dev.off()
  print(plot)
  
  #  make annual plot where all gcms go on same plot
  name = paste(ClimateSiteID, "Combined Plot Annual Streamflow")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  pdf(file=paste0(outLocationPath, "/", nameReduce, ".pdf"))
  par(mfrow = c(1,1))
  plot<-ggplot(data = annual_df) + geom_line(aes(x=date, y = total, colour= gcm))+
    facet_wrap(~rcp)+ ylab("Stream flow (mm)") + xlab("Year")+
    ggtitle(name)
  dev.off()
  print(plot)
  
  # scatter plot of annual magnitude vs daily standard deviation
  name = paste(ClimateSiteID, "Streamflow Climate Future Scatterplot")
  nameReduce = gsub(pattern = " ",replacement = "_", x = name)
  pdf(file=paste0(outLocationPath, "/", nameReduce, ".pdf"))
  par(mfrow = c(1,1))
  annual_df$delta_annual_mm <- annual_df$total - mean((annual_df %>% filter(gcm=='Historic'))$total)
  delta_plot <- daily_df %>% group_by(gcm, rcp, yr) %>% 
    dplyr::summarize(Period=first(Period), delta_daily_sd = sd(total)-sd((daily_df %>% filter(gcm=='Historic'))$total))
  delta_plot <- delta_plot %>% left_join(annual_df %>% select(gcm, rcp, yr, delta_annual_mm, Period), 
                                         by = c("gcm", "rcp", "yr","Period"))
  delta_plot <- delta_plot %>% filter(gcm!='Historic') %>% group_by(gcm,rcp,Period) %>% dplyr::summarise(delta_daily_sd=mean(delta_daily_sd), delta_annual_mm=mean(delta_annual_mm))
  annual_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                                xintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_annual_mm, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_annual_mm, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_annual_mm, 0.5)))
  sd_quantile <- data.frame(Period = c("Early", "Middle", "Late"),
                            yintercept = c(quantile((delta_plot %>% filter(Period=='Early'))$delta_daily_sd, 0.5), quantile((delta_plot %>% filter(Period=='Middle'))$delta_daily_sd, 0.5), quantile((delta_plot %>% filter(Period=='Late'))$delta_daily_sd, 0.5)))
  annual_zero <- data.frame(Period = c("Early", "Middle", "Late"), xintercept=c(0,0,0)); sd_zero <- data.frame(Period = c("Early", "Middle", "Late"), yintercept = c(0,0,0))
  plot <- ggplot(data=delta_plot, aes(x=delta_annual_mm,y=delta_daily_sd,color=rcp)) + geom_point() +
    geom_text_repel(aes(label = gcm), color = 'black', max.overlaps=Inf) +
    geom_hline(data = sd_quantile, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_quantile, aes(xintercept = xintercept), color = "black") +
    #geom_hline(data = sd_zero, aes(yintercept = yintercept), color = "black") + geom_vline(data = annual_zero, aes(xintercept = xintercept), color = "black") +
    facet_wrap(~factor(Period, levels = c('Early', 'Middle', 'Late'))) + labs(title=paste('Changes in streamflow at',ClimateSiteID), x='Change in annual flow magnitude [mm]', y='Change in daily flow standard deviation [mm]', color='RCP') + 
    scale_color_manual(values = c("45" = "orange", "85" = "red")) + nps_theme()
  dev.off()
  print(plot)
}