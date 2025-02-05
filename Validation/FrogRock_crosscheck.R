# ---------------------------------------------------------------------
# This file contains the code to run the water balance model and compare to a calibrated excel sheet
# Currently the code will only work for pre-calibrated streamflow projections at Frog Rock.
# 
# EDITS IN PROGRESS
# compare optimized code to spreadsheet model
# add crosscheck of streamflow model (WB + IHACRES)
# develop "worksheet" where it is easy for users to play with parameters and see impact on streamflow
# go through the rest of the crosscheck folder and figure out what is worth keeping/what's not
# ---------------------------------------------------------------------

#######################################################################
### Load libraries ###
library(sf); library(raster); library(ggplot2); library(dplyr); library(xts)
library(lubridate); library(hydroGOF); library(stringr); library(terra); library(glue); library(tidyverse)
library(climateR); library(EGRET); library(daymetr); library(here); library(ggrepel); library(gridExtra); 
library(httr); library(jsonlite); library(sf); library(grid); library(GA); library(GGally)

### Source in function files ###
path <- here() 
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())

#######################################################################
### Set model variables ###
# Water Balance variables
gw_add=0; vfm = 1; jtemp = 2.1555781; jrange = 3; hock = 4; hockros = 4; 
dro = 0; mondro = 0; aspect = 4.236395; slope= 1.9382907; shade.coeff= 1; 
SWC.Max = 104; Soil.Init = SWC.Max; Snowpack.Init = 0; T.Base = 0 
PETMethod='Oudin'

# Location data
lat=44.95354; lon=-110.54083; elev = 2136

# create start and end date objects of data collection.
startDate<- ymd(paste(1980, 01, 01))
endDate<-  ymd(paste(1981, 12, 31))


#######################################################################
### GET INPUT DATA ### 
DailyClimData <- read.csv(here("RExcelCrossCheckWetBeaver/FrogRock_Input.csv"))

# Match format of GridMET data
DailyClimData$date <- make_date(DailyClimData$year, 1) + (DailyClimData$yday-1)
DailyClimData$month<- as.numeric(format(as.Date(DailyClimData$date, format="%Y-%m-%d"),"%m"))
DailyClimData$year<- NULL
DailyClimData$yday<- NULL
DailyClimData$swe..kg.m.2.<- NULL
colnames(DailyClimData)<- c("dayl..s." , "pr","srad" ,"tmmx", "tmmn", "vp..Pa." ,"date", "month")


#######################################################################
### RUN WATER BALANCE ###
DailyWB_R <- WB(DailyClimData, gw_add, vfm, jrange, hock, hockros, dro, mondro, aspect,
             slope, shade.coeff, jtemp, SWC.Max, Soil.Init, Snowpack.Init, T.Base, PETMethod, lat, lon)

write.csv(DailyWB_R, here('RExcelCrossCheckWetBeaver/FrogRock_R_Output.csv'), row.names = FALSE) 


#######################################################################
### COMPARE TO EXCEL AND PYTHON OUTPUT ###
DailyWB_Excel <- read.csv(here('RExcelCrossCheckWetBeaver/FrogRock_Excel_Output.csv'))
DailyWB_Python <-read.csv(here('RExcelCrossCheckWetBeaver/FrogRock_Python_Output.csv'))

### PLOT ###
wb_vars <- c('SOIL','AET','D','PET','PACK','SNOW')
if(!dir.exists(here('Output', 'Validation'))) {dir.create(here('Output', 'Validation'))}

# Compare Excel and Python 
jpeg(here('Output','Validation','ExcelPython_Compare.jpg'), width=600, height=400); par(mfrow=c(2,3))
for(var in wb_vars){
  print(var)
  # calculate R2
  model <- lm(DailyWB_Python[[var]] ~ DailyWB_Excel[[var]])
  r2_value <- summary(model)$r.squared
  
  # plot
  plot(DailyWB_Excel[[var]], DailyWB_Python[[var]], 
       xlab = paste('Excel', var), ylab = paste('Python', var), 
       main = paste(var, "R² = ", round(r2_value, 2))) + nps_theme()
}
dev.off()

# Compare R and Excel
jpeg(here('Output','Validation','ExcelR_Compare.jpg'), width=600, height=450); par(mfrow=c(2,3))
for(var in wb_vars){
  print(var)
  # calculate R2
  model <- lm(DailyWB_R[[var]] ~ DailyWB_Excel[[var]])
  r2_value <- summary(model)$r.squared
  
  # plot
  plot(DailyWB_Excel[[var]], DailyWB_R[[var]], 
       xlab = paste('Excel', var), ylab = paste('R', var), 
       main = paste(var, "R² = ", round(r2_value, 2))) + nps_theme()
}
#mtext("Comparison of Excel vs R Water Balance", side = 3, line=-1, outer = TRUE, cex = 1.5)
dev.off()





