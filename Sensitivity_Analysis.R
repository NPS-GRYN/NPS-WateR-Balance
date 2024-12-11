library(sf)
library(raster)
library(ggplot2)
library(WaterBalance)
library(dplyr)
library(xts)
library(lubridate)
library(hydroGOF)
library(stringr)
library(terra)
library(glue)
library(tidyverse)
library(climateR)
library(EGRET)
library(daymetr)

### Set path variables and source in files ###
path<- "C:\\David\\Water balance\\Operational version\\Version 3\\USGS Basins comparison\\WaterBalance"
codePath <- file.path(path, 'Code')
outPath <- file.path(path,"Output")
dataPath <- file.path(path, "Data")

file.sources <- list.files(path=codePath,pattern="*.R")
setwd(codePath)
sapply(file.sources,source,.GlobalEnv)
setwd(path)

### Set specific location path and read in RDS files of previous run ###
outLocationPath = file.path(outPath, "Redwood", "DaymetOptimized0.1FullGrid")
results = readRDS(file = paste0(outLocationPath, "/results.rds"))
WBcoeffs = readRDS(file = paste0(outLocationPath, "/WBcoeffs.rds"))
IHcoeffs = readRDS(file = paste0(outLocationPath, "/IHcoeffs.rds"))

colnames(IHcoeffs)<- toupper(colnames(IHcoeffs))
colnames(WBcoeffs)<- toupper(colnames(WBcoeffs))

for(i in 1:(ncol(IHcoeffs)-1)){
  pdf(file=paste0(outLocationPath, "/NSEDaily",colnames(IHcoeffs)[i] ,"ScatterPlot.pdf"))
  plot(IHcoeffs[,i], IHcoeffs[,ncol(IHcoeffs)], xlab = colnames(IHcoeffs)[i], ylab = "NSE Daily", main = "NSE Daily")
  dev.off()
  plot(IHcoeffs[,i], IHcoeffs[,ncol(IHcoeffs)], xlab = colnames(IHcoeffs)[i], ylab = "NSE Daily", main = "NSE Daily")
  
}

for(i in 1:(ncol(WBcoeffs)-1)){
  pdf(file=paste0(outLocationPath, "/NSEMonthly",colnames(WBcoeffs)[i] ,"ScatterPlot.pdf"))
  plot(WBcoeffs[,i], WBcoeffs[,ncol(WBcoeffs)], xlab = colnames(WBcoeffs)[i], ylab = "NSE Monthly", main = "NSE Monthly")
  dev.off()
  plot(WBcoeffs[,i], WBcoeffs[,ncol(WBcoeffs)], xlab = colnames(WBcoeffs)[i], ylab = "NSE Monthly", main = "NSE Monthly")
  
}
