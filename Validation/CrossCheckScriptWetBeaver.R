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

###set path variables and source in Files ###
path<- "C:/Users/jcrane/repos/combined-water-balance/WaterBalance"
codePath <- file.path(path, 'Code')
outPath <- file.path(path,"Output")
dataPath <- file.path(path, "Data")

file.sources <- list.files(path=codePath,pattern="*.R")
setwd(codePath)
sapply(file.sources,source,.GlobalEnv)
setwd(path)

### define the data frame v to hold the variables ###
v = data.frame(
  #Overall variables
  FolderName = "ExcelValidation" , PETMethod = "Oudin", scrape = 0,optimization = 0 , 
  NonZeroDrainInitCoeff = 0, incompleteMonths = 1, GridMet = 1, 
  future = 0, fillLeapDays = 1, userSetJTemp = 0, plot = 1, delayStart = 0, percentRedGrid = 0.5, 
  
  #Climate Variables
  ClimateSiteID = "Wet Beaver Creek",lat=34.7,lon=-111.43,GaugeSiteID = "09505200",
  startY =1979, startM=01, startD=01, endY = 2023, endM = 06, endD = 24,
  
  #adjustment variables
  tmmx_slope = 1.013, tmmx_bias = -0.4422, tmmn_slope = 0.9647,tmmn_bias = -0.7425,
  p_slope = 0.8816,p_bias = 0.055,
  
  #Water Balance variables
  gw_add = 0.051,vfm = 0.453286145164077,jrange = 3, jtemp = 1.982841,hock =  4,hockros = 4,dro = 0,mondro = 0,
  aspect = 200,slope= 0,shade.coeff= 1,SWC.Max = 249,Soil.Init = 249, Snowpack.Init = 0,T.Base = 0,
  
  #Water Balance lower optimization limits
  l_gw_add = 0, l_vfm = 0.25, l_jrange = 2, l_hock = 0,l_hockros = 0,
  l_dro= 0, l_mondro = 0, l_aspect= 0, l_slope =  0, l_shade.coeff = 0.7,
  
  #Water Balance upper optimization limits
  u_gw_add = 0.1, u_vfm = 0.4, u_jrange = 4, u_hock =  8,u_hockros = 8,
  u_dro= 0.5,u_mondro= 0.5,u_aspect= 360, u_slope = 45, u_shade.coeff = 1)

#iterate over the variable data frame
#initialize counter for creating folders dynamically
counter = 1
with(v,WaterBalanceDrain(
  FolderName= FolderName, PETMethod = PETMethod,scrape = scrape,optimization = optimization, 
  NonZeroDrainInitCoeff = NonZeroDrainInitCoeff,incompleteMonths = incompleteMonths,GridMet = GridMet, 
  percentRedGrid = percentRedGrid, future = future, fillLeapDays = fillLeapDays, userSetJTemp = userSetJTemp, plot = plot, delayStart=delayStart,
  ClimateSiteID = ClimateSiteID, lat = lat, lon = lon,GaugeSiteID = GaugeSiteID,
  startY = startY, startM = startM, startD = startD,endY = endY, endM = endM, endD = endD, 
  tmmx_slope = tmmx_slope, tmmx_bias = tmmx_bias, tmmn_slope = tmmn_slope, tmmn_bias = tmmn_bias,
  p_slope = p_slope, p_bias =p_bias,
  gw_add = gw_add, vfm = vfm , jrange = jrange, jtemp=jtemp, hock = hock ,hockros = hockros,dro = dro,mondro = mondro ,
  aspect =  aspect,slope = slope,shade.coeff = shade.coeff ,SWC.Max = SWC.Max, Soil.Init = SWC.Max,
  Snowpack.Init =Snowpack.Init, T.Base = T.Base, 
  l_gw_add = l_gw_add, l_vfm = l_vfm,l_jrange = l_jrange, l_hock = l_hock,l_hockros = l_hockros,
  l_dro = l_dro, l_mondro = l_mondro, l_aspect = l_aspect, l_slope = l_slope,l_shade.coeff = l_shade.coeff, 
  u_gw_add = u_gw_add, u_vfm = u_vfm, u_jrange = u_jrange, u_hock = u_hock,u_hockros = u_hockros, 
  u_dro = u_dro,u_mondro = u_mondro,u_aspect = u_aspect, u_slope = u_slope, u_shade.coeff = u_shade.coeff,
  counter = counter))
  
#Read in the RDS file saved by the function above
outLocationPathEV = file.path(outPath, "WetBeaverCreek", "ExcelValidation")
results = readRDS(file = paste0(outLocationPathEV, "/results.rds"))
