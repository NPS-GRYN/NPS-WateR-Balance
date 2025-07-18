# ---------------------------------------------------------------------
# This script contains code for running and calibrating the streamflow model
# across multiple watersheds without additional user-provided input. 
# ---------------------------------------------------------------------

### Load libraries and function files ###
library('here'); lib_install <- FALSE
setwd(here('Code')); sapply(list.files(pattern="*.R"), source, .GlobalEnv); setwd(here())


#######################################################################
#######################################################################
### WATERSHED NAMES AND USGS GAGE NUMBERS ###
watershed_siteid <- c("Little River", "Cataloochee")
watershed_gageid <- c("03497300", "03460000")


#######################################################################
#######################################################################
### SET PARAMETERS FOR ALL WATERSHEDS ###
# These parameters will be used for all watersheds
# To experiement with the effects of parameter changes on calibration, use
# the Streamflow_Multi_Parameter.R script.