#working folder
setwd("/Users/marioguevara/Downloads")
#required libraries
library(raster)
#read the covariates 1km
covar <- stack("MEX_worldgridsCOVS.tif")
#define name for covariates
names(covar) <- readRDS('worldgridsCOVS_names.rds')
