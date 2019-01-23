
#data preparation for estimating SOC values for the specific depth of 0-5cm

#working environment
  
setwd("~/Downloads/")

#load datasets

dat <- read.csv('V1_Holmquist_2018_depth_series_data.csv')
core <- read.csv("V1_Holmquist_2018_core_data.csv")

#merge datasets 

dat <- merge(dat, core)

#generate an id by coordinates

dat$IDPROF <- paste0("IDPROF_", dat$core_latitude, "_", dat$core_longitude)

#algorthms for quantitative pedology

library(aqp)

#generate a soil profile collection object

depths(dat) <- IDPROF ~ depth_min + depth_max

#visualize some profiles

plot(dat[1:10], color='fraction_organic_matter')

# site definition 

site(dat) <- ~ core_latitude + core_longitude

#spatial objects with the rgdal package

library(rgdal)

#define spatial coordinates

coordinates(dat) <- ~ core_longitude + core_latitude

#define reference system

dat@sp@proj4string <- CRS("+proj=longlat +datum=WGS84 +no_defs ")

#visualize the point locations

plot(dat@sp)

#add a base map of national borders

library(maps)

map('world', add=TRUE)

#run splines for estimating the values for the 0-5cm using the GSIF package

library(GSIF)

#fOM is fraction of organic matter

try(fOM <- mpspline(dat, 'fraction_organic_matter', d = t(c(0,5))))#change here for other depths of interest

fOM <- data.frame(
                  x = dat@sp@coords[,1],
                  y = dat@sp@coords[,2],
                  fOM = fOM$var.std[,1])

#check the structure of the new dataset

str(fOM)

#it contains 3 columns, the first 2 are coordinates and the third column is the estimated value for the first 0-5 cm of soil depth, we will generate a new column adding the source

fOM$source = 'Holmquist'

#we save the dataset and continue with the other datasets and other soil depths of interest, 

write.csv(fOM, file='fractionOM_Holmquist_0_5cm.csv')

#end








