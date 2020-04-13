#GSSmap, approach for no-data countries
#prepared by MG-Pillar 4 Latin America GSP
#Software R
#external inputs: country-specific polygon maps (soils, land use, vegetation, geology)
#output: prediction map 
#04132020
#libraries required
library(raster)
library(elevatr)
library(clhs)	
library(RStoolbox)
library(automap)
library(mlbench)
library(caret)
library(spatialEco)
#import the soil map
soil <- shapefile ("Soils.shp")
#get elevation data (9=~175m grids)
elev <- get_elev_raster(soil, prj = projection(soil), z = 9, clip = "bbox")
#mask for the study area
elev <- mask(elev, soil)
#name the data within the raster
names(elev) <- 'elev'
#derive slope and aspect
slp_asp <- terrain(elev, opt=c('slope', 'aspect'), unit='degrees')
#TPI for different neighborhood size:
     tpiw <- function(x, w=5) {
             m <- matrix(1/(w^2-1), nc=w, nr=w)
             m[ceiling(0.5 * length(m))] <- 0
             f <- focal(x, m)
             x - f
     }
#derive topographic position index (TPI)
tpi <- tpiw(elev, w=5)
#define a name for the data in the tpi object
names(tpi) <- 'tpi'
#define field of interest in soil map
soil@data$DESCRIP <- as.factor(soil@data$DESCRIP)
#rasterize the field on interest to the grid reference of elevation data 
soil.r <- rasterize(x = soil, y = elev, field = "DESCRIP")
#function to convert to dummy variables each soil polygon
dummyRaster <- function(rast){
 rast <- as.factor(rast)
 result <- list()
 for(i in 1:length(levels(rast)[[1]][[1]])){
 result[[i]] <- rast == levels(rast)[[1]][[1]][i]
 names(result[[i]]) <- paste0(names(rast),
 levels(rast)[[1]][[1]][i])
 }
 return(stack(result))
 }
#apply the function to the rasterized map 
soilmap_dummy <- dummyRaster(soil.r$layer)
#re-name the soil map converted to dummy variables
names(soilmap_dummy) <- levels(as.factor(soil$DESCRIP))
#import the vegetation map 
veg <- shapefile('Vegetation09.shp')
#define as factor the column of interest
veg$Descript <- as.factor(veg$Descript)
#rasterize the se	lected field
veg.r <- rasterize(x = veg, y = elev, field = "Descript")
#converto to dummy
vegmap_dummy <- dummyRaster(veg.r$layer)
#re-name the dummy map
names(vegmap_dummy) <- levels(as.factor(veg$Descript))
#visuzlize dummy variables
plot(vegmap_dummy)
#import the geology map
geo <- shapefile("Geology.shp")
#define as factor the column of interest
geo$DESCRIPT <- as.factor(geo$DESCRIPT)
#rasterize the selected field
geo.r <- rasterize(x = geo, y = elev, field = "DESCRIPT")
#convert to dummy variables
geomap_dummy <- dummyRaster(geo.r$layer)
#re-name the map
names(geomap_dummy) <- levels(as.factor(geo$DESCRIPT))
#visualize the map
plot(geomap_dummy)
#repeat with land use map
lu <- shapefile('landcover2009.shp')
#define the field of interest
lu$Descrip <- as.factor(lu$Descrip)
#rasterize the selected field
lu.r <- rasterize(x = lu, y = elev, field = "Descrip")
#convert it to dummy variables
lumap_dummy <- dummyRaster(lu.r$layer)
#re-name each dummy map
names(lumap_dummy) <- levels(as.factor(lu$Descrip))
#visualize the result
plot(lumap_dummy)
#stack all covariates in the same R object (raster stack)
covs <- stack(elev, tpi, slp_asp, soilmap_dummy, vegmap_dummy, lumap_dummy, geomap_dummy)
#include more covariates (download bioclimate layers)
bio <- raster::getData('worldclim', var='bio',  res=0.5, lon=-60.9, lat=13.9)
#reproject to the geographical reference of the other covariates
bio.r <- projectRaster(bio, crs=projection(covs))    
#crop bioclimatic layers to area of interest
bio.r <- crop(bio.r, covs)
#conver them to spatial points to downscale
x <- as(bio.r, 'SpatialPointsDataFrame')
#reference grid for downscaling bioclimate layers using Kriging
ref <- as(covs[[1]], 'SpatialPixelsDataFrame')
#loop to automated Kriging applied to bioclimate layers 
s <- stack()
for (i in 1:19){
dat <- x[i]
names(dat) <- 'res'
a <- autoKrige(res~1., dat, ref)
ras <- raster(a$krige_output['var1.pred'])
s<- stack(ras, s)
print(i)
}
#define the name for the downscaled bioclimatic layers
names(s) <- paste0('bio_' , 1:19)
#add the downscaled bioclimate layers to the area of interest
covs <- stack(covs, s)
#visualize the covariates rapidly 
for (i in 1:101){plot(covs[[i]])}
#apply a feature extraction strategy (PCA) to the covariates
covs_pca <- rasterPCA(covs)
#converto to points the covariate file
covs_pts <- rasterToPoints(covs_pca$map, spatial=TRUE)
#define a seed for reprudicing results
set.seed(666)
#run a sampling strategy based on a latin hypercube
s.clhs <- clhs(covs_pts, size =1000, progress = FALSE, iter = 100, simple = FALSE)
#get the sampled points
pts <- s.clhs$sampled_data
#derive a regression matrix
train <- data.frame(data.frame(extract(covs, pts)))#, pts@coords)
#identify only zero columns
ind <- sapply(train, function(x) sum(x==0)) != nrow(train)
#remove only zero columns in regression matrix
train <- train[,ind]
#train a model applied to feature extraction
set.seed(7)
# define the control using a random forest selection function
control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# run the RFE algorithm
results <- rfe(train[,-14], train[,14], sizes=c(1:8), rfeControl=control)
# summarize the results
print(results)
# list the chosen features
predictors(results)
# plot the results
plot(results, type=c("g", "o"))
#perform a model prediction
 prediction <- predict(covs, results)
#standardize the map from 0 (no SAL) to 1 (SAL)
prediction_trans <- raster.transformation(prediction, trans = "stretch", smin = 0, smax = 1)
#get country limit
lim <- getData('GADM', country='LC', level=1)
#project it to the national reference system
lim <- spTransform(lim, CRS(projection(soil)))
#plot the prediction map
plot(prediction_trans, col=colorRampPalette(c("gray","gold", "red"))(5), main='Salinity affected soils in Saint Lucia')
#add the limit of the country
plot(lim, add=TRUE)
#save prediction in a raster file
writeRaster(prediction_trans, file='gss_st_lucia_0413mg.tif', overwrite=TRUE)
#Next: Improve covariates, add more terrain parameters and include your remote sensing imagery, repeat and compare.
#Next: Repeat using multiple datasets for assessing sensitivity and uncertainty of models to data variations
#Next: Repeat with other soil properties
