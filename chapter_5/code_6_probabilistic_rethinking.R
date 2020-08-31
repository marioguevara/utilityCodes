#libraries
library(rethinking)
library(raster)
library(openair)
library(parallel)
library(automap)
#prediction factor dataset
predictor_stack_test <- stack(readRDS('input/predictor_stack_test_sgd.rds'))
#training dataset (list)
DATA <- readRDS('input/train_test_validate_0827mg.rds')
#training data
train <- DATA[[1]]
#testing data
test <- DATA[[2]]
#ensemble dataset
ens <- DATA[[3]]
#predictors (terrain parameter columns 1:15 that we have at 10m and 1km grids)
x <-train[,c(1:15)]
#response (SOIL CARBON, NOTE THAT THERE ARE SIX DEPTHS)
y <- log1p(train$SOC1)
library(rethinking)
m1 <- glimmer(y~., cbind(x,y))
     ms1 <- map2stan( m1$f , data=m1$d , cores=5)
m2 <- glimmer(y~., cbind(x[3],y))
     ms2 <- map2stan( m2$f , data=m2$d , cores=5)
m3 <- glimmer(y~., cbind(x[7],y))
     ms3 <- map2stan( m3$f , data=m3$d, cores=5 )
m4 <- glimmer(y~., cbind(x[13],y))
     ms4 <- map2stan( m4$f , data=m4$d, cores=5 )
m5 <- glimmer(y~., cbind(x[c(3, 7, 13)],y))
     ms5 <- map2stan( m5$f , data=m5$d , cores=5)
ensample_test <- ensemble (ms1, ms2, ms3, ms4, ms5 , data=test[names(x)])
ensample_ens <- ensemble (ms1, ms2, ms3, ms4, ms5 , data=ens[names(x)])
test$rethinking <- apply(ensample_test$link , 2 , mean )
ens$rethinking <- apply(ensample_ens$link , 2 , mean )

(predictor_stack_test <- stack(readRDS('input/predictor_stack_test_sgd.rds')))

library(GSIF)
tl <- getSpatialTiles(readRDS('input/predictor_stack_test_sgd.rds'), block.x=2.5)

results <- data.frame()

for (i in 1:length(tl)){
data <- as.data.frame(crop(predictor_stack_test[[1:15]], tl[i]),xy =TRUE)
data <- na.omit(data)
#rm(predictor_stack_test)
ensamble_a <- ensemble (ms1, ms2, ms3, ms4, ms5 , data=data[names(x)])
data$rethinking <- apply( ensamble_a$link , 2 , mean )
results <- rbind(results, data)
print(i)
}

coordinates(results) <- ~ x+y
map_rethinking<- rasterize(data, raster(readRDS('input/predictor_stack_test_sgd.rds')), 'rethinking')


ens$rethinking <- extract(map_rethinking, ens[c('coords.x1', 'coords.x2')])


#prepare for Kriging of residuals
coordinates(ens) <- ~ coords.x1 + coords.x2
#define CRS
proj4string(ens) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#transform CRS
ens <- spTransform(ens, CRS('+proj=laea +lat_0=45.5 +lon_0=-114.125 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs '))
#get the residual
ens$residual <- ens@data$SOC1 - expm1(ens@data$rethinking)
#extract the variable
data <- ens['residual']
#residual response (generic)
names(data) <- 'response'
#au <- autoKrige(response~1., spSample)
#set n cores
no_cores <- 7
#get a metric CRS for reference dataset
grid <- projectRaster(predictor_stack_test, crs=CRS('+proj=laea +lat_0=45.5 +lon_0=-114.125 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs '))
#convert to SPDF
grid <- as(grid, 'SpatialPixelsDataFrame')
#split to krige in parallel
parts <- split(x = 1:length(grid), f = 1:no_cores)
#estimate variogram parameters

 
  ok.vgm <- autofitVariogram(formula = response ~ 1, input_data = data,
                             model = c("Sph", "Exp", "Gau", "Ste"))
#plot variogram
plot(ok.vgm)
#make cluster
    cl <- makeCluster(no_cores)
    clusterExport(cl = cl, varlist = c("data", "ok.vgm",  "grid", "parts"), envir =  environment())
  clusterEvalQ(cl = cl, expr = c(library('sp'), library('gstat')))
#run Kriging prediction and stop cluster
  parallelX <- parLapply(cl = cl, X = 1:no_cores, fun = function(x) krige(formula = data$response ~ 1, locations = data, newdata = grid[parts[[x]],], model = ok.vgm$var_model, nmax=100))
  stopCluster(cl)
#get the maps  
library(maptools)
  # Merge all the predictions
  mergeParallelX <- maptools::spRbind(parallelX[[1]], parallelX[[2]])
  for(i in 3:no_cores){
    mergeParallelX <- maptools::spRbind(mergeParallelX, parallelX[[i]])
  }
#Create SpatialPixelsDataFrame from mergeParallelX
mergeParallelX <- SpatialPixelsDataFrame(points = mergeParallelX, data = mergeParallelX@data)
kr_map <- raster(mergeParallelX)[[1]]
#save the maps
writeRaster(map_rethinking, file='output/RETHINKING_prediction.tif', overwrite=TRUE )
writeRaster(kr_map, file='output/RETHINKING_OK_residuals.tif', overwrite=TRUE )
#






