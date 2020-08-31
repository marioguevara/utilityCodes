#libraries
library(dplyr)
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
##model stacking and others using autoML in h2o
#data frame with predictors and response (x, and y respectively)
xy <- cbind(x,y)
#name of response
y <- "y"
#names of predictors
x <- setdiff(names(xy), y)
#h2o format
write.csv(xy, file='h2o_dataset.csv', row.names=FALSE)
#get library of interest
library(h2o)
#initialize the parallel environment
h2o.init()
#import the h2o data format
xy <- h2o.importFile("h2o_dataset.csv")
#test 50 models automatically
auto_ml <- h2o.automl(
  x =x, y = y, training_frame = xy, nfolds = 5, 
  max_models = 30, #max_runtime_secs = 60 * 120,
  keep_cross_validation_predictions = TRUE, sort_metric = "RMSE", seed = 123,
  stopping_rounds = 50, stopping_metric = "RMSE", stopping_tolerance = 0
)
#check best results
auto_ml@leaderboard %>% 
  as.data.frame() %>%
  dplyr::select(model_id, rmse) %>%
  dplyr::slice(1:25)	
###


#auto_ml_reduce <- h2o.automl(
# x =x, y = y, training_frame = xy, nfolds = 5, 
# include_algos = c("GLM", "DeepLearning", "DRF", 'XRF'),
# keep_cross_validation_predictions = TRUE, sort_metric = "RMSE", seed = 123,
# stopping_rounds = 50, stopping_metric = "RMSE", stopping_tolerance = 0
#)

#check best results
#auto_ml_reduce@leaderboard %>% 
 # as.data.frame() %>%
 #dplyr::select(model_id, rmse) %>%
 #dplyr::slice(1:25)	
###



write.csv(test, file='h2o_testing_dataset.csv', row.names=FALSE)
test  <- h2o.importFile("h2o_testing_dataset.csv")
test$aml <- h2o.predict(auto_ml , test) 
#test$aml_reduce <- h2o.predict(auto_ml_reduce , test) 
test$aml_leader <- h2o.predict(auto_ml@leader, test)



write.csv(ens, file='h2o_ens_dataset.csv', row.names=FALSE)
ens  <- h2o.importFile("h2o_ens_dataset.csv")
ens$aml <- h2o.predict(auto_ml , test) 
#ens$aml_reduce <- h2o.predict(auto_ml_reduce , test) 
ens$aml_leader <- h2o.predict(auto_ml@leader, test)


data <- as.data.frame(predictor_stack_test[[1:15]],xy =TRUE)
data <- na.omit(data)
write.csv(data, file='h2o_predictors_dataset.csv', row.names=FALSE)
data2  <- h2o.importFile("h2o_predictors_dataset.csv")
#data$aml <- h2o.predict(auto_ml , data) 
data2$aml_leader <- h2o.predict(auto_ml@leader , data2) 

data <- as.data.frame(data2)


coordinates(data) <- ~ x+y

 map_stacking <- rasterize(data, raster(readRDS('input/predictor_stack_test_sgd.rds')), 'aml_leader')

#ensemble dataset
ens <- DATA[[3]]
ens$aml_leader <- extract(map_stacking, ens[c('coords.x1', 'coords.x2')])


#prepare for Kriging of residuals
coordinates(ens) <- ~ coords.x1 + coords.x2
#define CRS
proj4string(ens) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#transform CRS
ens <- spTransform(ens, CRS('+proj=laea +lat_0=45.5 +lon_0=-114.125 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs '))
#get the residual
ens$residual <- ens@data$SOC1 - expm1(ens@data$aml_leader)
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
writeRaster(map_rethinking, file='output/STACK_AUTOML_prediction.tif', overwrite=TRUE )
writeRaster(kr_map, file='output/STACK_AUTOML_OK_residuals.tif', overwrite=TRUE )
#



