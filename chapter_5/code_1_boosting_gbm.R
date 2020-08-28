#libraries
library(gbm)
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
#predictors
x <-train[,c(1:15)]
#response
y <- log1p(train$SOC1)
###
###Boosting
hyper_grid <- expand.grid(
  shrinkage = c(.01, .05, .1 ),
  interaction.depth = c(3,4,5),
  n.minobsinnode = c(1,2,3),
  bag.fraction = c(.65, .8, 1), 
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)

# total number of combinations
nrow(hyper_grid)
## [1] 81
#empty data frame to store results
res <- data.frame()
#empty raster stack to see the variance of predictions
maps<- stack()
#loop  searching optimal parameter value
for(i in 1:nrow(hyper_grid)) {
 
  # reproducibility
  
# train GBM model
gbm.tune <- gbm(
  formula = y ~ .,
  distribution = "gaussian",
  data = cbind(x, y),
  n.trees = 1000,
    interaction.depth = hyper_grid$interaction.depth[i],
    shrinkage = hyper_grid$shrinkage[i],
    n.minobsinnode = hyper_grid$n.minobsinnode[i],
    bag.fraction = hyper_grid$bag.fraction[i],
    train.fraction = .5,
	 cv.folds = 5,
    n.cores = NULL, # will use all cores by default
    verbose = FALSE
  )  
print('error and parameters')
print(sqrt(min(gbm.tune$cv.error)))
print(hyper_grid$interaction.depth[i])
print(hyper_grid$shrinkage[i])
print(hyper_grid$n.minobsinnode[i])
print(hyper_grid$bag.fraction[i])
pred <- predict(predictor_stack_test, gbm.tune)
maps <- stack(maps, pred)
RES <- data.frame(err=sqrt(min(gbm.tune$cv.error)), depth=hyper_grid$interaction.depth[i], shrin=hyper_grid$shrinkage[i], minobs=hyper_grid$n.minobsinnode[i], bag=hyper_grid$bag.fraction[i])
res <- rbind(res, RES)
hyper_grid$optimal_trees[i] <- which.min(gbm.tune$valid.error)
hyper_grid$min_RMSE[i] <- sqrt(min(gbm.tune$valid.error))
}
#get the best fit based on error
print(bestfit <- res[res$err==min(res$err),])
#use the mdoel with all data and best parameters
# train GBM model
#random number generator

gbm.fit.final <- gbm(
  formula = y ~ .,
  distribution = "gaussian",
  data = cbind(x, y),
  n.trees = 1000,
  interaction.depth = bestfit[,2],
  shrinkage = bestfit[,3],
  n.minobsinnode = bestfit[,4],
  bag.fraction = bestfit[,5], 
  train.fraction = 1,
  cv.folds=5,
  n.cores = NULL, # will use all cores by default
  verbose = FALSE
)

#check variable importance by permitation
par(mar = c(5, 16, 1, 1))
summary(
  gbm.fit.final, 
  cBars = 10,
  method = permutation.test.gbm, # also can use permutation.test.gbm
  las = 2
  )

#build a map from model with optimal parameters
map_gbm <- predict(predictor_stack_test , gbm.fit.final )
#get the prediction variance (surrogate of uncertainty) of the map
var <- calc(maps, sd)*2
#convert to log in the testing set
test$log1p_SOC1 <- log1p(test$SOC1)
#make a prediction
test$gbm <- predict(gbm.fit.final, test )
#calculate the error
Metrics::rmse(log1p(test$SOC1), test$gbm)
#make a quantile plot
conditionalQuantile(test, obs = "log1p_SOC1", mod = "gbm")
#test$gbm_sc <- rescale(test$gbm, to = c(min(test$log1p_SOC1), max(test$log1p_SOC1)), from = range(test$log1p_SOC1, na.rm = TRUE, finite = TRUE))#scale data to min and max
#conditionalQuantile(test, obs = "log1p_SOC1", mod = "gbm_sc")#plot again

#predict to the ensemble dataset
ens$gbm <- predict(gbm.fit.final, ens )
#prepare for Kriging of residuals
coordinates(ens) <- ~ coords.x1 + coords.x2
#define CRS
proj4string(ens) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#transform CRS
ens <- spTransform(ens, CRS('+proj=laea +lat_0=45.5 +lon_0=-114.125 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs '))
#get the residual
ens$residual <- ens@data$SOC1 - expm1(ens@data$gbm)
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
    # Create SpatialPixelsDataFrame from mergeParallelX
mergeParallelX <- SpatialPixelsDataFrame(points = mergeParallelX, data = mergeParallelX@data)
kr_map <- raster(mergeParallelX)[[1]]
#save the maps
writeRaster(map_gbm, file='output/GBM_prediction.tif', overwrite=TRUE )
writeRaster(kr_map, file='output/GBM_OK_residuals.tif', overwrite=TRUE )
writeRaster(maps, file='output/GBM_variance.tif', overwrite=TRUE )
#



