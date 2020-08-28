#libraries
 library(ggplot2)
library(dplyr)
library(gbm)
library(raster)
library(openair)
library(parallel)
library(automap)
library(xgboost)
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
###
xy <- cbind(x, y)
names(xy)[16] <- 'label'
xy <- as.matrix(xy)
#X_train = xgb.DMatrix(as.matrix(xy , label = log1p(test$SOC1)))
#X_train = xgb.DMatrix(as.matrix(xy %>% select(-y)))
dtrain <- xgb.DMatrix(as.matrix(x), label = y)
#cross validation
xgb.fit <- xgb.cv(data = dtrain, nrounds = 1000, nthread = 2, nfold = 5, metrics = list("rmse"), objective='reg:linear', early_stopping_rounds = 10)
# plot error vs number trees
ggplot(xgb.fit$evaluation_log) +
  geom_line(aes(iter, train_rmse_mean), color = "red") +
  geom_line(aes(iter, test_rmse_mean), color = "blue") +   theme_minimal(base_size = 20) 
#evaluation
xgb.fit$evaluation_log %>%
  dplyr::summarise(
    ntrees.train = which(train_rmse_mean == min(train_rmse_mean))[1],
    rmse.train   = min(train_rmse_mean),
    ntrees.test  = which(test_rmse_mean == min(test_rmse_mean))[1],
    rmse.test   = min(test_rmse_mean),
  )


# create hyperparameter grid
hyper_grid <- expand.grid(
  eta = c(.01, .05, .1, .3),
  max_depth = c(1, 3, 5, 7),
  min_child_weight = c(0.1, 0.5, 1, 3, 5, 7),
  subsample = c(.65, .8, 1), 
  colsample_bytree = c(.8, .9, 1),
  optimal_trees = 0,               # a place to dump results
  min_RMSE = 0                     # a place to dump results
)



# grid search 
for(i in 1:nrow(hyper_grid)) {
  
  # create parameter list
  params <- list(
    eta = hyper_grid$eta[i],
    max_depth = hyper_grid$max_depth[i],
    min_child_weight = hyper_grid$min_child_weight[i],
    subsample = hyper_grid$subsample[i],
    colsample_bytree = hyper_grid$colsample_bytree[i]
  )
  
  # reproducibility
  #set.seed(123)
  
  # train model
  xgb.tune <- xgb.cv(
    params = params,
    data = dtrain,
   
    nrounds = 5000,
    nfold = 5,
    objective = "reg:linear",  # for regression models
    verbose = 0,               # silent,
    early_stopping_rounds = 10 # stop if no improvement for 10 consecutive trees
  )
  
  # add min training error and trees to grid
  hyper_grid$optimal_trees[i] <- which.min(xgb.tune$evaluation_log$test_rmse_mean)
  hyper_grid$min_RMSE[i] <- min(xgb.tune$evaluation_log$test_rmse_mean)
}

##
# parameter list
params <- list(
  eta = 0.1,
  max_depth = 5,
  min_child_weight = 5,
  subsample = 0.80,
  colsample_bytree = 0.8
)

# train final model
xgb.fit.final <- xgboost(
  params = params,
  data = dtrain,
  
  nrounds = 1576,
  objective = "reg:linear",
  verbose = 0
)


importance_matrix <- xgb.importance(model = xgb.fit.final)

# variable importance plot

xgb.plot.importance(importance_matrix, top_n = 10, measure = "Gain")

# predict values for test data
test$xgb <- predict(xgb.fit.final, xgb.DMatrix(as.matrix(test[1:15])))

# results
caret::RMSE(test$xgb, log1p(test$SOC1))

map_xgb <- predict( xgb.fit.final, xgb.DMatrix(as.matrix(predictor_stack_test[[1:15]])))
ref <- predictor_stack_test[[1]]
values(ref) <- map_xgb
map_xgb <- ref
######

# predict values for test data
ens$xgb <- predict(xgb.fit.final, xgb.DMatrix(as.matrix(ens[1:15])))
#prepare for Kriging of residuals
coordinates(ens) <- ~ coords.x1 + coords.x2
#define CRS
proj4string(ens) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#transform CRS
ens <- spTransform(ens, CRS('+proj=laea +lat_0=45.5 +lon_0=-114.125 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs '))
#get the residual
ens$residual <- ens@data$SOC1 - expm1(ens@data$xgb)
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
writeRaster(map_xgb, file='output/XGB_prediction.tif', overwrite=TRUE )
writeRaster(kr_map, file='output/XGB_OK_residuals.tif', overwrite=TRUE )
#writeRaster(maps, file='output/GBM_variance.tif', overwrite=TRUE )
#

