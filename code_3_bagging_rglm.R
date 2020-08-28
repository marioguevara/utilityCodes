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
#predictors (terrain parameter columns 1:15 that we have at 10m and 1km grids)
x <-train[,c(1:15)]
#response (SOIL CARBON, NOTE THAT THERE ARE SIX DEPTHS)
y <- log1p(train$SOC1)
###
###randomGLM

######
library(randomGLM)
# different kinds of accuracy measures for a continuous outcome
accuracyCor=function(y.predicted,y) cor(y.predicted, y, use="p")
accuracyMSE=function(y.predicted,y) mean((y.predicted-y)^2,na.rm=TRUE)
accuracyMedAbsDev=function(y.predicted,y) {median( abs(y.predicted-y),na.rm=TRUE)}
accuracyM = accuracyCor
# choose nFeaturesInBag and nCandidateCovariates at the same time
nFeatureInBagVector=ceiling(c(0.1, 0.2, 0.4, 0.6, 0.8, 1)*dim(x)[2])
nCandidateCovariatesVector=c(2:15)
# define vector that saves prediction accuracies
acc2=matrix(NA,length(nFeatureInBagVector),length(nCandidateCovariatesVector))
rownames(acc2) = paste("feature", nFeatureInBagVector, sep="")
colnames(acc2) = paste("cov", nCandidateCovariatesVector, sep="")
GLMmaps <-stack()
# loop over nFeaturesInBag and nCandidateCovariates values, calculate individual accuracy
for (i in 1:length(nFeatureInBagVector))
{
  cat("step", i, "out of ", length(nFeatureInBagVector), "entries from nFeatureInBagVector\n")
  for (j in 1:length(nCandidateCovariatesVector))
  {    
    cat("step", j, "out of ", length(nCandidateCovariatesVector), "entries from nCandidateCovariatesVector\n")
    RGLMtmp = randomGLM(x, y, 
	classify=FALSE, 
	nFeaturesInBag = nFeatureInBagVector[i],
	nCandidateCovariates = nCandidateCovariatesVector[j],
	nBags=50,
	keepModels=TRUE)
    predicted = RGLMtmp$predictedOOB
    acc2[i, j] = accuracyM(predicted, y)
GLMtmp <-  predict(predictor_stack_test[[names(x)]], RGLMtmp)
GLMmaps <- stack(GLMmaps, GLMtmp)
    rm(RGLMtmp, GLMtmp, predicted)
  }
}
round(acc2,2)

RGLMtmp = randomGLM(x, y, 
	classify=FALSE, 
	nFeaturesInBag = 15,
	nCandidateCovariates = 10,
	nBags=50,
	keepModels=TRUE)

test$randomGLM <-  predict(RGLMtmp, test[names(x)])
caret::RMSE(test$randomGLM, log1p(test$SOC1))
map_rglm <-  predict(predictor_stack_test[[names(x)]], RGLMtmp)

#var_rglm <- calc(GLMmaps, sd)^2
#map_rglm <- calc(GLMmaps, mean)

# predict values for test data
ens$rglm <- predict(RGLMtmp, ens[1:15])
#prepare for Kriging of residuals

#prepare for Kriging of residuals
coordinates(ens) <- ~ coords.x1 + coords.x2
#define CRS
proj4string(ens) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
#transform CRS
ens <- spTransform(ens, CRS('+proj=laea +lat_0=45.5 +lon_0=-114.125 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs '))
#get the residual
ens$residual <- ens@data$SOC1 - expm1(ens@data$rglm)
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
writeRaster(map_rglm, file='output/RGLM_prediction.tif', overwrite=TRUE )
writeRaster(kr_map, file='output/RGLM_OK_residuals.tif', overwrite=TRUE )
#

