
library(raster)
library(caret)
library(openair)

srdb <- read.csv('https://raw.githubusercontent.com/bpbond/srdb/master/srdb-data.csv')

#selecciona variables
  RS <-  srdb[c('Longitude' , 'Latitude', 'Rs_annual', 'Entry_date')]
  
  RS$yr <- 1


for (i in 1:nrow(RS)){RS$yr[i] <- unlist(strsplit(RS$Entry_date[i],'-'))[1]}

yrs <- as.numeric(levels(as.factor(RS$yr)))

lis <- list.files(pattern='45.tif')

RSsm <- data.frame()

for (i in 1:length(yrs)){

RSi <- RS[RS$yr == yrs[i],]

RSi <- RSi[is.na(RSi$Rs_annual) == FALSE,]

RSi <- RSi[is.na(RSi$Longitude) == FALSE,]

ras <- stack(lis[grep( as.character(yrs[i]), lis)])

names(ras) <- 'sm_mean'

coordinates(RSi) <- ~ Longitude + Latitude

RSi$sm <- extract(ras, RSi)

RSi$soc <- extract(raster("/home/mario/Downloads/GSOCmap1.5.0.tif"), RSi)
RSi$st <- extract(raster("/home/mario/Downloads/SBIO1_Annual_Mean_Temperature_0_5cm.tif"), RSi)

df <- as.data.frame(RSi)

RSsm <- rbind(RSsm, df)

}


RSsm <- RSsm[is.na(RSsm$soc) == FALSE, ]


RSsm <- RSsm[RSsm$Rs_annual > 0,]
RSsm <- RSsm[is.na(RSsm$sm_mean) == FALSE,]


#library(dplyr)
#df=RSsm[-4] %>% group_by(yr) %>% summarise_each(funs(mean))

df=RSsm
df$yr <- log1p(as.numeric(df$yr))
df <- data.frame(df)
df$logRs <- log1p(df$Rs_annual)
df$logSm <- log1p(df$sm_mean)
df$logSoc <- log1p(df$soc)
df$logSt <- log1p(df$st)

	#library(caret)

	# define the control using a random forest selection function
	#control <- rfeControl(functions=rfFuncs, method="repeatedcv", repeats=10)
	# run the RFE algorithm
	
	#library(doMC)
        #registerDoMC(10)

	#results <- rfe(df[,c(5, 9, 10)], df[,8], sizes=c(1:3), rfeControl=control)
	# summarize the results
	#print(results)
	# list the chosen features
	#predictors(results)
	# plot the results
	#plot(results, type=c("g", "o"))
	
	df_na <- na.omit(df)

	 fit <- train(logRs ~ logSm + logSoc + logSt + yr, data=df_na, method='ranger', trControl = trainControl(method='repeatedcv',savePredictions = TRUE, repeats=5))
res <- data.frame(mod=fit$pred$pred, obs=fit$pred$obs)
 conditionalQuantile(res, obs = "obs", mod = "mod")
 openair::modStats(res, obs='obs', mod='mod')


#fit_orig <- train(Rs_annual ~ sm_mean + soc + yr, data=df, method='ranger', trControl = trainControl(method='repeatedcv',savePredictions = TRUE, repeats=5))
#res_orig <-  data.frame(mod=fit_orig$pred$pred, obs=fit_orig$pred$obs)
#conditionalQuantile(res_orig, obs = "obs", mod = "mod")
# openair::modStats(res_orig, obs='obs', mod='mod')


library(ranger)

 ran <- ranger(logRs ~ logSm + logSoc + logSt + yr, data=df_na, importance = "impurity", num.threads=10, keep.inbag=TRUE )
 
 #ran <- ranger(Rs_annual ~ sm_mean + soc + yr, data=df, importance = "impurity", num.threads=10, keep.inbag=TRUE )


library(ggplot2)
 library(dplyr)

imps <- data.frame(var = ran$forest$independent.variable.names,
                   imps = ran$variable.importance/max(ran$variable.importance))
imps %>% 
  ggplot(aes(imps, x = reorder(var, imps))) +
  geom_point(size = 3, colour = "#ff6767") +
  coord_flip() +
  labs(x = "Predictors", y = "Importance scores") +
  theme_bw(18)


 SM <- stack(lis)
 soc <- raster("/home/mario/Downloads/GSOCmap1.5.0.tif")
 st <-raster("/home/mario/Downloads/SBIO1_Annual_Mean_Temperature_0_5cm.tif")
 soc <- projectRaster(soc, SM)
  st <- projectRaster(st, SM)
 PREDS <- stack()

NA2mean <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
for (i in 1:dim(SM)[3]){
 sm <- SM[[i]]
 
pr <- stack(soc, sm, st)
 names(pr)[1] <- 'soc'
 names(pr)[2] <- 'sm'
 names(pr)[3] <- 'st'
 
pr <- as(pr, 'SpatialPixelsDataFrame')
pr$yr <- 1990+i
#pr$soc <- over(pr, as(raster("/home/mario/Downloads/GSOCmap1.5.0.tif"), 'SpatialPixelsDataFrame'))

 pr$logSoc <- log1p(pr$soc)
 pr$logSm <- log1p(pr$sm)
  pr$logSt <- log1p(pr$st)
 
 pr@data[] <- lapply(pr@data, NA2mean)
 #pr <- spatialEco::sp.na.omit(pr)
 #pr@data <- na.omit(pr@data)
 
 pr$logRs_predicted <- predict(fit , pr)
 #spplot(pr['logRs_predicted'])


 
 pr$Rs_predicted_back <-  expm1(pr$logRs_predicted)
 
 pr <- raster(pr['Rs_predicted_back']) 

 
 PREDS <- stack(PREDS, pr)
 
 print(i)
 }
 
 library(miscset)

 library(greenbrown)

 trendsRS <-  TrendRaster(PREDS, start=c(1991, 1), freq=1)
 
 setwd("/home/mario/Downloads/curso/sm_kknn_eco_swc_terrain_15km.")

 s <- stack(list.files(pattern='45.tif'))
 
 library(greenbrown)	

 trendsSM <-  TrendRaster(s, start=c(1991, 1), freq=1)

meanRS <- calc(s, median)


library(raster)
library(rgdal)

dat <- read.csv("/home/mario/Downloads/XY_PRACTICAS.csv")

dat <- dat[is.na(dat$Long) == FALSE,]

dat$Lat <- as.numeric(dat$Lat)
	


dat$Long <- as.numeric(dat$Long)
dat <- dat[is.na(dat$Long) == FALSE,]


coordinates(dat) <- ~ Long + Lat

#meanRS <- calc(PREDS, median)


ref <- as(meanRS, 'SpatialPixelsDataFrame')

dat$meanRS <- extract(meanRS, dat)
plot(density(dat$meanRS), lwd=2, col='red')
num <-as.numeric()
for (i in 1:100){
xy <- spsample(ref,n=105,"random")
xy$meanRS <- extract(meanRS, xy)
num[i] <- median(xy$meanRS)
lines(density(xy$meanRS), col='gray', lty=2, lwd=0.5)
}
boxplot(num, dat$meanRS)
abline(v=median(num), col='gray')
abline(v=median(dat$meanRS), col='red')

confint(num)
confint(dat$meanRS)

dat$trmeanrs <- extract(trendsRS[[2]], dat)
plot(density(dat$trmeanrs), lwd=2, col='red')
num <-as.numeric()
for (i in 1:100){
xy <- spsample(ref,n=105,"random")
xy$trmeanrs <- extract(trendsRS[[2]], xy)
num[i] <- median(xy$trmeanrs)
lines(density(xy$trmeanrs), col='gray', lty=2, lwd=0.5)
}
boxplot(num, dat$trmeanrs)
abline(v=median(num), col='gray')
abline(v=median(dat$trmeanrs), col='red')


dat$meansoc <- extract(soc, dat)
plot(density(dat$meansoc), lwd=2, col='red')
num <-as.numeric()
for (i in 1:100){
xy <- spsample(ref,n=105,"random")
xy$meansoc <- extract(soc, xy)
num[i] <- median(xy$meansoc)
lines(density(xy$meansoc), col='gray', lty=2, lwd=0.5)
}
boxplot(num, dat$meansoc)
abline(v=median(num), col='gray')
abline(v=median(dat$meansoc), col='red')


dat$meansm <- extract(meanSM, dat)
plot(density(dat$meansm), lwd=2, col='red')
num <-as.numeric()
for (i in 1:100){
xy <- spsample(ref,n=105,"random")
xy$meansm <- extract(meanSM, xy)
num[i] <- median(xy$meansm)
lines(density(xy$meansm), col='gray', lty=2, lwd=0.5)
}
boxplot(num, dat$meansm)
abline(v=median(num), col='gray')
abline(v=median(dat$meansm), col='red')



dat$trmeansm <- extract(trendsSM[[2]], dat)
plot(density(dat$trmeansm), lwd=2, col='red')
num <-as.numeric()
for (i in 1:1000){
xy <- spsample(ref,n=105,"random")
xy$trmeansm <- extract(trendsSM[[2]], xy)
num[i] <- median(xy$trmeansm)
lines(density(xy$trmeansm), col='gray', lty=2, lwd=0.5)
}
boxplot(num, dat$trmeansm)
abline(v=median(num), col='gray')
abline(v=median(dat$trmeansm), col='red')


dat$st <- extract( raster("/home/mario/Downloads/SBIO1_Annual_Mean_Temperature_0_5cm.tif"), dat)
plot(density(dat$st), lwd=2, col='red', xlim=c(-10, 10))
num <-as.numeric()
for (i in 1:100){
xy <- spsample(ref,n=105,"random")
xy$st <- extract(raster("/home/mario/Downloads/SBIO1_Annual_Mean_Temperature_0_5cm.tif"), xy)
num[i] <-median(na.omit(xy$st))
lines(density(na.omit(xy$st)), col='gray', lty=2, lwd=0.5)
}

abline(v=median(num), col='gray')
abline(v=median(dat$st), col='red')
 
