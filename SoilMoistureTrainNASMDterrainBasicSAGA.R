#Prepared By Mario Guevara for the Soil Carbon Study Tour at UMBS, thanks to US Forest Service
#University of Delaware
#



#Directorio de trabajo y librerias


install.packages(c('dichromat','DAAG', 'Metrics', 'car', 'leaps',
 'RSAGA', 'GSIF', 'plotKML',  'FactoMineR'))



rm(list=ls())

###benchmark tool

CronometroON<- function(){
      tic<-proc.time()[3]
      assign(".tic", tic, envir=baseenv())
      invisible(tic)
      }
CronometroOFF<- function(){
      tic <- get(".tic", envir=baseenv())
      toc<-proc.time()[3]-tic
      hrs<-as.integer(toc/3600)
      minu<- as.integer(((toc/3600)-hrs)*60)
      seg<- ((((((toc/3600)-hrs)*60)))-minu)*60
      time<-paste(as.character(hrs),"hrs ",as.character(minu),"min ",as.character(round(seg,digit=2)),"seg",sep="")
      return(time)
      }

CronometroON()

library(dichromat)
library(DAAG)
library(Metrics)
library(car)
library(leaps)
library(rgdal)
library(raster)
library(RSAGA)
library(GSIF)
library(plotKML)
library(automap)
library(FactoMineR)
library(rasterVis)
library(RStoolbox)

#setwd('D:\\perf_cos_30cm')
s=stack('basicTerrain5kmMXUSsagaGis.tif')
#Importa shapefile con datos de CO

require("rgdal")
require("rgeos")
require("dplyr")
us=readRDS('USA_adm2.rds')
usa <- gUnaryUnion(us, id = us@data$NAME_1)
s=crop(s, usa)
s=mask(s,usa)

s=stack(list.files()[1])
r.scale=s[[1]]

r.min = cellStats(r.scale, "min")
r.max = cellStats(r.scale, "max")
r.scale <- ((r.scale - r.min) / (r.max - r.min) - 0.5 ) * 2

plot(r.scale*100, col=colorRampPalette(c("red", "white", "blue"))(255))
plot(s[[2]], zlim=c(0.05, 1), col='gray80', add=T, legend=FALSE)
plot(usa, add=T)

zoom(r.scale*100, col=colorRampPalette(c("red", "white", "blue"))(255))
plot(s[[2]], zlim=c(0.05, 1), col='gray80', add=T, legend=FALSE)
plot(usa, add=T)



set.seed(25)
rpc <- rasterPCA(r.scale)
rpc
summary(rpc$model)
ggRGB(rpc$map,1,2,3, stretch="lin", q=0)

wave <- data.frame(x = cov$x, y = cov$y)
     
p <- ggplot(wave, aes(x, y))+
		theme(text = element_text(size=20),
        	axis.text.x = element_text(angle=90, vjust=1)) 


     p + ggRGB(rpc$map,1,2,3, stretch="lin", q=0, ggLayer = TRUE)


shape=readRDS( "global5cm.rds"  )
coordinates(shape)=~Longitude+Latitude
#shape <- readOGR(dsn = ".", layer = "NASMD0_5cmSoilMoisture1km_Terrain5kmBasicSAGA")
e=extract(s, shape)
dat=as.data.frame(shape)
dat=cbind(dat, e)
dat=data.frame(x=dat[,7], y=dat[,8],sm=dat[,6], scale(dat[,10:24]), Y=dat$Y)
all_data=dat[complete.cases(dat),]
cov=as.data.frame(scale(r.scale), xy=TRUE)
cov=cov[complete.cases(cov),]
str(all_data)
str(cov)

positions <- sample(nrow(all_data),size=floor((nrow(all_data)/4)*3))
training<- all_data[positions,]
testing<- all_data[-positions,]

library(e1071)
svm_fit<-svm(sm~.,data=training[-c(1,2,19)])
#svm_predictions<-predict(svm_fit,newdata=testing)
error<-sqrt((sum((testing$y-svm_predictions)^2))/nrow(testing))
error




library(foreach)
library(doSNOW)
cl<-makeCluster(4) #change the 4 to your number of CPU cores
registerDoSNOW(cl)  
num_splits<-4
split_testing<-sort(rank(1:nrow(testing))%%4)

svm_predictions<-foreach(i=unique(split_testing),
.combine=c,.packages=c("e1071")) %dopar% {
as.numeric(predict(svm_fit,newdata=testing[split_testing==i,]))
}
stopCluster(c1)
#


#Importa DEM

dem=raster('AsterDEM_HN30m.img')
dem
dem30ag=aggregate(dem, 30, mean)
dem30ag
dem30agsd=aggregate(dem, 30, sd)

#writeRaster(dem30ag, file='dem30ag.tif', overwrite=TRUE)

#SAGA GIS environment
#env <- rsaga.env(path = "/usr/bin")

#SAGA GIS digital terrain parameters 

#rsaga.slope(in.dem = "dem30ag.tif", out.slope = "slope.sgrd", env = env)
#rsaga.sgrd.to.esri("slope", "slope", format = "ascii", env=env)
#slp=raster("slope.asc")

#rsaga.aspect(in.dem = "dem30ag.tif", out.aspect = "aspect.sgrd", env = env)
#rsaga.sgrd.to.esri("aspect", "aspect", format = "ascii", env=env)
#asp=raster("aspect.asc")

#rsaga.wetness.index("dem30ag.tif","swi.sgrd",env = env)
#rsaga.sgrd.to.esri("swi", "swi", format = "ascii", env=env)
#swi=raster("swi.asc")

#rsaga.profile.curvature("dem30ag.tif", "cprof.sgrd", method = "poly2zevenbergen", env = env)
#rsaga.sgrd.to.esri("cprof", "cprof", format = "ascii", env=env)
#cpr=raster("cprof.asc")

#rsaga.plan.curvature("dem30ag.tif", "cplan", method = "poly2zevenbergen", env = env)
#rsaga.sgrd.to.esri("cplan", "cplan", format = "ascii", env=env)
#cpl=raster("cplan.asc")

#rsaga.hillshade("dem30ag.tif", "hill", method = "standard", env = env)
#rsaga.sgrd.to.esri("hill", "hill", format = "ascii", env=env)
#hill=raster("hill.asc")

#plot(hill, col=gray.colors(100), legend=F)
#plot(shape, add=T, cex=0.1)


s=stack(dem30ag,dem30agsd)

evi=raster('EVMMOD3a.tif')
pr2 <- projectRaster(dem30ag, crs='+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0')
evi=crop(evi, pr2)
evi <- projectRaster(evi, crs='+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
s=resample(s, evi)
s=stack(s, evi)

geo=raster('GEAISG3a.tif')
geo=crop(geo, pr2)
geo <- projectRaster(geo, crs='+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
s=stack(s, geo)

t=raster('tmed_anual.tif')
pr3 <- projectRaster(dem30ag, crs='+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0')
t=crop(t, pr3)
t <- projectRaster(t, crs='+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
s=resample(s, t)
s=stack(s, t)

pac=raster('prec_anual_acumulada.tif')
pac=crop(pac, pr3)
pac <- projectRaster(pac, crs='+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
s=resample(s, pac)
s=stack(s, pac)

bh=raster('BalHid.tif')
bh=crop(bh, pr3)
bh <- projectRaster(bh, crs='+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
s=resample(s, bh)
s=stack(s, bh)

lu=raster('GLCESA3a.tif')
lu=crop(lu, pr2)
lu=as.factor(lu)
lu <- projectRaster(lu, crs='+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',method='ngb')
lu=resample(lu, s,method='ngb')
soil=raster('STGHWS1a.tif')
soil=crop(soil, pr2)
soil=as.factor(soil)
soil <- projectRaster(soil, crs='+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs',method='ngb')
soil2=resample(soil, s,method='ngb')
soil2=as.factor(soil2)
plot(soil2)
plot(lu)
preds=stack(soil2,lu)
plot(preds)
CronometroOFF()

#Harmonize data sets

e=extract(s, shape)
dat=as.data.frame(shape)
dat=cbind(dat, e)
dat=data.frame(x=dat[,4], y=dat[,5],CO=dat[,2], scale(dat[,10:16]))
dat=dat[complete.cases(dat),]
cov=as.data.frame(s, xy=TRUE)
cov=cov[complete.cases(cov),]
str(dat)
str(cov)

#Descriptive statistics

summary(dat)

apply(dat[,1:9], 2, max)
apply(dat[,1:9], 2, min)
apply(dat[,1:9], 2, median)
apply(dat[,1:9], 2, mean)
apply(dat[,1:9], 2, sd)

par(mfrow=c(1,2))
plot(density(dat$CO), main='CO available data')
abline(v=mean(dat$CO), col='blue')
abline(v=median(dat$CO), col='blue', lty=2)

plot(density(log(dat$CO)), main=' log(CO) available data')
abline(v=mean(log(dat$CO)), col='blue')
abline(v=median(log(dat$CO)), col='blue', lty=2)

#Spatial autocorrelation

datsp=dat
coordinates(datsp)=~x+y
fit=autoKrige(CO~1, datsp)
kr = stack(fit$krige_output)
#plot(hill,col=gray.colors(100), main='CO kriging', legend=F)
plot(kr[[1]], col=rainbow(255), add=F, alpha=0.75)
plot(datsp, add=T)

#correlation matrix and exploratory PCA 
round(cor(dat),2)

pca=princomp(scale(dat))
summary(pca)
plot(pca)
biplot(pca, col=c('lightgray','blue'))
pca$loadings


#MODELOS

separar1 <- function(dataframe, seed=NULL,trainsize=0.50) {#from 0-1
  if (!is.null(seed)) set.seed(seed)
  index <- 1:nrow(dataframe)
  trainindex <- sample(index, round(trainsize*nrow(dataframe)))
  trainset <- dataframe[trainindex, ]
  testset <- dataframe[-trainindex, ]
  list(trainset=trainset,testset=testset)
  }

dat2=dat
#set.seed(999)
datosNummodel <- separar1(dat2) # se ejecuta la función para separar los datos
entrenamiento <- datosNummodel$trainset  # se crea la base de entranamiento
prueba <- datosNummodel$testset  # se crea la base de prueba
dat2$CO=log(dat2$CO)
lm1=step(lm(CO~., data=entrenamiento[,3:10]))
pru1=predict(lm1, prueba)
cor.test(prueba$CO, pru1)
rmse(exp(prueba$CO), exp(pru1))

datosNummodel <- separar1(dat) # se ejecuta la función para separar los datos
entrenamiento <- datosNummodel$trainset  # se crea la base de entranamiento
prueba <- datosNummodel$testset
lm2=step(lm(CO~., data=entrenamiento[c(3,6,8)]))
pru2=predict(lm2, prueba)
cor.test(prueba$CO, pru2)
rmse(exp(prueba$CO), exp(pru2))

lm1p=predict(scale(s), lm1)
#plot(hill,col=gray.colors(100), legend=F,main='CO lm')
plot(lm1p, col=rainbow(255), add=F, alpha=0.75)
barplot(vif(lm1), las=2)

?regsubsets
b<-regsubsets(CO~., data=entrenamiento[,3:10],nbest=2)
plot(b,scale='adjr2')

#UNCERTAINTY

res=data.frame(r_sq_cv=numeric(), rmse_cv=numeric() , r_sq_ext=numeric(), rmse_ext=numeric())

p=preds
for (i in 1:10){
datosNummodel <- separar1(dat2) # se ejecuta la función para separar los datos
entrenamiento <- datosNummodel$trainset  # se crea la base de entranamiento
prueba <- datosNummodel$testset  # se crea la base de prueba
lm1=step(lm(CO~., data=entrenamiento[,3:10]))
lm1p=predict(scale(s), lm1)
pru=predict(lm1, prueba)
lm1cv <- cv.lm(entrenamiento,lm1,m=10, plotit=FALSE)
res[i,1]=cor(entrenamiento$CO, lm1cv$cvpred)^2
res[i,2]=rmse(exp(entrenamiento$CO), exp(lm1cv$cvpred))
res[i,3]=cor(prueba$CO, pru)^2
res[i,4]=rmse(exp(prueba$CO), exp(pru))
p=stack(p,lm1p)
}

medianLm1=calc(p[[3:12]], median)
sdLm1=calc(p[[3:12]], sd)

myTheme <- rasterTheme(region =  dichromat(c('white', 'green','darkgreen', 'black')))
levelplot(medianLm1,par.settings = myTheme, maxpixels=391989, main='median CO')
x11()
levelplot(sdLm1, par.settings = myTheme,maxpixels=391989,main='uncertainty')

writeRaster(exp(medianLm1), file='medianLM.tif')
writeRaster(exp(sdLm1), file='uncertaintyLM.tif',overwrite=TRUE)

#
coordinates(dat)=~x+y
proj4string(dat)=crs('+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')
covsp=scale(s)
covsp=as.data.frame(covsp, xy=T)
covsp=covsp[complete.cases(covsp),]
coordinates(covsp)=~x+y
gridded(covsp)=T
proj4string(covsp)=crs('+proj=utm +zone=16 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs')

f = CO ~ AsterDEM_HN30m.2 + EVMMOD3a + 
    tmed_anual + prec_anual_acumulada + BalHid

omm <- fit.gstatModel(dat, f, covsp, 
   method="ranger")

plot(omm)

#plot(omm@regModel)
omm@vgmModel
om.rk <- predict(omm, covsp, nfold=10)
om.rk
plot(om.rk@predicted, main='OC random forest regression kriging')
plotKML(om.rk)


