###
###NDVI 
path <- "/home/mguevara/Downloads/SALINITY/Task 1-20201007T192716Z-001/Task 1"
data<- readRDS('/home/mguevara/Downloads/SALINITY/coachella_015ECe_ECa.rds')
countyName <- 'Riverside'
####
data_ECe <- data[[1]]
data_ECa <- data[[2]]
(lis <- list.files(full.names=TRUE)[-1])

pas <- paste0(lis[1], '/' ,list.files(lis[1], pattern='shp'))
setwd(path)
library(raster)
sh <- shapefile(pas)
county <- sh[sh@data$County==countyName,]

county <- spTransform(county, CRS(projection(data_ECe)))
#county_subset <- county[data_ECe,]
ov <- over(data_ECe, county)
data_ECe@data  <- cbind(data_ECe@data, ov)
#data_df <- data_ECe@data
data_df <- na.omit(data_df)

data_df <- data_ECe@data[c(2, 4, 14)]


data_df$ECe0.15_log <- log1p(data_df$ECe0.15)
data_df$avg_ndvi_log <- log1p(data_df$avg_ndvi)
##

library( ggpubr)

ggscatter(
  data_df, x = "ECe0.15", y = "avg_ndvi",
  color = "Crop2014", palette = "jco",
  add = "reg.line"
  ) +
  facet_wrap(~Crop2014) #+
  #stat_cor(label.y = 4.4) +
  #stat_regline_equation(label.y = 4.2)

ggplot(data_df) +
  geom_jitter(aes(ECe0.15_log, avg_ndvi_log, colour=Crop2014),) + geom_smooth(aes(ECe0.15_log, avg_ndvi_log, colour=Crop2014), method=lm, se=FALSE) +
  facet_wrap(~Crop2014, scales="free_x") +
   theme_bw(20) 
###
ggplot(data_df) +
  geom_jitter(aes(ECe0.15, avg_ndvi, colour=Crop2014),) + geom_smooth(aes(ECe0.15_log, avg_ndvi_log, colour=Crop2014), method=lm, se=FALSE) +
  facet_wrap(~Crop2014, scales="free_x") +
   theme_bw(20) 
###

ggplot(data_df, aes(ECe0.15_log, avg_ndvi_log, shape=Crop2014, colour=Crop2014, fill=Crop2014)) +
   geom_smooth(method="lm") +
   geom_point(size=3) +
   theme_bw(20) 


###Coachella valley
path <-  "/home/mguevara/Downloads/SALINITY/GroundTruth_Data-20200827T035123Z-001/GroundTruth_Data/CoachellaValley/CSV files - soil _ ECa"
setwd(path)
library(rgdal)
ind <- 5
res<-data.frame()
for(i in 1:length(list.files(pattern='Pro.csv'))){
print(nam <- names(read.csv(list.files(pattern='Pro.csv')[i])[1:ind]))
res <- rbind(res, nam)
#dat <- read.csv(list.files(pattern='Pro.csv')[i])
#res <- rbind(res, dat[nam])
}

nam <- c("ID" ,     "X" ,    "Y" ,   "ECe0.15" )


res<-data.frame()
lis <- list()
for(i in 1:length(list.files(pattern='Pro.csv'))){
dat <- read.csv(list.files(pattern='Pro.csv')[i])
idx <- grep("ECe0.15", names(dat))
dat2 <-dat[idx]
s1 <-dim(dat2)[2]==0
if (s1==FALSE){
print(dim(dat))
res <- rbind(res, dat[nam])
 }
else{print(names(dat))
}

}

summary(res)

res$Y[res$Y==0] <- NA
res <- na.omit(res)

coordinates(res) <- ~ X+Y
proj4string(res) <-CRS("+proj=utm +zone=11+datum=WGS84")

library(automap)
 plot(variogram)
#variogram = autofitVariogram(ECe0.15 ~ 1, res2)





nam <- c("ID" ,     "X" ,    "Y" ,   "EMh" )
res2<-data.frame()
lis <- list()
for(i in 1:length(list.files(pattern='ECa'))){
dat <- read.csv(list.files(pattern='ECa.csv')[i])
idx <- grep("EMh", names(dat))
dat2 <-dat[idx]
s1 <-dim(dat2)[2]==0
if (s1==FALSE){
print(dim(dat))
res2 <- rbind(res2, dat[nam])
 }
else{print(names(dat))
}
}

summary(res2)

res2$Y[res2$Y==0] <- NA
res2 <- na.omit(res2)

coordinates(res2) <- ~ X+Y
proj4string(res2) <-CRS("+proj=utm +zone=11+datum=WGS84")
library(automap)
#variogram2 = autofitVariogram(EMh ~ 1, res2)



saveRDS(list(res, res2), file='/home/mguevara/Downloads/SALINITY/coachella_015ECe_ECa.rds')


library(soilDB)
library(aqp)
library(elevatr)
library(rgdal)
library(raster)
r <- readRDS('/home/mguevara/Downloads/SALINITY/coachella_015ECe_ECa.rds')
r1 <- r[[1]]
r1 <- spTransform(r1, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs "))
x <- mapunit_geom_by_ll_bbox(bbox(r1))

X <- x[r2,]
plot(X)

invisible(text(getSpPPolygonsLabptSlots(X),
labels = as.character(X$musym), cex = 1), col='red')

r1@data <- cbind(r1@data, soil=over(r1, x)$musym)
r2 <- r[[2]]
r2 <- spTransform(r2, CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs "))
#r2@data <- cbind(r2@data, soil=over(r2, x)$musym)

#get elevation data for that area
elev <- get_elev_raster(r2, prj = projection(r2), z = 10, clip = "bbox")

x@data$soilmap <- as.factor(x@data$musym)

soilmap <- rasterize(x = x, y = elev, field = "soilmap")

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

soilmap_dummy <- dummyRaster(soilmap)
#library(RStoolbox)
#rpc <- rasterPCA(soilmap_dummy)
sr <- sampleRandom(soilmap_dummy, 50000)
pca <- prcomp(sr, scale=TRUE, retx=FALSE) 
pcas <- predict(soilmap_dummy, pca, index=1:10)
pcas <- resample (pcas, elev)
covs <- stack(pcas, elev)

names(covs)[1:10] <- paste0('pca_', 1:10) 
names(covs)[11] <- 'elev'

r2@data <- cbind(r2@data, data.frame(extract(covs, r2)))
names(r2)



#r1@data$EC <- log1p(r1@data$ECe0.15
#r2@data$EC <- r2@data$EMh
#r1@data$method <- 'FLD'
#r2@data$method <- 'EMI'


data <- rbind(r1@data[3:5], r2@data[3:5])



levels(as.factor(x$musym))


#####
#####
#####

library(rgp)

newFuncSet <- functionSet("+","-","*","/","sqrt","exp","ln")

result1 <- symbolicRegression(ECe0.15~EMh,
                              data=dataset, functionSet=newFuncSet,
                              stopCondition=makeStepsStopCondition(300))

pred <- predict(result1, dataset, model = "BEST",  detailed = TRUE)

res <- spTransform(res, CRS('+proj=longlat +datum=WGS84 +no_defs  '))

names(dat)

s1 <- nam==names(dat[nam])
s2 <- s1==TRUE
s3 <- grep('TRUE', s2)


if (length(s3)==4){
 res <- rbind(res, dat[nam])
  } else {
print(names(dat))
idx[i] <- i
}

Y <- res@data[sample(nrow(res@data), 500), ][2]
X  <- res2@data[sample(nrow(res2@data), 500), ][2]
dataset <- cbind(X, Y)

x <- log1p(dataset$EMh)
y <- log1p(dataset$ECe0.15)

library()

#
#San Joaquin

#path <- '/home/mguevara/Downloads/SALINITY/GroundTruth_Data-20200827T035123Z-001/GroundTruth_Data/SanJoaquinValley'
path <- "/home/mario/Downloads/salinity_datasets/GroundTruth_Data/SanJoaquinValley"

setwd(path)
lis <- list.files(list.files(), full.name=TRUE)
nam <- c("ID" , "X" ,  "Y"  , "EMv" , "EMh")
idx <- grep('EMI', lis)
library(readxl)
res <- data.frame()
for (i in 1:length(lis[idx])){
#for (i in 7:7){
data <- read_excel(lis[idx][i])
names(data)[1:5]<- nam
print(i)
#print(summary(data[c('X','Y')]))
data <- data[nam]
res <- rbind(res, data)
}
library(rgdal);library(mapview)
coordinates(res) <- ~ X+Y
plot(res)
x1 <- crop(res, drawExtent())
x2 <- crop(res, drawExtent())
proj4string(x1) <-CRS("+proj=utm +zone=10+datum=WGS84")
proj4string(x2) <-CRS("+proj=utm +zone=11+datum=WGS84")
x1 <- spTransform(x1, CRS='+proj=longlat +datum=WGS84')
x2 <- spTransform(x2, CRS='+proj=longlat +datum=WGS84')
x <- rbind(x1, x2)
dim(x)
mapview(x)
###
  ####
  path <- "/home/mario/Downloads/salinity_datasets/GroundTruth_Data/SanJoaquinValley"
  setwd(path)
  lis <- list.files(list.files(), full.name=TRUE)
  #nam <- c("ID" , "X" ,  "Y"  , "EMv" , "EMh")
  #nam2 <- c("Field ID"  ,  "X (UTM)"   ,  "Y (UTM)" ,    "EMv"  ,    "EMh"     ,"Depth (cm.)","ECe (dS/m)")
  #idx <- grep('EMI', lis)
  library(readxl)
  length(lis[idx])
  lis2 <- lis[-idx]
  #lis2 <- lis[c(2,4)]
  #res <- data.frame()
  res1 <- read_excel(lis2[1])
  res1 <- data.frame(Field=res1$FIELD, X=res1$x, Y=res1$y, EMv=res1$EMv, EMh=res1$EMh, ECe=res1$ECe, SP=res1$SP, DEPTH=res1$Depth)

  res3 <- read_excel(lis2[3])
  res3 <- data.frame(Field=res[1], X=res3[,3], Y=res3[,4], EMv=res3$Emv, EMh=res3$EMh, ECe=res3[,10], SP=res3[,9], DEPTH=res3[2])
  names(res3) <-names(res1)


  res <- rbind(res1, res3)

  lis2 <- lis2[-c(1, 3)]

  for (i in 1:length(lis2)){

    #for (i in 7:7){

    data <- read_excel(lis2[i])
    Field <- data[grep('Field' , names(data)) ] 
    #X <- data[grep('x' , names(data)) ] # i = 1 = x, y i > = X, Y
    #Y <- data[grep('y' , names(data)) ] # i = 1 = x, y i > = X, Y
    X <- data[grep('X' , names(data)) ]
    Y <- data[grep('Y' , names(data)) ]
    EMv <- data[grep('EMv' , names(data)) ]


    EMh <- data[grep('EMh' , names(data)) ]
    ECe <- data[grep('ECe' , names(data)) ]
    SP <- data[grep('SP' , names(data)) ]
    DEPTH <- data[grep('Depth' , names(data))[1]]

  data <- cbind(Field, X,Y, EMv, EMh, ECe, SP, DEPTH)
  names(data) <- names(res)
  #names(data)[1:5]<- nam
  print(i)
  print(names(data) )
  #print(head(data))
  #print(summary(data[c('X','Y')]))
  #data <- data[nam2]
  res <- rbind(res, data)
  }

x1 <- crop(res, drawExtent())
x2 <- crop(res, drawExtent())
proj4string(x1) <-CRS("+proj=utm +zone=10+datum=WGS84")
proj4string(x2) <-CRS("+proj=utm +zone=11+datum=WGS84")
x1 <- spTransform(x1, CRS='+proj=longlat +datum=WGS84')
x2 <- spTransform(x2, CRS='+proj=longlat +datum=WGS84')
x <- rbind(x1, x2)
dim(x)
mapview(x)

res <- spTransform(x, CRS='+proj=longlat +datum=WGS84')
plot(res)
 
#saveRDS(res, 'San_Joaquin.rds')

res <-readRDS('San_Joaquin.rds')
res$DEPTH <- as.character(res$DEPTH)
res1 <- res[grep("\"" , as.character(res$DEPTH)),] ###300 points contain this pattern in depth column "\""

strsplit(res1$DEPTH[1], "\"")

res <- res[-grep("\"" , as.character(res$DEPTH)),]

#lisDEPTH <- strsplit(res1$DEPTH[i], "\"")
DEPTH <- character()
for(i in 1:length(res1$DEPTH)){
  DEPTH[i] <- unlist(strsplit(res1$DEPTH[i], "\""))
  }
 res1$DEPTH <- DEPTH

str(res); str(res1)

res <- rbind(res, res1)

# some issues with depth columns: "0-12\""  "12-24\"" "24-36\"" "36-48\""[22] "48-60\"" "48-54\"" "48-59\""

ch <- as.character(res$DEPTH)

top <- numeric()
bot <- numeric()
for(i in 1:length(ch)){
  ch_i <- ch[i]
  top[i] <- strsplit(ch_i[1], '-')[[1]][[1]]
  bot[i] <- strsplit(ch_i[1], '-')[[1]][[2]]
  }

res$top <- as.numeric(top)
res$bot <- as.numeric(bot)

summary(res@data)

res_ft <- res[c(1527:1556),]
res_ft$top <- res_ft$top*100
res_ft$bot <- res_ft$bot*100
res <- res[-c(1527:1556),]
res <- rbind(res, res_ft)


(readRDS('sanJoaquin_v0.rds'))

########
#####
###SJRIP

path <- "/home/mario/Downloads/salinity_datasets/GroundTruth_Data/SJRIP data_Amninder-20201020T174419Z-001/SJRIP data_Amninder"
setwd(path)
lis2 <- list.files(paste0(list.files(full.name=TRUE),'/', list.files(list.files(full.name=TRUE))), full.names=TRUE)
ECa <- data.frame()
ECe <- data.frame()
lis <-lis2[-c(19,27)]#19 and 27 files have inconsistency in first row 
for (i in 1:length(lis)){ 
dat <- read.delim(lis[i], header=F)
len <-  length(names(dat))
if ((len==1)==TRUE) {
dat <- read.delim(lis[i], header=F, sep=',')
ECa <- rbind(ECa, dat)
} else {
ECe <- rbind(ECe, dat)}
}

#1st row inconsistency 
dat19 <- dat <- read.delim(lis2[19], header=F, sep=',')[-1,]
head(dat19)
dat19$V5 <- 'NA'
names(dat19) <- c('X', 'Y','EMv', 'EMh', 'id' )

dat27 <- dat <- read.delim(lis2[27], header=F, sep=',')[-1,]
head(dat27)
dat27$V5 <- 'NA'
names(dat27) <- c('X', 'Y','EMv', 'EMh', 'id' )

#dat <- read.delim(list.files(paste0(list.files(full.name=TRUE)[1],'/', list.files(list.files(full.name=TRUE)[1])[1]), full.names=TRUE)[2], header=F)
names(ECe) <- c('id', 'DEPTH', 'ECe', 'PS', 'GWC')
#dat_a <- read.delim(list.files(paste0(list.files(full.name=TRUE)[1],'/', list.files(list.files(full.name=TRUE)[1])[1]), full.names=TRUE)[1], header=F, sep=',')
names(ECa) <- c('X', 'Y','EMv', 'EMh', 'id' )

ECa <- rbind(ECa, dat19)
ECa <- rbind(ECa, dat27)
ECa$X <- as.numeric(ECa$X)
dim(ECa)
ECa <- na.omit(ECa)
coordinates(ECa) <- ~ X+Y
proj4string(ECa) <-CRS("+proj=utm +zone=10+datum=WGS84")
ECa <- spTransform(ECa, CRS='+proj=longlat +datum=WGS84')

ECe$bot <- 0
ECe$bot[ECe$DEPTH==60] <-30
ECe$bot[ECe$DEPTH==90] <-60
ECe$bot[ECe$DEPTH==120] <-90

ECe$top <- ECe$DEPTH
ECe$top <- as.numeric(ECe$DEPTH)
ECe <- na.omit(ECe)

ECa$id <- row.names(ECa)
##
##



path <- "/home/mario/Downloads/salinity_datasets/GroundTruth_Data/SJRIP data_Amninder-20201020T174419Z-001/SJRIP data_Amninder"
setwd(path)
library(rgdal)
library(mapview)
lis3 <- paste0(list.files(full.name=TRUE),'/', list.files(list.files(full.name=TRUE)))
#list.files(lis3[1], full.name=TRUE)[1]
#list.files(lis3[1], full.name=TRUE)[2]
lis3 <- lis3[-c(12, 11, 9, 8, 5)] #inverse order due specific syntax
#lis3 <- lis3[-6]
ECe <- data.frame()
ECa <- data.frame()
for (i in 1:length(lis3)){
#for (i in 1:10){
#i=6
dat <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[1],header=F,)
len <-  length(names(dat))
if ((len==1)==TRUE) {
dat <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[1], header=F, sep=',')}
dat2 <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[2], header=F, sep=',')
len <-  length(names(dat2))
if ((len==1)==TRUE) {
dat2 <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[2], header=F)}
print(i)
print(lis3[i])
names(dat) <- c('X', 'Y','EMv', 'EMh', 'id' )
names(dat2) <- c('id', 'DEPTH', 'ECe', 'PS', 'GWC')
print(head(dat))
print(head(dat2))
dat$idx <- as.factor(row.names(dat))
dat2$idx <- as.factor(dat2$id)
  dat$field_season <- lis3[i]
  dat2$field_season <- lis3[i]
mer <- merge(dat2, dat, by='idx')
ECe <- rbind(ECe, mer)
ECa <- rbind(ECa, dat)
}











lis3 <- paste0(list.files(full.name=TRUE),'/', list.files(list.files(full.name=TRUE)))
lis3 <- lis3[c(12, 11, 9, 8, 5)] #inverse order due specific syntax
#lis3 <- lis3[-6]
lis3 <- lis3[-4]# inconsistency in line 1
ECe_rest <- data.frame()
ECa_rest <- data.frame()
for (i in 1:length(lis3)){
#i=4
dat <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[1],header=F)
len <-  length(names(dat))
if ((len==1)==TRUE) {
dat <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[1], header=F, sep=',')}
dat2 <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[2], header=F, sep=',')
len <-  length(names(dat2))
if ((len==1)==TRUE) { 
dat2 <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[2], header=F)}
print(i)
print(lis3[i])
names(dat2) <- c('X', 'Y','EMv', 'EMh', 'id' )
names(dat) <- c('id', 'DEPTH', 'ECe', 'PS', 'GWC')
print(head(dat))
print(head(dat2))
  print(dim(dat))
  dat2$idx <- as.factor(row.names(dat2))
dat$idx <- as.factor(dat$id)
  dat$field_season <- lis3[i]
  dat2$field_season <- lis3[i]
mer <- merge(dat, dat2, by='idx')
ECe_rest <- rbind(ECe_rest, mer)
ECa_rest <- rbind(ECa_rest, dat2)
}

ECe <- rbind(ECe, ECe_rest)
ECa <- rbind(ECa, ECa_rest)


coordinates(ECa) <- ~ X+Y
proj4string(ECa) <-CRS("+proj=utm +zone=10+datum=WGS84")
ECa <- spTransform(ECa, CRS='+proj=longlat +datum=WGS84')
mapview(ECa)

coordinates(ECe) <- ~ X+Y
proj4string(ECe) <-CRS("+proj=utm +zone=10+datum=WGS84")
ECe <- spTransform(ECe, CRS='+proj=longlat +datum=WGS84')
mapview(ECe)

lis3 <- paste0(list.files(full.name=TRUE),'/', list.files(list.files(full.name=TRUE)))
lis3 <- lis3[c(12, 11, 9, 8, 5)] #inverse order due specific syntax
i=4
dat <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[1],header=F)
len <-  length(names(dat))
if ((len==1)==TRUE) {
dat <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[1], header=F, sep=',')}
dat2 <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[2], header=F, sep=',')
len <-  length(names(dat2))
if ((len==1)==TRUE) {
dat2 <- read.delim(paste0(lis3[i], '/', list.files(lis3[i]))[2], header=F)}
print(i)
print(lis3[i])
dat <- dat[-1,] # remove the row error
names(dat2) <- c('X', 'Y','EMv', 'EMh', 'id' )
#names(dat) <- c('id', 'DEPTH', 'ECe', 'PS', 'GWC')
print(head(dat))
print(head(dat2))
##same thing dat and dat2
dat2$idx <- as.factor(row.names(dat2))
dat2$field_season <- lis3[i]
coordinates(dat2) <- ~ X+Y
proj4string(dat2) <-CRS("+proj=utm +zone=10+datum=WGS84")
dat2 <- spTransform(dat2, CRS='+proj=longlat +datum=WGS84')

ECa <- rbind(ECa, dat2)
mapview(ECa)



library(aqp)

ECe_aqp <- ECe
depths(ECe_aqp) <- id ~ top + bot

agg <- slab(ECe_aqp, fm= ~ ECe + PS)
 xyplot(top ~ p.q50 | variable, data=agg, ylab='Depth',
              xlab='median bounded by 25th and 75th percentiles',
              lower=agg$p.q25, upper=agg$p.q75,
              panel=panel.depth_function,
              alpha=0.25, sync.colors=TRUE,
              par.settings=list(superpose.line=list(col='RoyalBlue', lwd=2)),
              prepanel=prepanel.depth_function,
              cf=agg$contributing_fraction, cf.col='black', cf.interval=5, 
              layout=c(2,1), strip=strip.custom(bg=grey(0.8)),
              scales=list(x=list(tick.number=4, alternating=3, relation='free'))
              )


#> names(dat19)
#[1] "V1" "V2" "V3" "V4"
#> head(dat)
#             V1      V2      V3      V4
#1 LINE 1             NA      NA      NA
#2    713849.873 4085427 275.234 210.742
#3    713849.804 4085429 281.602 217.109
#4    713849.693 4085431 279.805 215.039
#5    713849.549 4085434 281.523 219.922
#6    713849.315 4085436 293.516 235.195
#> dat <- read.delim(lis[27], header=F, sep=',')
#> head(dat)
#             V1      V2      V3      V4
#1 LINE 1             NA      NA      NA
#2    713159.670 4084598 323.242 262.500
#3    713159.708 4084600 324.922 263.086
#4    713159.909 4084602 329.648 265.859
#5    713160.039 4084605 327.578 264.023
#6    713160.011 4084607 326.680 266.211



