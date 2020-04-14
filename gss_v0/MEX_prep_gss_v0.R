#libraries
library(raster)
library(rgdal)
#read dataset in a shapefile 
shape <- readOGR(dsn = ".", layer = "edaf_puntos_sii")
proj4string(shape)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
#Mexico limit 
#country limit from the global administrative areas project\
#https://gadm.org/
lim <- readRDS('gadm36_MEX_2_sp.rds')
lim <- spTransform(lim, CRS(projection(shape)))
#lista de covariables
lis <- readRDS('covariateList.rds')
#visualize the points
#figure 1
x <- shape

#x <- x[is.na(x@data$PH)==FALSE,]
bubble(x, "PH",
        panel=function(...) {
          sp.polygons(lim, fill='white')
          sp:::panel.bubble(...)
        }) 

#and histograms
hist(shape$PH)
#log transformed
hist(log1p(shape$PH))
#generate a year only colum
year <- strsplit(as.character(shape$FECHA), '/')
y <- numeric()
for (i in 1:length(year)){
y[i] <- year[[i]][3]} 
shape$year <- as.numeric(y)
#select only those years after 1998 
shape@data <- shape@data[shape@data$year > 1998,]
shape@data <- na.omit(shape@data)
#algorthms for quantitative pedology
library(aqp)
sp4=shape@data
sp4$IDPROF <- paste0("IDPROF_", sp4$COORD_Y, "_", sp4$COORD_X)
#generate a soil profile collection object
depths(sp4) <- IDPROF  ~ LIM_SUP + LIM_INF
site(sp4) <- ~ COORD_X + COORD_Y
coordinates(sp4) <- ~ COORD_X + COORD_Y
proj4string(sp4)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
library(GSIF)
try(ph <- mpspline(sp4, 'PH', d = t(c(0,100))))

dat <- data.frame(id = sp4@site$IDPROF,
X = sp4@sp@coords[,1],
Y = sp4@sp@coords[,2],
ph = ph$var.std[,1])

#empty data frame
m <- readRDS(paste0( lis[1]))
	m1 <- raster::as.data.frame(m[1,][-1,])
	m1$ph <- numeric()
#prepare for extraction
d1 <- dat
coordinates(d1) <- ~X+Y
proj4string(d1)<-crs('+proj=lcc +lat_1=17.5 +lat_2=29.5 +lat_0=12 +lon_0=-102 +x_0=2500000 +y_0=0 +ellps=WGS84 +units=m +no_defs')
d1 <- spTransform(d1, CRS(projection(m)))
#loop for extraction
for(i in 1:length(lis)){
m <- readRDS(lis[i])
 proj4string(m) <- CRS( "+proj=longlat +ellps=GRS80 +no_defs")
ov <- over(d1, m)
	  ov$x <-raster::as.data.frame(d1)$X
	  ov$y <- raster::as.data.frame(d1)$Y
	  ov$ph <- d1$ph
	  ov <- ov[complete.cases(ov),]
	  m1 <- rbind(m1, ov)
	print(i)
}

