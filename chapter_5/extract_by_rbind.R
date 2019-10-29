#...
#empty data frame
m <- readRDS(paste0( lis[1]))
	m1 <- raster::as.data.frame(m[1,][-1,])
	m1$SOC <- numeric()
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
#por si usan stack
#ov <- extract(m, d1)
	  ov$x <-raster::as.data.frame(d1)$X
	  ov$y <- raster::as.data.frame(d1)$Y
	  ov$SOC <- d1$OCSKGM 
	  ov <- ov[complete.cases(ov),]
	  m1 <- rbind(m1, ov)
	print(i)
}

