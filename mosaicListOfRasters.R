ListRasters <- function(list_names) {
  raster_list <- list() # initialise the list of rasters
   for (i in 1:(length(list_names))){ 
    grd_name <- list_names[i] # list_names contains all the names of the images in .grd format
    raster_file <- stack(grd_name, bands=2)
   }
  raster_list <- append(raster_list, raster_file) # update raster_list at each iteration
}

wgs84.tif.list <- list.files(path=getwd(), pattern=glob2rx("*.tif"), full.names=T,recursive=F)

list_names <- NULL
for (i in 1:length(wgs84.tif.list)) {
  list_names <- c(list_names, wgs84.tif.list[i])
}

raster.list <-sapply(list_names, FUN = ListRasters)

raster.list$fun <- mean

names(raster.list) <- NULL

raster.list$fun <- mean
mos <- do.call(mosaic, raster.list)

writeRaster(mos, file='/home/work/spac/mguevara/UNCERTAINTYmodel3varQRFsocUSMX.tif')




