#get the libraries
library(raster)
 library(GSIF)
 #load the DEM
 DEM_TODO <- raster('DEM_T1.tif   ')
 #get the extent of the DEM
 e <- extent (DEM_TODO)
 #convert the extent in spatial polygons
 e <- as(e, 'SpatialPolygons')
 #set the projection of the polygon extent 
 proj4string(e) <- CRS(projection(DEM_TODO))
 #get the tiles of 0.5 degrees squares (in the units of the map)
tiles <- getSpatialTiles(e, 0.5, 0.5, 
         overlap.percent = 20, limit.bbox = TRUE, return.SpatialPolygons = TRUE)
#visualize the tiles
plot(tiles)
#overlap tile 1
plot(tiles[1], col='red', add=T)

###aqui empiezas el loop 

for ...

#crop to tile i
dem_tile <- crop(DEM_TODO, tile[i])

#generate your terrain parameters using your RSAGA code
...
#stack them and converto to spatial pixel data frame
...
#run the prediction
...
#save prediction for each tile
...
#end the loop

#Mosaic predited maps (code later)...






