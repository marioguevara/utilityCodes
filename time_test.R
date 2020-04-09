#FUNCTION TO MEASURE TIME (CONABIO)
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
      
library(raster)
library(RSAGA)
R <- getData('alt', country='AUS', mask=TRUE)

###aqui cambio 60, 70 Y 100

R <- aggregate(R, 70, mean)

#R100
#class      : RasterLayer 
#dimensions : 55, 56, 3080  (nrow, ncol, ncell)
#resolution : 0.8333333, 0.8333333  (x, y)
#extent     : 112.8, 159.4667, -54.93333, -9.1  (xmin, xmax, ymin, ymax)
#crs        : +proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0 
#source     : memory
#names      : AUS_msk_alt 
#values     : 1, 1039.671  (min, max)


writeRaster(R, "a_AUS_DEM.tif", overwrite=TRUE)
#myenv <- rsaga.env(workspace =".", path="E:/saga-7.0.0_x64")
rsaga.import.gdal("a_AUS_DEM.tif", 'a_conus_albers_elev_.sgrd', env = rsaga.env())

CronometroON()
#LINEAS DE RICARDO APRA GENERACION DE ATRIBUTOS DE TERRENO EN CONUS
rsaga.geoprocessor("ta_compound", module = "Basic Terrain Analysis",
                   param = list(ELEVATION = "a_conus_albers_elev_.sgrd",
                                SHADE = "a_conus_albers_rsaga_analytical_hillshading_1km.sgrd",
                                SLOPE = "a_conus_albers_rsaga_slope_1km.sgrd",
                                ASPECT = "a_conus_albers_rsaga_aspect_1km.sgrd",
                                HCURV = "a_conus_albers_rsaga_plan_curvature_1km.sgrd",
                                VCURV = "a_conus_albers_rsaga_profile_curvature_1km.sgrd",
                                CONVERGENCE = "a_conus_albers_rsaga_convergence_index_1km.sgrd",
                                SINKS = "a_conus_albers_rsaga_closed_depressions_1km.sgrd",
                                FLOW = "a_conus_albers_rsaga_total_catchment_area_1km.sgrd",
                                WETNESS = "a_conus_albers_rsaga_topographic_wetness_index_1km.sgrd",
                                LSFACTOR = "a_conus_albers_rsaga_ls_factor_1km.sgrd",
                                CHNL_BASE = "a_conus_albers_rsaga_channel_network_base_level_1km.sgrd",
                                CHNL_DIST = "a_conus_albers_rsaga_channel_network_distance_1km.sgrd",
                                VALL_DEPTH = "a_conus_albers_rsaga_valley_depth_1km.sgrd",
                                RSP = "a_conus_albers_rsaga_relative_slope_position_1km.sgrd"),
                   show.output.on.console = TRUE, env = rsaga.env())
CronometroOFF()

