#librerias y directorio de trabajo

setwd("~/Documentos/LST_3band_emissivity_Daily_1km_night_v61/LST_1KM")
#install.packages("remotes")
#library(remotes)
#install_github("ropensci/MODIStsp")

library(MODIStsp)
library(raster)
library(terra)

MODIStsp_get_prodnames()

lim <- gadm(country="GTM", level=1, path=tempdir())
box = bbox(as(lim, 'Spatial'))


MODIStsp(gui             = FALSE,
         out_folder      = '/home/mario/Descargas/GTM',
         out_folder_mod  = '/home/mario/Descargas/GTM',
         selprod         = "LST_3band_emissivity_Daily_1km_night (M*D21A1N)",
         bandsel         = 'LST_1KM',
         sensor          = 'Terra',
         user            = 'mguevara@udel.edu' , # your username for NASA http server
         password        = 'Projectohel@das2024',  # your password for NASA http server
         start_date      = '2024.01.01',
         end_date        = '2024.06.15',
         verbose         = TRUE,
         bbox            =  box, #bbox of Latam
         spatmeth        = 'bbox',
         out_format      = 'GTiff',
         compress        = 'LZW',
         out_projsel     = 'User Defined',
         output_proj     = "+proj=longlat +datum=WGS84 +no_defs",
         delete_hdf      = TRUE,
         scale_val       = TRUE,
         parallel        = TRUE
)



