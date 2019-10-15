#Organizing covariate files 

####
####COVARIATES SIB

lis <- list.files("SIB3_CARBON_FLUX_909/data")[seq(91, 107, by=2)]
lis2 <- paste0(getwd(),'/',"SIB3_CARBON_FLUX_909/data",'/',lis)
lis3 <- list.files(lis2, full.names=TRUE)

RES <- stack()
for (i in 1998:2006){
lis4 <-lis3[grep(i,  lis3)]
minSIB <- calc( stack(lis4) , min)
meaSIB <- calc( stack(lis4) , mean)
medSIB <- calc( stack(lis4)  , median)
maxSIB <- calc( stack(lis4) , max)
sdSIB <- calc( stack(lis4)  , sd)
res <- stack(minSIB, meaSIB, medSIB, maxSIB, sdSIB)
names(res) <-c(paste0('min', i), paste0('mea', i), paste0('med', i), paste0('max', i), paste0( 'sd', i))
RES <- stack(RES, res)
print(i)
}

saveRDS(as(RES, 'SpatialPixelsDataFrame'), file='SIB.rds')

####
####COVARIATES LAND USE 

lis <- list.files('LAND_USE_HARMONIZATION_V1_1248/data')
lis2 <- paste0(getwd(),'/','LAND_USE_HARMONIZATION_V1_1248/data','/',lis)[25:48]

RES2 <- stack()
for (i in 1:length(lis2)){
nam <- unlist(strsplit(lis2[i] , 'LUHa_t1.v1_'))[2]
nam1 <- unlist(strsplit(nam , '.nc4'))
x <- stack(lis2[i])[[299:306]]
names(x) <- c(paste0('lu_', nam1,'1998'), paste0('lu_', nam1,'1999'), paste0('lu_', nam1,'2000'), paste0('lu_',nam1, '2001'), paste0('lu_',nam1, '2002'), paste0('lu_',nam1, '2003'), paste0('lu_',nam1, '2004'), paste0('lu_',nam1, '2005'))
RES2 <- stack(RES2, x)
print(i)
}

saveRDS(as(RES2, 'SpatialPixelsDataFrame'), file='LU1.rds')

####
####COVARIATES BIOMASS

s <-stack("Global_Biomass_1950-2010_1296/data/historical_global_1-degree_forest_biomass.nc4", varname='Total_Biomass')


