
# working directory
setwd("C:/Julian/133_errores_Mario")

# packages
library(raster)
library(randomForest)
library(quantregForest)

# read training data
dat<- read.csv(file="dat.csv",head=TRUE,sep=",")

# read complete data
region <- read.table("COV1.txt",sep="\t",header=TRUE)

# data structure
head(dat)

# only keep the data I give a shit about
dat=regM[-c(2,3)]

# quantile random forest model (Warning: quantile rf is slow)
quant_rf_model <- quantregForest(y=dat$NumEspecies, x=dat[,2:7])

# plot out-of-bag predictions (for training data)
plot(quant_rf_model)

# assuming a normal distribution which is obviously wrong:

### for 1 standard deviation you would have
pnorm(1, mean = 0, sd = 1)-pnorm(-1, mean = 0, sd = 1)
### this is the probability of the interval -1 sd (---mean---) +1 sd
### is 0.682

### for 2 standard deviations you would have
pnorm(2, mean = 0, sd = 1)-pnorm(-2, mean = 0, sd = 1)
### this is the probability of the interval -2 sd (---mean---) +2 sd
### is 0.954

### althought this was an example based on a standard normal distribution
### is is the same for normal distributions in general
### 68% of data is contained at 1 standard distribution of the mean
### 95% of data is contained at 2 standard distributions of the mean

# prediction for quantile rf model for 1 standard deviation (HORRIBLY SLOW)
pred_1sd <- predict(quant_rf_model,
                    xx, # ***
                    quantiles=c((1-.682)/2, 1-(1-.682)/2))

# *** this function sucks and covariables need to be in the same 
#     order as trainig data

# L?mites
Lim_sup_1sd <- pred_1sd[,3]
Lim_inf_1sd <- pred_1sd[,1]

# "incertidumbre"
var_1sd <- abs(Lim_sup_1sd-Lim_inf_1sd) # absolute difference

# usual random forest model
rf_model <- randomForest(y=dat$NumEspecies, x=dat[,2:7],ntree=1000,nodesize=10)

# OOB absolut error
oob_abs <- abs(dat$NumEspecies-rf_model$predicted) # absolute difference

# associate uncertainty with OOB squared error
oob_abs_database <- data.frame(oob_abs,x=regM$x,y=regM$y)
var_1sd_database <- data.frame(var_1sd,x=xx$x,y=xx$y)

# merge databases
merged_oob_var <- merge(oob_abs_database,var_1sd_database,by=c("x","y"))

# were they merged correctly?
head(merged_oob_var)
              # no, meaning projections suck and points dont match

# attempt TWO 

### turn var_1sd_database into a raster and extract info
var_1sd_rast <- var_1sd_database
coordinates(var_1sd_rast)=~x+y
gridded(var_1sd_rast)=TRUE
var_1sd_rast <- raster(var_1sd_rast)

### visualize and save if wanted
#plot(var_1sd_rast)
#writeRaster(var_1sd_rast, filename="uncertainty.tif", format="GTiff", overwrite=TRUE)

### turn oob to spatial points
oob_abs_sp <- oob_abs_database 
coordinates(oob_abs_sp) = ~x+y

### extract info from raster
extraction <- extract(var_1sd_rast,oob_abs_sp)

### structure of extraction
length(extraction)

### correlation between uncertainty and prediction errors
cor(extraction,oob_se,use="pairwise.complete.obs" )
# J: 0.1822 correlation small but positive

# create the rest of the rasters (original,lim_sup, lim_inf)

#### lim_sup
Lim_sup_1sd_database <- data.frame(Lim_sup_1sd,x=xx$x,y=xx$y)
Lim_sup_1sd_rast <- Lim_sup_1sd_database
coordinates(Lim_sup_1sd_rast)=~x+y
gridded(Lim_sup_1sd_rast)=TRUE
Lim_sup_1sd_rast <- raster(Lim_sup_1sd_rast)

### visualize and save if wanted
plot(Lim_sup_1sd_rast)
writeRaster(Lim_sup_1sd_rast, filename="Lim_sup.tif", format="GTiff", overwrite=TRUE)

#### lim_inf
Lim_inf_1sd_database <- data.frame(Lim_inf_1sd,x=xx$x,y=xx$y)
Lim_inf_1sd_rast <- Lim_inf_1sd_database
coordinates(Lim_inf_1sd_rast)=~x+y
gridded(Lim_inf_1sd_rast)=TRUE
Lim_inf_1sd_rast <- raster(Lim_inf_1sd_rast)

### visualize and save if wanted
plot(Lim_inf_1sd_rast)
writeRaster(Lim_inf_1sd_rast, filename="Lim_inf.tif", format="GTiff", overwrite=TRUE)

#### Original
original_database <- data.frame(predict(rf_model,xx),x=xx$x,y=xx$y)
original_rast <- original_database
coordinates(original_rast)=~x+y
gridded(original_rast)=TRUE
original_rast <- raster(original_rast)

### visualize and save if wanted
plot(original_rast)
writeRaster(original_rast, filename="original.tif", format="GTiff", overwrite=TRUE)
