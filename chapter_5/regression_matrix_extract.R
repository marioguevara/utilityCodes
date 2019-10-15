#This code is to generate a sampling strategy for a prediction domain.
#import datasets

#inputs


setwd("~/Downloads/chapter5")
time1 = 1997
time2 = 2008
train <- readRDS("soc_wosis_oct_2_2019.rds")
x <- readRDS("before_spline_soc_wosis_oct_2_2019.rds")
x <- x[c('id', 'year', 'date' )]

library(raster)
topo <- readRDS('topographic_predictors_15km_grids.rds')
sib <- readRDS('SIB.rds')
lu <- readRDS('LU1.rds')
bio <- stack("Global_Biomass_1950-2010_1296/data/historical_global_1-degree_forest_biomass.nc4", varname='Total_Biomass')
sm <- stack(list.files()[grep( 'sm_global', list.files() )])
names(sm)

dat <- merge(unique(x), train)
dat <- data.frame(dat, data.frame(extract(stack(topo), dat[c("X", "Y")])))
dat <- data.frame(dat, data.frame(extract(stack(sib), dat[c("X", "Y")])))
dat <- data.frame(dat, data.frame(extract(stack(lu), dat[c("X", "Y")])))
dat <- data.frame(dat, data.frame(extract(stack(bio), dat[c("X", "Y")])))
dat <- data.frame(dat, data.frame(extract(stack(sm), dat[c("X", "Y")])))
dat$year <- as.numeric(dat$year)
library(maps)

DAT <- data.frame()
for (i in 2005:2000){
train1 <- dat[dat$year==i,]
static <- train1[names(topo)]
nam <- train1[grep(i, names(train1))]
nam1 <- unlist(strsplit(names(nam) , i))
names(nam) <- nam1
d1 <- data.frame(scale(nam))
d1 <- data.frame(static, d1)
d1$SOC <- log1p(train1$SOC)
d1$year <- i
#d1 <- d1[ , colSums(is.na(d1)) == 0]
DAT <- rbind(DAT, d1)
}
#DAT <- DAT[ , colSums(is.na(DAT)) == 0]

e <- read.csv('prediction_domain_example.csv')
mp <- data.frame()
mp <- rbind(mp, data.frame(extract(stack(topo), e)))
mp <- data.frame(mp, data.frame(extract(stack(sib), e)))
mp <- data.frame(mp, data.frame(extract(stack(lu), e)))
mp <- data.frame(mp, data.frame(extract(stack(bio), e)))
mp <- data.frame(mp, data.frame(extract(stack(sm), e)))

preds <- data.frame()
for (i in 2005:2000){
static <- mp[names(topo)]
nam <- mp[grep(i, names(mp))]
nam1 <- unlist(strsplit(names(nam) , i))
names(nam) <- nam1
d1 <- data.frame(scale(nam))
d1 <- data.frame(static, d1)
d1$year <- i
preds <- rbind(preds, d1)
print(i)
}

is.nan.data.frame <- function(x)
do.call(cbind, lapply(x, is.nan))
preds[is.nan(preds)] <- -9999
preds[is.na(preds)] <- -9999
DAT[is.nan(DAT)] <- -9999
DAT[is.na(DAT)] <- -9999

DAT <- DAT[vapply(DAT, function(x) length(unique(x)) > 1, logical(1L))]

y <- DAT$SOC
x <- DAT[1:41]
	
###

#



