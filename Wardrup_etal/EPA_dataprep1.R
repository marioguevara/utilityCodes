
#read data

dat <- read.csv('nwca2011_soilchem.csv')

#generate columns to fill

dat$top <- dat$LAYER
dat$bottom <- dat$DEPTH

#first layer (1) will be 0 cm depth

dat$top[dat$top==1] <- 0 

#run a loop for each row 

for(i in 1:nrow(dat)){

#if not equal to 0

if (dat$top[i]!=0){

#top will be the value of bottom at the row i-1
  
dat$top[i] <- dat$bottom[i-1] 

}}

#verify values for the rows 21:40, explorre other rows

(dat)[111:112][21:40,]

#export as a csv file to verify in excel if needed

write.csv(dat, file='fixedTopEPA_v1.csv')

#import the coordinate site...

