# Laura Jimenez; March, 2020

# Exploring Humminbirds data
setwd("C:/Users/l215j162/Documents/KU/Doc-project/NicheEstimation")

# Read file with all hummingbird species
hummer <- read.csv("./2019-20/Hummers-Cooper/IOC_hummers_Cooper2.csv",header=T)
head(hummer)
dim(hummer)
hummer <- na.omit(hummer)
head(hummer)
dim(hummer)
# Select subset of species
hummer.set <- subset(hummer,(hummer[,5]>=300)&(hummer[,5]<1000))
head(hummer.set)
dim(hummer.set)
# create vector of species names
sp.names <- paste(hummer.set[,1],hummer.set[,3],sep = "_")
write.csv(cbind(sp.names,hummer.set$OccPnts),"./2019-20/Hummers-Cooper/AnalysisTwo/occs-300-1000/spnames.csv",row.names = F)
### END