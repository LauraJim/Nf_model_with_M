# Laura Jimenez
# First version: May, 2020
# Last modification: November, 2020
# Project: Second model for the fundamental niche, Vespa mandarinia exercise

# PART 1: Occurrence data --------------------------------
### R code to create a matrix with the following columns: long, lat, bio1, bio12 (or more variables)
### Input data: an occurrence csv file and as many raster files as environmental variables

# Upload packages
library(dismo)
library(rgdal)

# Read file with occurrence points -----------
occ <- read.csv("./Wasps/occs_noEurope_50kmthin.csv",header=T)
head(occ)

# Read files with environmental layers
bio1 <- raster("C:\\Users\\l215j162\\Documents\\KU\\Doc-project\\ENM-group\\Wasps\\merra_world\\bio1.asc")
bio2 <- raster("C:\\Users\\l215j162\\Documents\\KU\\Doc-project\\ENM-group\\Wasps\\merra_world\\bio12.asc")
bio.names <- c("Bio1","Bio12")

# Create a single raster with as many layers as environmental variables
bios1_2 <- stack(bio1, bio2)
# Plotting the environmental layers
x11()
plot(bios1_2)

# Extract environmental values for each occurrence point
occ1 <- extract(bios1_2,occ)
occ2 <- na.omit(cbind(occ,occ1))
colnames(occ2) <- c("long","lat",bio.names)
head(occ2)
dim(occ2)

# Save the resulting matrix
write.csv(occ2,file=paste0("./Wasps/Models-set3/occs_noEurope_50kmthin_",paste(bio.names[1],bio.names[2],sep = "_"),".csv"),row.names = F)


# PART 2: M region -------------------------
# Get a random sample of points inside the polygon that delimits M and extract their environmental values

# Read shapefile that contains M polygon
M.shp <- readOGR("./Wasps/M","M_b4000")
# crop and mask the M polygon with the environmental layers
clip.M <- mask(crop(bios1_2,M.shp),M.shp)

x11()
plot(clip.M)

# Take a random sample of size N of all the cells in the raster
N <- 5000
sam.M <- sampleRandom(clip.M,N,xy=T)
colnames(sam.M) <- c("long","lat",bio.names)
head(sam.M)

# check sampled points in G-space
x11()
plot(sam.M[,1:2],pch=20)
points(occ2[,1:2],pch=20,col="red")
# check sampled points in E-space
x11()
plot(sam.M[,3:4],pch=20)
points(occ2[,3:4],pch=20,col="red")

# Save the matrix of geographic coordinates and environmental values
write.csv(sam.M,file=paste0("./Wasps/Models-set3/M4000_with_",paste(bio.names[1],bio.names[2],sep = "_"),".csv"),row.names = F)

# PART 3: presences not used in estimation
occ.omit <- read.csv("./Wasps/occs_noEurope.csv",header=T)
head(occ.omit)
dim(occ.omit)

# Extract environmental values for each occurrence point
occo1 <- extract(bios1_2,occ.omit[2:3])
occo2 <- na.omit(cbind(occ.omit[2:3],occo1))
colnames(occo2) <- c("long","lat",bio.names)
head(occo2)
dim(occo2)

# Save the resulting matrix
write.csv(occo2,file=paste0("./Wasps/occs_noEurope_",paste(bio.names[1],bio.names[2],sep = "_"),".csv"),row.names = F)

# END #
