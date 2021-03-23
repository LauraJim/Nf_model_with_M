# Laura Jimenez
# First version: May, 2020
# Last modification: November, 2020
# Project: Second model for the fundamental niche, Vespa mandarinia exercise

# Upload packages -----------------------
library(spocc)
library(rgbif)
library(maps)
data("wrld_simpl", package = "maptools")

# Data cleaning -----------------------------

# Read file with occurrence records from GBIF
wasp <- read.csv("./Wasps/occurrences_GBIF.csv",header=T)
head(wasp)
dim(wasp)

# Keep only columns of interest and get rid of NAs
wasp1 <- na.omit(wasp[, c("scientificName", "longitude", "latitude", "year")])
head(wasp1)
dim(wasp1)

# Excluding records older than 1990
wasp2 <- wasp1[wasp1$year >= 1990, -4]
head(wasp2)
dim(wasp2)

# Changing names column
unique(wasp2$scientificName)
wasp2 <- wasp2[wasp2$scientificName == "Vespa mandarinia Smith, 1852", ]
wasp2$scientificName <- "Vespa mandarinia"
dim(wasp2)

# Excluding (0,0) coordinates 
wasp2 <- wasp2[(wasp2$longitude != 0) & (wasp2$latitude != 0), ]
dim(wasp2)

# Excluding duplicate records
wasp2 <- wasp2[!duplicated(wasp2,by=c("longitude","latitude")), ]
dim(wasp2)

# Excluding records without decimal precision
wasp2 <- wasp2[(wasp2$longitude - floor(wasp2$longitude))!=0,]
wasp2 <- wasp2[(wasp2$latitude - floor(wasp2$latitude))!=0,]
dim(wasp2)
write.csv(wasp2,"./Wasps/cleaned_occs.csv",row.names = F)

# Visualization of presences ------------------------
library(sp)
sp_occs <- SpatialPoints(coords=wasp2[,2:3], proj4string=crs(wrld_simpl))

x11()
plot(wrld_simpl)
plot(sp_occs,add=T,col="red",axes=T)

library(raster)
extent_eurasia <- extent(0,180,0,100)
eurasia <- crop(wrld_simpl,extent_eurasia)

# plot cropped polygons and species points
x11()
plot(eurasia,col="gray",axes=T,xaxs="i",yaxs="i",xlim=c(0,180),ylim=c(0,100))
plot(sp_occs,pch=20,col="red",add=T)

# Excluding records in Europe -----------------------
wasp3 <- wasp2[wasp2$longitude > 30, ]
dim(wasp3)
write.csv(wasp3,"./Wasps/occs_noEurope.csv",row.names = F)

# visualization
sp_occs2 <- SpatialPoints(coords=wasp3[,2:3], proj4string=crs(wrld_simpl))
extent_asia <- extent(60,180,0,85)
asia <- crop(wrld_simpl,extent_asia)

# plot cropped polygons and species points
x11()
plot(asia,col="gray",axes=T,xaxs="i",yaxs="i",xlim=c(60,180),ylim=c(0,85))
plot(sp_occs2,pch=20,col="blue",add=T)

# Data thinning
sp_occs3 <- sp::remove.duplicates(sp_occs2, zero = 50) # WGS 84 units are meters
sp_occs3 <- sp::spTransform(sp_occs3,crs(wrld_simpl))@data # to reproject and get data

x11()
plot(asia,col="gray",axes=T,xaxs="i",yaxs="i",xlim=c(60,180),ylim=c(0,85))
plot(sp_occs3,pch=20,col="blue",add=T)

write.csv(sp_occs3@coords,"./Wasps/occs_noEurope_50kmthin.csv",row.names = F)



# END #