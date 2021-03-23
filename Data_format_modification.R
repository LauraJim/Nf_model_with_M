# Laura Jimenez; November, 2019

# Exploring Humminbirds data
setwd("C:/Users/l215j162/Documents/KU/Doc-project/NicheEstimation")

# Read packages and world map
library(raster)
library(rgdal)
library(maptools)
data("wrld_simpl")

# Read background data ----------------------------------------------
# These file contains 20K points in the Americas
env <- read.csv("./2017-18/VarsCrdsBackground.csv",header=T)
dim(env)
head(env)
x11()
plot(env$x,env$y,pch=".")

# Convert points into SpatialPoints
env_M <- SpatialPoints(coords = env[,2:3],proj4string = crs(wrld_simpl))

# Read environmental layers ---------------------------------------------
bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio12 <- raster("./ClimateData10min/bio12WH.asc")
st1_12 <- stack(bio1,bio12)

# Get list of names of all the species -------------------------------------
# Method 1:
filelst <- list.files("./2019/Hummers-Cooper/occurrences")
spnames <- vector("character",length(filelst))
# Method 2:
# Read species names
spp.nms <- read.csv("./2019-20/Hummers-Cooper/AnalysisTwo/occs-300-1000/spnames.csv",header = T)
spnames <- as.vector(spp.nms[,1],mode = "character")

# Use shapefiles to select points inside M ----------------------------------
# Also, add occurrence points and extract their environmental combinations
for(i in 1:length(spnames)){
  # getting rid of the file extension
  #spnames[i] <- as.character(strsplit(filelst[i],"\\.")[[1]])
  # read the shapefile that contains the M polygon
  M.shp <- readOGR("./2019-20/Hummers-Cooper/AnalysisTwo/M-Hummers",spnames[i])
  # set coordinate system from world map to M
  M.tr <- sp::spTransform(M.shp,crs(wrld_simpl))
  # crop the continent with extent of polygon
  env_M1 <- crop(env_M,extent(M.tr))
  # select and save the points from 'env_M1' that are inside the polygon
  pts.in.M <- env_M1[M.tr,]
  M_pnts <- as.data.frame(pts.in.M)
  names(M_pnts) <- c("long","lat")
  # read occurrence points and extract their environmental combinations
  sp.occ <- read.csv(paste0("./2019-20/Hummers-Cooper/AnalysisTwo/occs-300-1000/",spnames[i],".csv"),header=T)[2:3]
  # adding the corresponding environmental measurements
  clim_M <- extract(st1_12,M_pnts) # M
  clim_occ <- extract(st1_12,sp.occ) # occurrences
  # plot points in the whole continent and save image
  png(paste0("./2019-20/Hummers-Cooper/AnalysisTwo/XYBio1Bio12/",spnames[i],".png"),width = 1200,height = 700)
  par(mfrow=c(1,2)) # two plots
    # Plot 1: G-space
    plot(env_M,pch=".",main=spnames[i],)
    # add polygon
    plot(M.tr,add=T,col=rgb(0.35,0,0.2,alpha = 0.6))
    # add points inside M
    points(pts.in.M,pch=20,col="slateblue4",cex=0.5)
    # add occurrence points
    points(sp.occ,pch=20,col="yellow",cex=0.6)
    # Plot 2: make equivalent plot in E-space
    plot(env$Bio1WH,env$Bio12WH,pch=".",col="gray",main="Environmental space",xlab="Bio1",ylab="Bio12")
    # add points inside M
    points(clim_M,pch=20,col="slateblue4",cex=0.5)
    # add occurrence points
    points(clim_occ,pch=8,col="orange",cex=0.6)
    legend("topleft",legend = c("Background","Points inside M","Occurrences"),pch = c(19,19,8),
           col = c("gray","slateblue4","orange")) #,bg="gray"
    # close plot-window
    dev.off()
  # save table with geographic coordinates and climatic combinations of M
  write.csv(cbind(M_pnts,clim_M),paste0("./2019-20/Hummers-Cooper/AnalysisTwo/XYBio1Bio12/",spnames[i],"_M.csv"),row.names = F)
  # save table with geographic coordinates and climatic combinations of occurrence points
  write.csv(cbind(sp.occ,clim_occ),paste0("./2019-20/Hummers-Cooper/AnalysisTwo/XYBio1Bio12/",spnames[i],"_occ.csv"),row.names = F)
}


### END