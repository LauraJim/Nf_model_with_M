# Laura Jimenez
# First version: August 2020
# Last modification: November, 2020

# Calling packages
library(raster)
library(rgdal)
library(rgeos)
library(tools)

# FUNCTIONS -------------------------------------------------------------

# Get a random sample of points inside the polygon that delimits M and extract their environmental values
## M.shp -- should be a shapefile containing the polygon that defines the area of study
## N -- sample size
## bios -- must be a stack of environmental raster layers
sam.polyM <- function(M.shp,N,bios){
  # crop and mask the environmental layers with the M polygon
  clip.M <- mask(crop(bios,M.shp),M.shp)
  # get ride of cells with NA values
  ind <- which(!is.na(clip.M[[1]][]))
  # get a random sample of indices
  sam <- sample(ind,N,replace = T)
  # choose the points corresponding to the selected indices
  Mpnts <- clip.M[][sam,]
  return(Mpnts)
}

# MAIN --------------------------------------------------------------------------------

# Read environmental layers
bio1 <- raster("C:/Users/l215j162/Documents/KU/Doc-project/NicheEstimation/ClimateData10min/bio1WH.asc")
bio12 <- raster("C:/Users/l215j162/Documents/KU/Doc-project/NicheEstimation/ClimateData10min/bio12WH.asc")
# create a single raster with as meny layers as environmental variables
stck_bios <- stack(bio1, bio12)

# Set parameters particular to the variables and species used
# labels for plot's axis
xy_labs <- c("Annual mean temperature (°C x 10)","Annual Precipitation (mm)") ###
# directory where the results will be saved
ext1 <- "./Hummers-Cooper/Analysis3/" ###
# name for file that contains all the MLEs mvalues
mle.name <- "Estimated_parameters_Hummers.csv" ###

# Read species names
spp.nms <- spp <- read.csv(paste0(ext1,"spp_names.csv"),header=T)
spnames <- as.vector(spp.nms[,1],mode = "character")
nsp <- length(spnames)

# Read table with MLE of all the species
mle.all <- read.csv(paste0(ext1,mle.name),header=T)
head(mle.all)

# Colors for different species
Mcol <- c("chocolate4","firebrick","hotpink","chartreuse4","cadetblue","goldenrod")
# colors for test occs, Mahalanobis model, weighted model
colpal <- c("darkorchid2","darkolivegreen3","darkorchid4")


# list with occurrence matrices
locc <- vector("list", nsp)
# list with random samples from Ms
lM <- vector("list", nsp)
# list with ellipses
lel <- vector("list", nsp)
# list with Mahalanobis ellipses
lmaha <- vector("list", nsp)
for (i in 1:nsp) {
  # Read fitting and testing points
  # presence points
  occpnts <- read.csv(paste0(ext1,"Data/",spnames[i],"_occ.csv"),
                      header=T)[3:4]
  locc[[i]] <- occpnts

  # Read M polygon
  M.shp <- readOGR(paste0(ext1,"Data"),spnames[i])

  # Get a random sample of points in M and extract its corresponding environmental values
  N <- 8000
  sampleM <- sam.polyM(M.shp = M.shp,N = N,bios = stck_bios)
  lM[[i]] <- sampleM

  # Select the corresponding MLE values from the matrix 'mle.all'
  mle.mu <- c(mle.all$mle.mu1[i],mle.all$mle.mu2[i])
  mle.A <- matrix(c(mle.all$mle.A11[i],rep(mle.all$mle.A12[i],2),
                    mle.all$mle.A22[i]),nrow = 2, byrow = T)
  mle.Sig <- chol2inv(chol(mle.A))
  
  maha.mu <- c(mle.all$maha.mu1[i],mle.all$maha.mu2[i])
  maha.A <- matrix(c(mle.all$maha.A11[i],rep(mle.all$maha.A12[i],2),mle.all$maha.A22[i]),
                  nrow = 2, byrow = T)
  maha.Sig <- chol2inv(chol(maha.A))
  
  # Define the ellipse of the corresponding multivarite normal models
  lel[[i]] <- ellipse::ellipse(x=mle.Sig, centre=mle.mu, level=0.99)
  lmaha[[i]] <- ellipse::ellipse(x=maha.Sig, centre=maha.mu, level=0.99)
}

# PLOT results for all species ------------------------------
library(scales)
# name for the plot that shows the results
plotname <- paste0(ext1,"Nf_allspp_occs2.png") ###

# plot will be saved as .png
#x11()
png(filename=plotname, width=15, height=15,
    units="cm", res=300, pointsize = 6)
par(mar=c(5, 5, 5, 2), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5,
    cex.main=2)
# plot points in M
plot(lM[[1]], xlim=c(-20,350), ylim=c(-1000,8000), type="n",
     xlab=xy_labs[1], ylab=xy_labs[2], main="Environmental space")
abline(h=0,col="grey70")
# for (j in 1:nsp) {
#   points(lM[[j]],col=alpha("grey70",0.3),pch=1)
# }
for (l in c(2,4,5,3,6,1)) {
  # add presence points
  points(locc[[l]],col=alpha(Mcol[l],0.4),pch=19,cex=1.3)
}
for (k in c(1,2,4,5,3,6)) {
  # add ellipse - Nf model
  lines(lel[[k]],col=Mcol[k],lwd=4)
  points(mle.all[k,2],mle.all[k,3],pch=19,cex=2)
}
# add legend
elements <- c("Estimated fundamental niches",gsub("_"," ",spnames),
              "Center of the niche")
legend("topleft", legend = elements, pch=c(rep(NA,nsp+1),19), cex=1.2,
       col = c("white",Mcol,"black"), lwd=c(NA,rep(3,nsp),NA), bty = "n")
#close plot-window
dev.off()

# PLOT results for single species ------------------------------
# select species
sp <- 1
# name for the plot that shows the results
plotname <- paste0(ext1,spnames[sp],"_Nf.png") ###

# plot will be saved as .png
#x11()
png(filename=plotname, width=15, height=15,
    units="cm", res=300, pointsize = 6)
par(mar=c(5, 5, 5, 2), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5,
    cex.main=2)
# plot points in M
plot(lM[[sp]], xlim=c(-50,350), ylim=c(0,7800), col=alpha("grey70",0.3),
     pch=1, xlab=xy_labs[1], ylab=xy_labs[2], main="")
# add presence points
points(locc[[sp]],col=alpha(Mcol[sp],0.4),pch=19,cex=1.3)
# add Mahalanobis model
lines(lmaha[[sp]],col=colpal[3],lwd=3)
points(mle.all[sp,7],mle.all[sp,8],pch=15,cex=2.5,col=colpal[3])
# add ellipse - Nf model
lines(lel[[sp]],col=Mcol[sp],lwd=3)
points(mle.all[sp,2],mle.all[sp,3],pch=15,cex=2.5)
# # add legend
# elements <- c("Estimated fundamental niches",gsub("_"," ",spnames[sp]),
#               "Center of the niche")
# legend("topleft", legend = elements, pch=c(NA,NA,19), cex=1.2,
#        col = c("white",Mcol[sp],"black"), lwd=c(NA,2,NA), bty = "n")
#close plot-window
dev.off()


# END #
