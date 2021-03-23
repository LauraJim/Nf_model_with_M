# Laura Jimenez
# First version: March, 2019
# Last modification: November, 2020
# TWO-STAGE SAMPLING FROM A MULTIVARIATE NORMAL DISTRIBUTION AND THE DENSITY OF E_{M}
# 
# Functions --------------------------------------------------------------------------

# Function that determines which points from a given matrix are inside
# an ellipse

## Parameters:
## cloud == matrix that contains the coordinates of the points to evaluate
## centroid == vector with the coordinates of the centroid of the ellipse
## sigma == covariance matrix that defines de ellipse
## alpha == confidence level

el.in <- function(cloud,centroid,sigma,alpha){
  # step 1: calculate de Mahalanobis distance
  maha <- mahalanobis(x=cloud,center=centroid,cov=sigma)
  # step 2: a point is inside the confidence region (1-alpha=confidence%) if
  # the distance divided by the quantile of a Chi-square variable with k d.f. is less than 1
  chi.q <- qchisq(alpha,ncol(cloud))
  check <- (maha/chi.q) <= 1
  cloud.inout <- cbind(rep(1,nrow(cloud))*check,cloud)
  colnames(cloud.inout) <- c("in-out",colnames(cloud))
  return(as.matrix(cloud.inout))
}

# Pre-processing ----------------------------------------------------------------------

# Species fundamental niche --------------
# we are using virtual sp B from Jimenez et al. (2019)
muB <- c(-0.38,1.85)
sigmaB <- matrix(c(0.041227508,0.036665561,0.036665561,0.075), nrow = 2)
elB <- ellipse::ellipse(x=sigmaB,centre=muB,level=0.99)

# Background ---------------
# Existing environmental combinations in a grid covering North America
G <- read.csv("./VirtualSpB-Example/VarsCrdsBackgroundOnly.csv",header=T)
colnames(G)
# Identifying background points inside the fundamental niche
G.inout <- el.in(G[,7:8],muB,sigmaB,0.99)
head(G.inout)
G.inNf <- G[which(G.inout[,1]==1),c(3,4,7,8)]
G.outNf <- G[which(G.inout[,1]==0),c(3,4,7,8)]

# Points inside the M  ------------
## M was delimited using ecoregions and a sample of spB
pntsM <- read.csv("./VirtualSpB-Example/spB_points_inside_ecoregions_M.csv",header = T)
head(pntsM)
dim(pntsM)
# check which points are inside Nf
M.inout <- el.in(pntsM[,3:4],muB,sigmaB,0.99)
head(M.inout)
M.inNf <- pntsM[which(M.inout[,1]==1),]
M.outNf <- pntsM[which(M.inout[,1]==0),]

# Points from background and M that are inside the Nf (existing niche, Ne)
Ne <- rbind(G.inNf,M.inNf)
head(Ne)
# Points in M 
# Adding a column in NApnts.spBNf that indicates if a point is inside de M (==1)
GinM <- as.vector(G.inNf[,1] %in% pntsM[,1] *1)
# Separate spB.inMNf into two groups
G.inNfinM <- G.inNf[GinM==1,]
G.inNfoutM <- G.inNf[GinM==0,]

# Color palette ---------------------------
library(scales)
# background
col.b <- alpha("grey70",0.8)
p.b <- 1
# points inside both Nf and M
col.in <- alpha("darkmagenta",0.8)
# points inside Nf but outside M
col.out <- alpha("green4",0.8)
# points inside M and outside Nf
col.M <- alpha("sienna2",0.6)
p.m <- 16
p.om <- 19

# Plot different subsets of points ----------------
#x11()
png('./VirtualSpB-Example/Figure1_4600x2300px600dpi4pnt.png',
    width = 4600, height = 2300, res = 600, pointsize = 4)
par(mfrow=c(1,2),cex.lab=1.6,cex.main=1.8)
### ENVIRONMENTAL SPACE
# background
plot(G.outNf$Bio1WHStnd,G.outNf$Bio12WHStnd,pch=p.b,cex=0.8,col=col.b,
     xlab="Bio1 (standardized)",ylab="Bio12 (standardized)",main="Environmental space")
# env. combinations inside M but outside the niche
points(M.outNf$Bio1WHStnd,M.outNf$Bio12WHStnd,col=col.M,pch=p.m)
# env. combinations inside M and inside the niche
points(M.inNf$Bio1WHStnd,M.inNf$Bio12WHStnd,col=col.in,pch=p.om)
# env. combinations outside M but inside the niche
points(G.inNfoutM$Bio1WHStnd,G.inNfoutM$Bio12WHStnd,col=col.out,pch=p.om)
# fundamental niche
lines(elB,lwd=2)
legend(x=-3,y=7.5,legend = c("Fundamental niche","Existing environmental conditions in the Americas"),
       lwd=c(2,NA),pch = c(NA,19),col=c("black",col.b),bty = "n",cex = 1.6)
### GEOGRAPHIC SPACE
# background
plot(G.outNf$Long,G.outNf$Lat,pch=p.b,cex=0.8,col=col.b,xlab="Longitude",
     ylab="Latitude",main="Geographic space")
# sites inside M but outside the niche
points(M.outNf$Long,M.outNf$Lat,col=col.M,pch=p.m)
# sites inside M and inside the niche
points(M.inNf$Long,M.inNf$Lat,col=col.in,pch=p.om)
# sites outside M but inside the niche
ind <- which(G.inNfoutM$Long > -100) # eliminate red stars in Nf that were not identified in M
points(G.inNfoutM$Long[ind],G.inNfoutM$Lat[ind],col=col.out,pch=p.om)
## Explain each type of point
legend(x=-180,y=-30,legend = c("Sample of sites inside the Americas","Inside of both the fundamental niche and M",
                               "Inside M but outside the fundamental niche","Outside of M but inside the fundamental niche"),
       pch=rep(19,4),col=c(col.b,col.in,col.M,col.out),bty = "n",cex = 1.6)
dev.off()

# STAGE 1: sample of size N from mvtnorm --------------------------------------------------------
N <- 100
# generate a sample from the corresponding mntnormal 
mtv.sam <- tmvtnorm::rtmvnorm(n=N, mean=muB, sigma=sigmaB)
colnames(mtv.sam) <- c("Bio1-std","Bio12-std")

# STAGE 2: subsample of size n using weights ----------------------------------------------------
# Select points from original sample according to the values from the kernel
n <- 25
subsam <- sample(1:N,n,prob = ws)
bias.sam <- mtv.sam[subsam,]

# Calculate and draw the kernel using the points inside M (define M carefully)
# kernel for contour plot
fhat.M <- ks::kde(x=cbind(pntsM$Bio1WHStnd,pntsM$Bio12WHStnd))
# kernel evaluated in the sample
fhat.Mval <- ks::kde(x=cbind(pntsM$Bio1WHStnd,pntsM$Bio12WHStnd),eval.points=mtv.sam)
# weights
ws <- as.vector(fhat.Mval$estimate)
ws <- ws/sum(ws)
#x11()
#hist(ws)

# Save samples
#write.csv(mtv.sam,"./VirtualSpB-Example/mvt_normal_sample_100.csv",row.names = F)
#write.csv(bias.sam,"./VirtualSpB-Example/mvt_normal_subsample_25.csv",row.names = F)
# or
mtv.sam <- read.csv("./VirtualSpB-Example/mvt_normal_sample_100C.csv",header = T)
bias.sam <- read.csv("./VirtualSpB-Example/mvt_normal_subsample_25C.csv",header = T)

# PLOTTING DIFFERENT SAMPLE SCHEMES -----------------------------------------------------
# define the plot limits from the range of the points in M
xlims <- c(round(min(pntsM$Bio1WHStnd)-0.2,1),round(max(pntsM$Bio1WHStnd)+0.2,1))
ylims <- c(round(min(pntsM$Bio12WHStnd)-0.2,1),round(max(pntsM$Bio12WHStnd)+0.2,1))

# Environmental space ---------------------------

#x11()
png('./VirtualSpB-Example/Figure2_7500x2500px600dpi10pnt.png',
    width = 7500, height = 2500, res = 600, pointsize = 10)
par(mfrow=c(1,3),cex.lab=1.5,xaxs="i",yaxs="i")
### PLOT 1: Probability of selecting a point, given the kernel of M
lvls2 <- c(0.25,0.5,0.75,0.9,0.95,0.99)
Nf.cols <- colorRampPalette(c("aquamarine4","aquamarine"))
Nfc <- Nf.cols(length(lvls2))
# plot multivariate normal contours
plot(xlims,ylims,type="n",xlab="",ylab="Bio12 (standardized)",main="")
for (j in length(lvls2):1) {
  car::dataEllipse(mtv.sam[,1],mtv.sam[,2],levels=lvls2[j],grid=F,col=Nfc[j],lwd=2,
                   plot.points = F,add = T,fill = T,fill.alpha=0.2,ellipse.label = lvls2[j])
}
# add points used for kernel estimation
points(mtv.sam,col=col.out,pch=17,cex=1.2)
points(colMeans(mtv.sam)[1],colMeans(mtv.sam)[2],pch=19,cex=2)
### PLOT 2: Probability of selecting a point, given the fundamental niche
# plot estimated kernel of M with points inside M
lvls1 <- c(10,25,50,75,95,99)
M.cols <- colorRampPalette(c("white",col.M))
plot(fhat.M,display="filled.contour",cont=lvls1,main="",xlab="Bio1 (standardized)",
     ylab="",xlim=xlims,ylim=ylims,col=M.cols(length(lvls1)+1))
# add points used for kernel estimation
points(pntsM$Bio1WHStnd,pntsM$Bio12WHStnd,col=col.M,pch=19)
# points from M inside fundamental niche
points(M.inNf$Bio1WHStnd,M.inNf$Bio12WHStnd,col=col.in,pch=19)
### PLOT 3: Sample from normal distribution vs. sample from two-stage model
col.sub <- col.in
plot(xlims,ylims,type="n",xlab="",ylab="", main="")
# points from M outside fundamental niche
points(M.outNf$Bio1WHStnd,M.outNf$Bio12WHStnd,col=col.M,pch=19)
# add ellipse that represents the Nf
car::dataEllipse(mtv.sam[,1],mtv.sam[,2],levels=0.99,grid=F,col=Nfc[length(lvls2)],lwd=2,
                 plot.points = F,fill = T,fill.alpha=0.2,ellipse.label = 0.99,add=T)
# add ellipse estimated from subsample
car::dataEllipse(bias.sam[,1],bias.sam[,2],levels=0.99,grid=F,col=col.sub,lwd=2,
                 plot.points = F,fill = T,fill.alpha=0.15,ellipse.label = 0.99,add=T)
# # add whole sample from multivariate normal distribution
# points(mtv.sam,col="brown",pch=17,cex=1.2)
# add subsample selected with the weights calculated from the kernel of M
points(bias.sam,col=col.sub,pch=17,cex=1.5)
# emphasize the centers of the ellipses
points(colMeans(mtv.sam)[1],colMeans(mtv.sam)[2],pch=19,cex=2)
points(colMeans(bias.sam)[1],colMeans(bias.sam)[2],pch=19,cex=2)
dev.off()

# Geographical space -------------------
x11()
# background
plot(env$Long,env$Lat,pch=".",xlim=c(-180,-100),ylim=c(20,80),xlab="Longitude",ylab="Latitude",main="Geographical space")
# sites inside M (the one selected with ecoregions)
points(pntsM$Long,pntsM$Lat,col=colB,pch=8)
# sites inside existing (i.e. inside fudamental niche and rectaular M)
points(spBinM$x,spBinM$y,col="steelblue2",pch=8)
# select points from blue set using the kernel
# or, extract values from Bio1 and Bio12 for the set of points (x[subsam,1],x[subsam,2])

# Defining M by using Ecoregions ------------------------------------------------------
# library(raster)
# library(rgdal)
# library(maptools)
# library(sp)
# data("wrld_simpl")
# # we'll use ecoregions to define M
# eco.shp <- readOGR(".\\2019\\na_cec_eco_l2","NA_CEC_Eco_Level2")
# extent(eco.shp)
# # transform the ecoregion's projection to datum WGS84
# eco.tr <- sp::spTransform(eco.shp,crs(wrld_simpl))
# extent(eco.tr)
# # plot all the ecoregions in eco.shp with the points of interest on top
# x11()
# plot(eco.tr)
# points(occ_sp,pch=19,col=3)
# # make a spatial object from coordinates of points inside Nf
# occ_sp <- SpatialPoints(coords = NApnts.spBNf,proj4string = crs(wrld_simpl))
# # select ecoregions (polygons) that overlap with points
# eco_occsp <- eco.tr[!is.na(over(eco.tr,occ_sp)), ]
# 
# # new plot with selected ecoregions and sites
# x11()
# plot(eco_occsp)
# # now, find all the points in North America that are inside of the selected ecoregions
# # select points inside the same extent of the ecoregions shapefile
# env_sp <- SpatialPoints(coords = env[,3:4],proj4string = crs(wrld_simpl))
# env_sp1 <- crop(env_sp,extent(eco.tr))
# # select points inside the ecoregions identified to be in M
# eco_pnts <- env_sp1[eco_occsp]
# points(eco_pnts,pch=19,col=2)
# points(occ_sp,pch=19,col=3)
# M_pnts <- as.data.frame(eco_pnts)
# # adding the corresponding environmental measurements
# mns <- colMeans(env[,5:6])
# sds <- c(sd(env[,5]),sd(env[,6]))
# bio1 <- raster(".\\ClimateData10min\\bio1WH.asc")
# bio12 <- raster(".\\ClimateData10min\\bio12WH.asc")
# st <- stack(bio1,bio12)
# stz <- scale(st,center=mns,scale=sds) # careful with the values used for scaling!!!
# clim_pnts <- extract(stz,M_pnts)
# write.csv(cbind(M_pnts,clim_pnts),".\\2019\\spB_points_inside_ecoregions_M.csv",row.names = F)

# # Special case: Truncated multivariate normal distribution --------------
# library(tmvtnorm)
# up <- c(max(env.sp[,1]),max(env.sp[,2]))
# low <- c(min(env.sp[,1]),min(env.sp[,2]))
# ind <- which(env[,6] > low[1] & env[,6] < up[1] & env[,7] > low[2] & env[,7] < up[2])
# env.box <- env[ind,]
# x <- rtmvnorm(n=50, mean=mus.fn[spn,], sigma=sigma.fn, upper=up,lower=low)
# x11()
# plot(env.box[,6:7], main="samples from truncated bivariate normal distribution",xlim=c(0,1.5), ylim=c(-1,4),
#      xlab=expression(x[1]), ylab=expression(x[2]), pch=19)
# points(x,pch=19,col=el.col[spn])
# lines(el.sp,col=el.col[spn],lwd=2)
