# Laura Jimenez / February 2020
# Maximum likelihood approach of the fundamental niche estimation problem using a weighted distribution
# where the weights represent the availability of environmental combinations inside M
# 
# Setting working directory
setwd("C:/Users/l215j162/Documents/KU/Doc-project/NicheEstimation")
# Calling packages
#library(sp)
library(raster)
library(rgdal)
library(rgeos)
#library(maptools)
#data("wrld_simpl")

# FUNCTIONS -------------------------------------------------------------

# Get a random sample of points inside the polygon that delimits M and extract their environmental values
## M.shp should be a shapefile containing the polygon that defines the area of study
## N sample size
## bios must be a stack of environmental raster layers
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
  # set coordinate system from world map to M
  #M.tr <- sp::spTransform(M.shp,crs(wrld_simpl))
  # clip the continent with extent of polygon
  #clip <- gIntersection(wrld_simpl, M.tr, byid = TRUE, drop_lower_td = TRUE) #clip polygon 2 with polygon 1
  # get sample of points inside intersection
  #sam <- spsample(clip, n = N, "random")
  #clim_M <- extract(stck_1_12,M_pnts)
}

# Negative log-likelihood function for theta=(mu,A)
## guess -- is a vector of length 5 when d=2, it contains the mu and A values as elements
## sam1 -- matrix containing the original sample of environmental combinations that correspond to presences
## sam2 -- matrix containing a second random sample of environmental combinations which come from the area of study (M)
negloglike <- function(guess,sam1,sam2){
  # define the parameters of interest from the guess parameter
  mu <- guess[1:2]
  A <- matrix(c(guess[3],guess[4],guess[4],guess[5]),nrow=2,ncol=2)
  # original sample size: number of presence points in the sample
  n <- nrow(sam1)
  # function that calculates quadratic terms
  quad <- function(xi) { ax<-as.matrix(xi - mu); t(ax) %*% A %*% ax }
  q1 <- apply(sam1, 1, quad) # quadratic terms of presence points
  q2 <- apply(sam2, 1, quad) # quadratic terms of M points
  # negative log-likelihood value
  S <- 0.5*sum(q1) + n*log(sum(exp(-0.5*q2)))
  return(S)
}

# MAIN --------------------------------------------------------------------------------

# Read datasets ---------------------------------------------------------------
ext1 <- "./2019-20/Nf-second model/Hummers-Cooper/Analysis3/"
# Read environmental layers
bio1 <- raster("./ClimateData10min/bio1WH.asc")
bio12 <- raster("./ClimateData10min/bio12WH.asc")
stck_1_12 <- raster::stack(bio1,bio12)

# Read species names
spp.nms <- spp <- read.csv(paste0(ext1,"spp_names.csv"),header=T)
spnames <- as.vector(spp.nms[,1],mode = "character")

# Colors for different species
Mcol <- c("chocolate4","firebrick","hotpink","chartreuse4","cadetblue","goldenrod")
# colors for test occs, Mahalanobis model, weighted model
colpal <- c("darkorchid2","darkolivegreen3","darkorchid4")
  
# Set parameters --------------------------------------------------------------
# Set the sample size (>= 10,000) that will be used to approximate expected value
N <- 10000

# Choose the species to work with
i <- 1

# Perform estimation for all the species --------------------------------------
  # Read presence points and M polygon
  sp.occpnts0 <- read.csv(paste0(ext1,"Data/",spnames[i],"_occ.csv"),header=T)
  sp.occpnts1 <- na.omit(sp.occpnts0[,3:4])
  dim(sp.occpnts1)
  n <- nrow(sp.occpnts1)
  M.shp <- readOGR(paste0(ext1,"Data"),spnames[i])
  
  # Get a subsample of presence points to evaluate model robustness
  # # 1) set the proportion of occurrence points that will be used for model estimation
  # per <- 0.80 # change name of file 
  # # 2) calculate the amount of points to be included in the sample
  # nsub <- floor(n*per)
  # # 3) get a subsample of occurrence points
  # ind.sub <- sort(sample(x=1:n,size=nsub,replace=F))
  # # 4) proceed with the estimation using the new subsample
  # sp.occpnts <- sp.occpnts1[ind.sub,]
  # dim(sp.occpnts)
  # write.csv(sp.occpnts0[-ind.sub,],
  #           paste0(ext1,"test-occ/",spnames[i],"_test.csv"),row.names = F)
  sp.occpnts <- sp.occpnts1
  
  # get a random sample of points in M and extract its corresponding environmental values
  sam.Mpnts <- sam.polyM(M.shp = M.shp,N = N,bios = stck_1_12)
  
  # Lookig for the MLE of mu and A --------------------------------------------
  # get initial values for parameters
  # for mu:
  mu.ini <- colMeans(sp.occpnts)
  # for A:
  Sig.ini <- cov(sp.occpnts)
  A.ini <- chol2inv(chol(Sig.ini))
  # whole vector of inicial values
  vals.ini <- c(mu.ini,A.ini[1,1],A.ini[1,2],A.ini[2,2])#c(mu.ini,A.ini[1,1],A.ini[1,2],A.ini[2,2])
  #negloglike(vals.ini,sp.occpnts,sam.Mpnts)
  
  # fix the values of the samples used to evaluate the neg-log-likelihood
  like.fn <- function(theta){ negloglike(theta,sam1=sp.occpnts,sam2=sam.Mpnts) }
  # apply optimization function to minimize the negative log-likelihood
  #find.mle <- optim(par=vals.ini, fn=like.fn, method="L-BFGS-B", lower=c(-50,0.001,0.001,-999999,0.001),
  #                  upper=c(50,10000,rep(-999999,3)))
  find.mle <- optim(par=vals.ini, fn=like.fn, method="Nelder-Mead")
  mle <- find.mle$par
  mle.mu <- mle[1:2]
  mle.A <- matrix(c(mle[3:4],mle[4:5]),nrow=2,ncol=2)
  mle.Sig <- chol2inv(chol(mle.A))
  #mle.Sig <- solve(mle.A)
  
  # get the ellipse defined by the ml estimators
  el<-ellipse::ellipse(x=mle.Sig, centre=mle.mu, level=0.99)
  
  # get the ellipse from a multivarite normal model / mahalanobis distance method
  el.ml <- ellipse::ellipse(x=Sig.ini, centre=mu.ini, level=0.99)
  
  # plot will be saved as .png
  png(paste0(ext1,spnames[i],"_mle.png"),width = 800, height = 800)
  # plot points inside M
  plot(sam.Mpnts,col="grey70",pch=1,xlim=c(0,350),ylim=c(0,8200),
       xlab="Annual mean temperature (°C*10)", ylab="Annual precipitation (mm)")
  # add presence points to the plot
  points(sp.occpnts1[-ind.sub,],col=colpal[1],pch=19,cex=1.5) # test presences
  points(sp.occpnts,col=Mcol[i],pch=19,cex=1.5) # presences used in model
  # ellipse maha
  lines(el.ml,col=colpal[2],lwd=2)
  # ellipse mle
  lines(el,col=colpal[3],lwd=2)
  sp.leg <- paste(spnames[i],"(",nrow(sp.occpnts),")")
  legend("topleft",legend = c(sp.leg,"Points inside M","Presences",
                              "Testing presences",
                              "Ellipse from Mahalanobis method",
                              "Ellipse from weighted-normal model"),
         pch=c(NA,1,19,19,NA,NA),col = c("white","black",Mcol[i],colpal),
         lwd=c(NA,NA,NA,NA,2,2),bty = "n")
  # close plot-window
  dev.off()
  
  # save MLE values for both models, Weighted-Normal and Mahalanobis
  (MLE.vec <- c(spnames[i],mle,vals.ini))

# Save MLE matrix --------------------------------------------
#MLE.summary <- MLE.vec
(MLE.summary <- rbind(MLE.summary,MLE.vec))
colnames(MLE.summary) <- c("SpName","mle.mu1","mle.mu2","mle.A11","mle.A12","mle.A22",
                        "maha.mu1","maha.mu2","maha.A11","maha.A12","maha.A22")
head(MLE.summary)
write.csv(MLE.summary,paste0(ext1,"Estimated_parameters_Hummers.csv"),row.names = F)

# END