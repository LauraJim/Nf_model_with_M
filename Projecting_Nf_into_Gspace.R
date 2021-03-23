# Laura Jimenez
# First version: June 2020
# Last version: December 2020
# Project resulting ellipses back into G-space

setwd("C:\\Users\\l215j162\\Documents\\KU\\Doc-project\\NicheEstimation")
library(raster)

# Read environmental layers cropped to the area of interest
bio1.proj <- raster(".\\ClimateData10min\\bio1WH.asc")
bio2.proj <- raster(".\\ClimateData10min\\bio12WH.asc")
stck_1_2.proj <- stack(bio1.proj,bio2.proj)

# Set the values of the MLEs
mle.mu <- c(177.4347015,1640.714552)
mle.A <- matrix(c(0.000650966,-1.37E-05,-1.37E-05,1.83E-06), nrow = 2, byrow = T)
mle.Sig <- chol2inv(chol(mle.A))

# Calculate suitabilities in each cell
max.val <- mvtnorm::dmvnorm(x=mle.mu,mean = mle.mu, sigma = mle.Sig)
# function that calculates the log(suitability)
sui.fun <- function(cell){log(mvtnorm::dmvnorm(x=c(cell[1],cell[2]),mean = mle.mu, sigma = mle.Sig))-log(max.val)}
# apply this function to the whole raster layer
suit.rast <- calc(stck_1_2.proj,fun=sui.fun)
# take exponential to go back to the original scale
suit.rast1 <- calc(suit.rast,fun = exp,
                   filename="C:\\Users\\l215j162\\Dropbox\\Nf-second model\\Hummers-Cooper\\Analysis3\\thalassinus_suit_map_maha.asc",
                   overwrite=T)

# save a TIFF
writeRaster(suit.rast1,"C:\\Users\\l215j162\\Dropbox\\Nf-second model\\Hummers-Cooper\\Analysis3\\thalassinus_suit_map_maha.tiff", overwrite = T)

x11()
plot(suit.rast1)
#points(occ.geo,pch=15,col="red",cex=0.5)

# END