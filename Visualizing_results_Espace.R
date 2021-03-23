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
bio1 <- raster("C:\\Users\\l215j162\\Documents\\KU\\Doc-project\\ENM-group\\Wasps\\merra_world\\bio1.asc")
bio12 <- raster("C:\\Users\\l215j162\\Documents\\KU\\Doc-project\\ENM-group\\Wasps\\merra_world\\bio12.asc")
# create a single raster with as meny layers as environmental variables
stck_bios <- stack(bio1, bio12)

# Set parameters particular to the variables and species used
# labels for plot's axis
xy_labs <- c("Annual mean temperature (°C x 10)","Annual Precipitation (mm)") ###

# Read table with MLE of all the species
mle <- c(170.2141,1012.1963,2223.452,8031.546,39554.330)
maha <- c(167.5217,992.9783,2348.477,9813.123,52770.377)

# colors for test occs, Mahalanobis model, weighted model
colpal <- c("darkorchid2","darkolivegreen3","darkorchid4")

# presence points
occpnts <- read.csv("./Wasps/Models-set3/occs_noEurope_50kmthin_Bio1_Bio12.csv", header=T)[3:4]
# invasions
invade <- read.csv("./Wasps/coord_invasion.csv", header=T)[4:5]

# Read M polygon
M.shp <- readOGR("./Wasps/M","M_b4000")

# Get a random sample of points in M and extract its corresponding environmental values
N <- 8000
sampleM <- sam.polyM(M.shp = M.shp,N = N,bios = stck_bios)

# Select the corresponding MLE values from the matrix 'mle.all'
mle.mu <- mle[1:2]
mle.Sig <- matrix(c(mle[3],rep(mle[4],2),mle[5]), nrow = 2, byrow = T)

maha.mu <- maha[1:2]
maha.Sig <- matrix(c(maha[3],rep(maha[4],2),maha[5]), nrow = 2, byrow = T)

# Define the ellipse of the corresponding multivarite normal models
el.mle <- ellipse::ellipse(x=mle.Sig, centre=mle.mu, level=0.99)
el.maha <- ellipse::ellipse(x=maha.Sig, centre=maha.mu, level=0.99)
# el1.mle <- ellipse::ellipse(x=mle.Sig, centre=mle.mu, level=0.95)
# el1.maha <- ellipse::ellipse(x=maha.Sig, centre=maha.mu, level=0.95)
# el2.mle <- ellipse::ellipse(x=mle.Sig, centre=mle.mu, level=0.9)
# el2.maha <- ellipse::ellipse(x=maha.Sig, centre=maha.mu, level=0.9)

# PLOT results for single species ------------------------------
# name for the plot that shows the results
plotname <- "./Wasps/Models-set3/Vespa_Nf.png" ###

# plot will be saved as .png
library(scales)
x11()
#png(filename=plotname, width=15, height=15,
#    units="cm", res=300, pointsize = 6)
par(mar=c(5, 5, 5, 2), xaxs="i", yaxs="i", cex.axis=1.5, cex.lab=1.5,
    cex.main=2)
# plot points in M
plot(sampleM, xlim=c(-60,320), ylim=c(200,2000), col=alpha("grey70",0.3),
     pch=1, xlab=xy_labs[1], ylab=xy_labs[2], main="")
# add presence points
points(occpnts,col=alpha("darkred",0.5),pch=19,cex=1.3)
# add Mahalanobis model
lines(el.maha,col=colpal[3],lwd=3)
# lines(el1.maha,col=colpal[3],lwd=3)
# lines(el2.maha,col=colpal[3],lwd=3)
points(maha.mu[1],maha.mu[2],pch=15,cex=2,col=colpal[3])
# add ellipse - Nf model
lines(el.mle,col="red",lwd=3)
# lines(el1.mle,col="red",lwd=3)
# lines(el2.mle,col="red",lwd=3)
points(mle.mu[1],mle.mu[2],pch=15,cex=2,col="red")
# add invasions
points(invade,col="royalblue3",pch=19,cex=2)
# add legend
elements <- c(expression(italic("Vespa mandarinia")), "Points inside M",
              "Presence points", "Recorded Invasions",
              "Ellipse from weighted-normal model",
              "Ellipse from Mahalanobis model")
legend("topleft", legend = elements, pch=c(NA,1,19,19,NA,NA), cex=1.2,
       col = c("white","black","darkred","royalblue3","red",colpal[3]),
       lwd=c(rep(NA,4),3,3), bty = "n")
#close plot-window
#dev.off()

# other <- read.csv("./Wasps/occs_noEurope_Bio1_Bio12.csv",header = T)[,3:4]
# points(other,col=alpha("darkgreen",0.5),pch=19,cex=1.3)
