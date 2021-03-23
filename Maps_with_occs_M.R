# Laura Jimenez
# First version: November, 2020
# Last modification: November, 2020
###
# Maximum likelihood approach of the fundamental niche estimation problem using a weighted distribution
# where the weights represent the availability of environmental combinations inside M
###

# Packages
library(rgdal)
library(raster)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
#library(ggspatial)

# Creating nice maps 

### Part 1: Vespa mandarinia -----------------------
# Read presence points
sp.occpnts1 <- read.csv("./Wasps/Models-set3/occs_noEurope_50kmthin_Bio1_Bio12.csv",header=T)
head(sp.occpnts1)
dim(sp.occpnts1)
n <- nrow(sp.occpnts1)

occs <- data.frame(Longitude=sp.occpnts1[,1],Latitude=sp.occpnts1[,2])

# Read M polygon
M.shp <- readOGR("./Wasps/M","M_b4000")
M <- spTransform(M.shp, CRS("+proj=longlat +datum=WGS84"))
M <- fortify(M)
Mcol <- rgb(0.35,0,0.2,alpha = 0.6)

# Extent of map
Mext <- extent(M.shp)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# plot vespa -----------------
x11()
ggplot(data = world) +
  geom_sf( ) +
  theme_bw() +
  geom_map(map=M, data=M, aes(map_id=id), color="royalblue3", alpha=0.1) +
  geom_point(data = occs, aes(x = Longitude, y = Latitude), size = 2, 
             shape = 23, fill = "darkred") +
#  annotation_scale(location = "br", width_hint = 0.2) +
#  annotation_north_arrow(location = "br", which_north = "true", 
#                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
#                         style = north_arrow_fancy_orienteering) +
  coord_sf(xlim = c(Mext[1]-10,Mext[2]+10), ylim = c(Mext[3]-10,Mext[4]+10),
           expand = FALSE) +
  theme(panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "aliceblue"))

ggsave('./Wasps/Models-set3/map_vespa_M_occs1.png',  width = 24, height = 12, units = "cm",
       dpi = 600, pointsize = 6)

### Part 2: Hummingbirds -----------------------
spp <- read.csv("./Hummers-Cooper/Analysis3/spp_names.csv",header=T)
occfile <- paste0("./Hummers-Cooper/Analysis3/Data/",spp[,1],"_occ.csv")
spnum <- 6

# Read presence points
spocc <- read.csv(occfile[spnum],header=T)
# transform coords to data frame
occ <- data.frame(Longitude=spocc[,1],Latitude=spocc[,2])

# Read M polygon
M.shp <- readOGR("./Hummers-Cooper/Analysis3/Data",
                 as.character(spp[spnum,1]))
M <- spTransform(M.shp, CRS("+proj=longlat +datum=WGS84"))
M <- fortify(M)

# Colors for polygons
Mcol <- c("chocolate4","firebrick","hotpink","chartreuse4","cadetblue","goldenrod")
# country polygons
world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)

# plot hummingbirds -----------------
x11()
ggplot(data = world) +
  geom_sf( ) +
  geom_map(map=M, data=M, aes(map_id=id), color=Mcol[spnum], size = 1.2,
           alpha=0.1) +
  geom_point(data = occ, aes(x = Longitude, y = Latitude),
             shape = 21, fill = Mcol[spnum], alpha=0.4) +
  coord_sf(xlim = c(-140,-30), ylim = c(-60,60),
           expand = FALSE) +
  theme(panel.grid.major = element_line(color = "white"),
        panel.background = element_rect(fill = "aliceblue"))

ggsave(paste0("./Hummers-Cooper/Analysis3/",as.character(spp[spnum,1]),
              "_map_M_occs.png"),
       width = 12, height = 12, units = "cm", dpi = 600, pointsize = 6)

### Part 3: Suitability maps -----------------------------------
# Vespa
# Read raster with output from weighted model
outp <- raster("./Wasps/Models-set3/Suitability_bios_1_12_M4000.tif")
# emap <- extent(-170, 179, -60, 80) # whole world
# emap <- extent(-140, -110, 30, 65) # USA-CAN
 emap <- extent(70, 150, 10, 55) # Asia
# emap <- extent(-15, 30, 35, 60) # Europe

outp1 <- crop(outp, emap)
outpp <- rasterToPoints(outp1)
outppd <- data.frame(outpp)
colnames(outppd) <- c("Longitude","Latitude","Suitability")

# Read presence points and invasion coords
occ1 <- read.csv("./Wasps/Models-set3/occs_noEurope_50kmthin_Bio1_Bio12.csv",header=T)
occ2 <- read.csv(paste0("./Wasps/coord_invasion.csv"),header=T)
  
x11()
ggplot() +
  geom_tile(data = outppd,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='slateblue1', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=3) +
  coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE) +
  geom_point(data = occ1,aes(x=long, y=lat), shape = 23, fill = "yellowgreen") +
  geom_point(data = occ2,aes(x=lon_init, y=lat_init), shape = 23, fill = "darkred")

ggsave('./Wasps/Models-set3/Suitmap_vespa_Asia.png',  width = 24, height = 12, units = "cm",
       dpi = 600, pointsize = 6)

# Hummingbirds
# Read raster with output from weighted model
outp <- raster("./Hummers-Cooper/Analysis3/thalassinus_suit_map.tif")
 emap <- extent(-170, -20, -60, 80) # America

outp1 <- crop(outp, emap)
outpp <- rasterToPoints(outp1)
outppd <- data.frame(outpp)
colnames(outppd) <- c("Longitude","Latitude","Suitability")

# Read presence points and invasion coords
occ1 <- read.csv("./Hummers-Cooper/Analysis3/Data/Colibri_thalassinus_occ.csv",header=T)

x11()
ggplot() +
  geom_tile(data = outppd,aes(x=Longitude, y=Latitude, fill=Suitability)) +
  theme_bw() +
  #borders("world", xlim = c(-179, 179), ylim = c(-60, 80)) +
  scale_fill_gradient2("Suitability",limits=c(0,1), low = 'grey80',
                       mid='cadetblue', high = 'slateblue4',na.value = NA,
                       midpoint = 0.5, n.breaks=3) +
  coord_sf(xlim = emap[1:2], ylim = emap[3:4], expand = FALSE)


ggsave('./Hummers-Cooper/Analysis3/Colibri_thalassinus_suitmap.png',
       width = 24, height = 12, units = "cm", dpi = 600, pointsize = 6)

## END ##
