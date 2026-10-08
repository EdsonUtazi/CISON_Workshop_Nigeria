set.seed(500)


####################################################################################
library(geosphere)
library(fields)
library(gstat)
library(sf)
library(terra)

#Set working directory
setwd("path-to-directory")

#loading the data
dat      <- read.csv("Practice_data_outcome.csv", header=TRUE)
dat.cov  <- read.csv("Covariates_selected.csv", header=TRUE)

data.merge <- merge(dat, dat.cov[,1:11], by="ID")
coords <- data.frame(LONGNUM = data.merge$LONGNUM, LATNUM = data.merge$LATNUM)

#####Some distance calculations - Place first

#Euclidean distances in degrees. To convert to km, multiply by 111. 
#This is because 1 degree longitude = 111.32 km at the equator
#(the equator is divided into 360 degrees of longitude).
#However, as you move towards the poles, this distance decreases.
hh <- dist(coords)
dist1 <- as.numeric(dist(coords)*111)
summary(dist1)
max.dist <- max(dist(coords))  
min.dist <- min(dist(coords))

##Great circle or spherical distance or geodetic distance in meters, 
#geosphere package. Divide by 1000 to convert to km.
#For Nigeria, this gives results similar to the Euclidean distances.
dd <- distm(coords, fun=distHaversine)
dist2 <- as.numeric(dd[lower.tri(dd, diag = FALSE)])/1000
summary(dist2)


#The Earth's shape is often approximated using an oblate ellipsoid model that bulges at the equator and flattens at
#the poles. The most common ellipsoid is the World Geodetic System (WGS84) used by GPS devices. 
#A coordinate reference system (CRS) specifies how coordinates are related to locations on the Earth.

#Convert the distances from lon-lat to UTM coordinates

# Create spatial points from longitude and latitude
SP_longlat <- st_as_sf(data.merge, coords = c("LONGNUM", "LATNUM"), crs = 4326)
#EPSG:4326 is the the standard identifier for the WGS 84 CRS

# Transform to UTM Zone 32N
SP_UTM <- st_transform(SP_longlat, crs = 32632) #Nigeria is in UTM zone 32N.
#EPSG:32632 is specifically WGS 84 / UTM Zone 32N, covering longitudes 6°E–12°E

#View both lon-lat and UTM coordinates
head(st_coordinates(SP_longlat))
head(st_coordinates(SP_UTM))

#Calculate Euclidean distances in meters between the UTM coordinates
#and divide this by 1000 to convert to km
dist_UTM <- dist(st_coordinates(SP_UTM))/1000

#Observe that this gives very similar results as dist2
summary(dist_UTM)

#Some plots
#Compare Euclidean distance and great circle distance, both calculated using lon-lat coordinates
plot(dist1, dist2)

#Compare great circle distance calculated using lon-lat coordinates and Euclidean distance calculated 
#using the UTM coordinates
plot(dist2, dist_UTM)

#Plot the lon-lat coordinates and the UTM coordinates
par(mfrow=c(1,3))
plot(st_coordinates(SP_longlat), main = "lon-lat")
plot(st_coordinates(SP_UTM), main = "UTM")
###################################################################################