set.seed(500)


####################################################################################
library(geosphere)
library(fields)
library(gstat)
library(sf)
library(terra)

options(scipen=999) #make R prefer displaying numbers in standard decimal notation.

#Set working directory
setwd("path_to_directory")

#loading the data
dat      <- read.csv("Practice_data_outcome.csv", header=TRUE)
head(dat)

coords <- data.frame(LONGNUM = dat$LONGNUM, LATNUM = dat$LATNUM)
coords <- na.omit(coords) #Remove rows with missing coordinates

#Euclidean distances in decimal degrees. To convert to km, multiply by 111. 
#1 degree longitude = 111.32 km at the equator
#However, as you move towards the poles, this distance decreases.
hh <- dist(coords) #matrix of distances between pairs of locations in the data
dim(hh) 
dist1 <- as.numeric(dist(coords)*111) #convert to a vector
summary(dist1)


##Great circle or spherical distance or geodetic distance in meters, 
#geosphere package. Divide by 1000 to convert to km.
#For Nigeria, this gives results similar to the Euclidean distances.
dd <- distm(coords, fun=distHaversine)
dist2 <- as.numeric(dd[lower.tri(dd, diag = FALSE)])/1000
summary(dist2)


#Convert the distances from lon-lat to UTM coordinates

# Convert the data (coords) to an sf object
SP_longlat <- st_as_sf(coords, coords = c("LONGNUM", "LATNUM"), crs = 4326)
#EPSG:4326 identifies the WGS 84 geographic CRS, which uses latitude and 
#longitude coordinates in degrees
#EPSG is a defunct organization that originally developed a registry of codes 
#identifying coordinate reference systems (CRSs)


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

