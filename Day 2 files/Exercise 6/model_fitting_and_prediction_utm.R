####Spatial analysis using UTM coordinates

#Load libraries
library(INLA)
library(inlabru)
library(terra); library(maptools)
library(gtools); library(sp); library(spdep)
library(sf)
library(ggplot2)

set.seed(500)

#Set working directory (folder for the current exercise)
setwd("path_to_directory")

#loading the data
vaxdata <- read.csv("Data/Practice_data_outcome.csv", header=TRUE)
vaxcov  <- read.csv("Data/Covariates_selected.csv", header=TRUE)

#Align both data sets
data.merge <- merge(vaxdata, vaxcov[,1:11], by="ID")

#Delete clusters where <=1 individual was sampled
zero.clust <- which(is.na(data.merge$total_sampled)|data.merge$total_sampled<=1)
if (length(zero.clust)>0){
  data.merge <- data.merge[-zero.clust,]
}


#Coordinates - transform to UTM
pts <- st_as_sf(data.merge, coords = c("LONGNUM", "LATNUM"), crs = 4326)  # WGS84 geographic
utm_zone <- floor((mean(data.merge$LONGNUM) + 180) / 6) + 1 #Nigeria is mostly in UTM zone 32N
epsg <- 32600 + utm_zone        # Northern Hemisphere; use 32700 + zone for Southern 
pts_utm <- st_transform(pts, crs = epsg)

coords_utm <- st_coordinates(pts_utm)
head(coords_utm)

data.merge$UTM_Easting <- coords_utm[, 1]
data.merge$UTM_Northing <- coords_utm[, 2]


#Relabelling variables - optional
data.merge <- data.merge %>%
  mutate(
    Numvacc = vax_count,
    weights = total_sampled,
    
    xp1 = l_Pigs_density,
    xp2 = Dist_to_edge_cult_areas,
    xp3 = l_Travel_time_urban_areas,
    xp4 = l_Distance_to_conflicts,
    xp5 = factor(Urban_rural1, levels = c(0, 1),
                 labels = c("rural", "urban")),
    xp6 = Avg_precipitation / 1000,  # scaling
    xp7 = Avg_EVI / 1000,            # scaling
    xp8 = l_Land_surface_temp,
    xp9 = Proximity_national_borders
  )

#Create a dummy variable for the urban category  
#in xp5
data.merge <- data.merge %>%
  mutate(
    xp5_urban = as.numeric(xp5 == "urban")
  )


#Convert data.merge to an sf object
data.merge <- st_as_sf(data.merge, coords = c("UTM_Easting", "UTM_Northing"), 
                       crs = epsg)


#meshfit: fine triangulated mesh
#Read in admin0 shapefile for Nigeria 
shp_ng  <- st_read("Shapefiles/gadm41_NGA_0.shp", crs = 4326) # WGS84 geographic

shp_ng <- st_transform(shp_ng, crs = epsg) #convert to UTM

c.bnd <- st_coordinates(shp_ng)[,1:2]

meshfit <- inla.mesh.2d(loc=coords_utm, loc.domain=c.bnd, 
                        max.edge=c(22264 + 6000, 66792 + 1000), 
                        cutoff=11132 + 6000)	#mesh		

#Plot the mesh and add coordinates and boundary points
plot(meshfit); plot(st_geometry(shp_ng), add=TRUE)
points(coords_utm, pch=21, bg=1, col="white", cex=1) 


#For priors
alpha <- 2 #This implies that the smoothness parameter nu=1
r0    <- 53433.6 #This is 5% of the extent of Nigeria in the north-south direction 
#(i.e. 0.05*(ymax-ymin) of any of the rasters above)

#Matern SPDE model object using inla.pcmatern
spde <- inla.spde2.pcmatern(mesh=meshfit, alpha=alpha, prior.range=c(r0, 0.01), 
                            prior.sigma=c(3, 0.01)) 


# Define additional prior distributions 
hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01))) #IID variance
control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, 
                     prec.intercept=1/1000) #Regression coefficients


#Model fitting using inlabru
cmp <- ~ -1 + Intercept(1) + xp1 + xp2 + xp3 + xp4 + xp5_urban +
  xp6 + xp7 + xp8 + xp9 +
  field(geometry, model = spde) + #spatial random effect
  f.iid(ID, model="iid", hyper = hyper.prec, constr = TRUE) #iid random effect

likelihood <- bru_obs(
  Numvacc ~ .,
  family = "binomial",
  Ntrials = data.merge$weights,
  data = data.merge)

mod.fit <- bru(
  components = cmp,
  likelihood,
  options = list(
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE, config = TRUE),
    control.fixed = control.fixed,
    control.inla = list(
      int.strategy = "eb")) #int.strategy = "eb", strategy = 'simplified.laplace'
)


summary(mod.fit)

#To view the estimates of xp5
#mod.fit$summary.random$xp5


###Prediction at 10 x 10 km resolution

#Read in prediction covariate raster data at 10 km resolution
#Prefix "l_" means that the variable should be log-transformed
l_Pigs_density  	        <- rast("Raster_files/Pigs_density.tif")
Dist_to_edge_cult_areas  	<- rast("Raster_files/Dist_to_cult_areas.tif") 
l_Travel_time_urban_areas <- rast("Raster_files/Travel_time_to_urban_areas.tif")
l_Distance_to_conflicts  	<- rast("Raster_files/Dist_to_conflicts.tif")
Avg_precipitation	        <- rast("Raster_files/Avg_precipitation.tif")
Avg_EVI	                  <- rast("Raster_files/Avg_modis_EVI.tif")
l_Land_surface_temp	      <- rast("Raster_files/Avg_modis_daytime_land_surfaceTemp.tif")
Proximity_national_borders <- rast("Raster_files/Proximity_nat_borders.tif")

#Read in urban-rural covariate at 10 km resolution
Urban_rural <- rast("Raster_files/NGA_urban_rural_10km.tif")

#Relabel the covariates and apply log transformation
xp1 <- l_Pigs_density; xp1 <- log(xp1 + 0.05)
xp2 <- Dist_to_edge_cult_areas
xp3 <- l_Travel_time_urban_areas; xp3 <- log(xp3 + 0.05)
xp4 <- l_Distance_to_conflicts; xp4 <- log(xp4 + 0.05)
xp5 <- Urban_rural
xp6 <- Avg_precipitation; xp6 <- xp6/1000 # Note scaling 
xp7 <- Avg_EVI; xp7 <- xp7/1000 #Note scaling
xp8 <- l_Land_surface_temp; xp8 <- log(xp8 + 0.05)
xp9 <- Proximity_national_borders

#Population data for population-weighted aggregation
#pop <- rast("Raster_files/Worldpop_under5s_2018.tif")

#Remove original raster files to save memory
rm(l_Pigs_density); rm(Dist_to_edge_cult_areas); rm(l_Travel_time_urban_areas)
rm(l_Distance_to_conflicts); rm(Avg_precipitation); rm(Avg_EVI)
rm(l_Land_surface_temp); rm(Proximity_national_borders); rm(Urban_rural)


#Combine prediction grid and covariates 
pred.dat <- as.data.frame(c(xp1, xp2, xp3, xp4, xp5, xp6, xp7, xp8, xp9), 
                          xy = TRUE, na.rm=FALSE)
head(pred.dat)

colnames(pred.dat) <- c("x", "y", "xp1", "xp2", "xp3", "xp4", "xp5", "xp6", 
                        "xp7", "xp8", "xp9")

#Convert prediction coordinates to an sf object
pts.grid_sf <- st_as_sf(pred.dat[ ,1:2], coords = c("x", "y"),
                           crs = crs(xp1)) # WGS84 geographic

#Transform prediction coordinates to UTM
pts.grid_utm <- st_transform(pts.grid_sf, crs = epsg) #Convert to UTM

Pred_coords <- st_coordinates(pts.grid_utm)
head(Pred_coords)

colnames(Pred_coords) <- c("UTM_Easting", "UTM_Northing")

#Add converted coordinates to data frame
pred.dat <- cbind(Pred_coords, pred.dat[,-c(1:2)])


#Identify all grid cells with missing covariate values and remove these
ind     <- apply(pred.dat, 1, function(x) any(is.na(x)))
miss    <- which(ind==TRUE)
nonmiss <- which(ind==FALSE)

#Complete cases
pred.dat.comp <- pred.dat[nonmiss,]

#Add new id variable for the iid RE to the prediction data frame
pred.dat.comp$ID <- 1:nrow(pred.dat.comp)

#Declare xp5 as a factor variable within the prediction data frame
pred.dat.comp$xp5 <- factor(pred.dat.comp$xp5, levels = c(0, 1), 
                         labels = c("rural", "urban"))

pred.dat.comp <- pred.dat.comp %>%
  mutate(
    xp5_urban = as.numeric(xp5 == "urban"))


#Convert prediction data frame to an sf object
pred.dat.comp <- st_as_sf(pred.dat.comp, coords = c("UTM_Easting", "UTM_Northing"), 
                          crs = epsg)

#Prediction samples
n.samples <- 1000 #Can change to 100 for trial runs

#Generate samples from iid term
iid.sd <- sqrt(1/mod.fit$summary.hyperpar["Precision for f.iid", 1]) 
iid.samp <- matrix(rnorm(nrow(pred.dat.comp)*n.samples, 0, iid.sd), 
                   nrow=nrow(pred.dat.comp), ncol=n.samples)

#Generate samples for other model components
s.samp <- generate(mod.fit, pred.dat.comp,
                  ~  Intercept + xp1 + xp2 + xp3 + xp4 + xp5_urban +
                    xp6 + xp7 + xp8 + xp9 + field, 
                  n.samples = n.samples) 

#Combine both samples and take the inverse logit
samp <- inv.logit(s.samp + iid.samp)

pred.obs <- data.frame(t(apply(samp, 1, FUN=function(x){ 
  c(mean(x), sd(x), quantile(x, probs=c(0.025,0.5,0.975), 
                             na.rm = TRUE))}))) 

colnames(pred.obs) <- c("mean", "sd", "low", "median", "up")

#Sample raster
#Note that the sample raster is in geographic coordinates
#You can project this to UTM if you wish
samp.rast  <- rast("Raster_files/Avg_modis_EVI.tif")

#mean
ll <- rep(NA, length(ind)); ll[nonmiss] <- pred.obs$mean
rr.mean <- rast(samp.rast); 
values(rr.mean) <- ll

#median
ll <- rep(NA, length(ind)); ll[nonmiss] <- pred.obs$median
rr.median <- rast(samp.rast); 
values(rr.median) <- ll

#sd
ll <- rep(NA, length(ind)); ll[nonmiss] <- pred.obs$sd
rr.sd <- rast(samp.rast); 
values(rr.sd) <- ll

#low
ll <- rep(NA, length(ind)); ll[nonmiss] <- pred.obs$low
rr.low <- rast(samp.rast); 
values(rr.low) <- ll

#up
ll <- rep(NA, length(ind)); ll[nonmiss] <- pred.obs$up
rr.up <- rast(samp.rast); 
values(rr.up) <- ll

#95% CI width
rr.95CIwidth <- rr.up - rr.low


# save output raster layers
writeRaster(rr.mean, paste0("pred_mean_utm_data.tif"), overwrite=TRUE)
writeRaster(rr.sd,   paste0("pred_sd_utm_data.tif"), overwrite=TRUE)
writeRaster(rr.low,  paste0("pred_low_utm_data.tif"), overwrite=TRUE)
writeRaster(rr.up,   paste0("pred_up_utm_data.tif"), overwrite=TRUE)
writeRaster(rr.med,  paste0("pred_median_utm_data.tif"), overwrite=TRUE)
writeRaster(rr.95CIwidth,  paste0("pred_95_perc_CI_width_utm_data.tif"), overwrite=TRUE)



