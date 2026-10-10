
library(INLA)
library(INLAspacetime) #for fmesher
library(sf)
library(abind)
library(splancs)
library(patchwork)

#Data relates to PM10 concentration in the Piemonte region of Northern Italy
#measured between October 2005 and March 2006
#The daily data were measured at 24 monitoring stations


setwd("path_to_folder")


Piemonte_data <- read.csv("Piemonte_data_byday.csv")
head(Piemonte_data)

coordinates <- read.csv("coordinates.csv")

borders <- read.csv("Piemonte_borders.csv")

Piemonte_shp <- st_read("ASL_Piemonte.shp")
plot(st_geometry(Piemonte_shp),
  col = NA, border = "black")

#Number of stations
n_stations <- length(unique(Piemonte_data$Station.ID))

#Number of days
n_days <- nrow(Piemonte_data)/n_stations 

#Create a new time index
Piemonte_data$time <- rep(1:n_days, each = n_stations)
#View(Piemonte_data)

head(coordinates)

#Extract unique coordinates from the data
nrow(coordinates)



#Log transform the outcome variable to encourage normality
Piemonte_data$logPM10 <- log(Piemonte_data$PM10)

#Standardize the covariates prior to model fitting
mean_covs <- colMeans(Piemonte_data[ ,c(3,6:10)])
sd_covs   <- apply(Piemonte_data[ ,c(3,6:10)], 2, sd)
Piemonte_data[ ,c(3,6:10)] <- scale(Piemonte_data[ ,c(3,6:10)], 
                              center = mean_covs,
                              scale = sd_covs)


#Create mesh
Piemonte_mesh <- inla.mesh.2d(loc=cbind(coordinates$UTMX, coordinates$UTMY),
                              loc.domain=borders,
                              max.edge=c(50,1000), #coarse mesh to speed up computation
                              offset=c(10, 140))

plot(Piemonte_mesh)
lines(borders, col = "blue")
points(cbind(coordinates$UTMX, coordinates$UTMY), col="red")

#Set prior distributions
st_bbox(Piemonte_shp) #Extent of Piemonte

r0 <- 13 #0.05*(5145805.5-4879335.0)/1000 #This is 5% of the extent of Piemonte in the 
#north-south direction (i.e. 0.05*(ymax-ymin))
#Note the division by 1000 as coordinates are already given in km

rprior <- list(theta = list(prior = "pccor1", 
                            param = c(0.5, 0.9))) #autocorrelation parameter
control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, 
                     prec.intercept=1/1000) #Fixed effects
hyper.prec.surv = list(theta = list(prior="pc.prec", 
                                    param=c(5,0.01))) #iid variance

#Create SPDE object
Piemonte_spde <- inla.spde2.pcmatern(mesh=Piemonte_mesh, alpha=2,
                                   prior.range = c(r0, 0.01),
                                   prior.sigma = c(3, 0.01))


#Convert data to an sf object
Piemonte_data_sf <- st_as_sf(Piemonte_data, coords = c("UTMX", "UTMY"), 
                             crs = st_crs(Piemonte_shp))


time_mapper <- bru_mapper(fm_mesh_1d(sort(unique(Piemonte_data$time))), 
                          indexed = TRUE)


#Components
cmp <- ~ -1 + Intercept(1) + A + WS + TEMP + HMIX + PREC + EMI +
  field(geometry, model = Piemonte_spde, 
    group=time, group_mapper=time_mapper, 
    control.group=list(model="ar1", hyper = rprior))

#Likelihood
likelihood <- bru_obs(
  logPM10~.,
  family = "gaussian",
  data = Piemonte_data_sf
)

mod.fit <- bru(
  components = cmp,
  likelihood,
  options = list(
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE),
    control.fixed = control.fixed,
    control.inla = list(
      int.strategy = "eb")) #int.strategy = "eb", strategy = 'simplified.laplace'
)

summary(mod.fit)



#Load prediction covariates
#These functions are from Blangiardo and Cameletti (2015)
# Load the covariate arrays (each array except for A is 56x72x182)
load(paste("Covariates/Altitude_GRID.Rdata",sep="")) #A; AltitudeGRID
load(paste("Covariates/WindSpeed_GRID.Rdata",sep="")) #WS; WindSpeedGRID
load(paste("Covariates/HMix_GRID.Rdata",sep="")) #HMIX; HMixMaxGRID
load(paste("Covariates/Emi_GRID.Rdata",sep="")) #EMI; EmiGRID
load(paste("Covariates/Temp_GRID.Rdata",sep="")) #TEMP; Mean_Temp
load(paste("Covariates/Prec_GRID.Rdata",sep="")) #PREC; Prec


# Load the Piemonte grid c(309,529),c(4875,5159),dims=c(56,72)
load(paste("Covariates/Piemonte_grid.Rdata",sep=""))

# Extract the standardized covariates for day i_day (you get a 56X72X8 matrix)
i_day <- 122
which_date <- unique(Piemonte_data$Date)[i_day]
print(paste("**---- You will get a prediction for ", which_date, "---**"))


# IMPORTANT NOTE: In practice, this standardization should be done the other 
#way round. The prediction covariate data/layers should be standardized first, 
#and then the standardization parameters applied to the input data.

# Standardise the covariates for the selected day
WindSpeedGRID_i	= (WindSpeedGRID[,,i_day] - mean_covs[2]) /sd_covs[2]
HMixMaxGRID_i	= (HMixMaxGRID[,,i_day] - mean_covs[4]) / sd_covs[4]	
EmiGRID_i		= (EmiGRID[,,i_day] - mean_covs[6]) / sd_covs[6]
Mean_Temp_i		= (Mean_Temp[,,i_day] - mean_covs[3]) / sd_covs[3]
Prec_i			= (Prec[,,i_day] - mean_covs[5]) / sd_covs[5]
Prec_i[is.na(Prec_i)]=0 #when all prec for day i are 0 the scale proc generates NA

#UTMX_std = (unique(Piemonte_grid[,1]) - mean_covs[2]) / sd_covs[2]
#UTMX_GRID = matrix(rep(UTMX_std,72),56,72) #constant for each row

#UTMY_std = (unique(Piemonte_grid[,2]) - mean_covs[3]) / sd_covs[3]
#UTMY_GRID = t(matrix(rep(UTMY_std,56),72,56)) #constant for each column

AltitudeGRID = (AltitudeGRID - mean_covs[1]) / sd_covs[1]

#--- Create the array (56x72x8) of standardized covariates for day i
covariate_array_std <- abind(AltitudeGRID,
                            #UTMX_GRID,
                            #UTMY_GRID,	
                            WindSpeedGRID_i,
                            Mean_Temp_i,
                            HMixMaxGRID_i,
                            Prec_i,
                            EmiGRID_i,
                            along=3)

dim(covariate_array_std) #[1] 56 72  6


# Set to NA the (standardized) altitude values >7 (1000 n)
elevation <- covariate_array_std[,,1]
index_mountains <- which(elevation > 7)
elevation[elevation > 7] <- NA
covariate_array_std[,,1] <- elevation

# Reshape the 3D array (56x72x6) into a dataframe (4032x6)
covariate_matrix_std <- data.frame(apply(covariate_array_std,3,
                                         function(X) c(t(X))))

colnames(covariate_matrix_std) <- colnames(Piemonte_data[,c(3,6:10)])
head(covariate_matrix_std)
#head(covariate_matrix_std1)

#mean(Prec_i)
#mean(covariate_matrix_std$PREC)

#Add time and prediction coordinates and convert to an sf object
covariate_matrix_std$time <- rep(i_day, nrow(covariate_matrix_std))
covariate_matrix_std <- cbind(covariate_matrix_std, Piemonte_grid)

covariate_matrix_std_sf <- st_as_sf(covariate_matrix_std, 
                                    coords = c("UTMX_km", "UTMY_km"), 
                                    crs = st_crs(Piemonte_shp))

#Prediction samples
n.samples <- 500 #Can change to 100 for trial runs

post.samp <- generate(mod.fit, covariate_matrix_std_sf,
                  ~  Intercept + A + WS + TEMP + HMIX + PREC + EMI +    
                    field, 
                  n.samples = n.samples) 

#Take exponent to back transform to the original scale
#This can also be done inside the generate() function
post.samp <- exp(post.samp) 

pred.out <- data.frame(t(apply(post.samp, 1, FUN=function(x){ c(mean(x), sd(x), 
                  quantile(x, probs=c(0.025,0.5,0.975), na.rm = TRUE))}))) 

colnames(pred.out) <- c("mean", "sd", "low", "median", "up")
dim(pred.out)
head(pred.out)

#Calculate the exceedance probabilities - prob. of exceeding 50 units
ff1 <- function(x) length(which(x>=50))/n.samples
pred.exceed <- apply(post.samp, 1, ff1) 
length(pred.exceed)


###For plotting
#Convert shapefile coords to km
Piemonte_shp_km <- Piemonte_shp
st_geometry(Piemonte_shp_km) <- st_geometry(Piemonte_shp_km)/1000

# Identify points inside Piemonte
Piemonte_grid1 <- Piemonte_grid; colnames(Piemonte_grid1) <- c("x", "y")
borders1 <- borders; colnames(borders1) <- c("x", "y")

inside_Piemonte <- matrix(inout(Piemonte_grid1, borders1),
                          56, 72, byrow = TRUE)
inside_Piemonte[inside_Piemonte == 0] <- NA

# Create the coordinate sequences
seq.x.grid <- seq(min(Piemonte_grid[, 1]), max(Piemonte_grid[, 1]),
                  length.out = 56)

seq.y.grid <- seq(min(Piemonte_grid[, 2]), max(Piemonte_grid[, 2]),
                  length.out = 72)



#Mean
pred.mean      <- pred.out[, "mean"]
pred.mean.grid <- matrix(pred.mean, 56, 72, byrow=T)

# Select only points inside Piemonte and set NA to the outer points 
pred.mean.grid[index_mountains] <- NA

# Keep only values inside Piemonte
inside_pred.mean.grid <- inside_Piemonte * pred.mean.grid

# Convert matrix to a terra raster
r <- rast(nrows = 72, ncols = 56,
  xmin = min(seq.x.grid), xmax = max(seq.x.grid),
  ymin = min(seq.y.grid), ymax = max(seq.y.grid),
  crs = st_crs(Piemonte_shp)$wkt)

# Reorient matrix for raster format
r[] <- as.vector(inside_pred.mean.grid[,72:1])

# Convert raster to data frame for ggplot
r_df <- as.data.frame(r, xy = TRUE, na.rm = FALSE)
names(r_df)[3] <- "PM10_mean"

p1 <- ggplot() +
  geom_raster(data = r_df, aes(x = x, y = y, fill = PM10_mean)) +
  geom_sf(data = Piemonte_shp_km, fill = NA, colour = "black", linewidth = 0.4) +
  geom_point(data = coordinates, aes(x = UTMX, y = UTMY), shape = 21,
    fill = "white", colour = "black", size = 2) +
  scale_fill_viridis_c(
    name = "Predicted \nPM10",
    option = "viridis",
    na.value = "transparent", limits=c(0, 70)
  ) +
  coord_sf(xlim = c(309, 529), ylim = c(4875, 5159), expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
    axis.title = element_blank())


#Plot exceedance probabilities
pred.exceed.grid <- matrix(pred.exceed, 56, 72, byrow=T)

# Select only points inside Piemonte and set NA to the outer points 
pred.exceed.grid[index_mountains] <- NA

# Keep only values inside Piemonte
inside_pred.exceed.grid <- inside_Piemonte * pred.exceed.grid

# Convert matrix to a terra raster
r1 <- rast(nrows = 72, ncols = 56,
          xmin = min(seq.x.grid), xmax = max(seq.x.grid),
          ymin = min(seq.y.grid), ymax = max(seq.y.grid),
          crs = st_crs(Piemonte_shp)$wkt)

# Reorient matrix for raster format
r1[] <- as.vector(inside_pred.exceed.grid[,72:1])

# Convert raster to data frame for ggplot
r1_df <- as.data.frame(r1, xy = TRUE, na.rm = FALSE)
names(r1_df)[3] <- "Exc_prob"

p2 <- ggplot() +
  geom_raster(data = r1_df, aes(x = x, y = y, fill = Exc_prob)) +
  geom_sf(data = Piemonte_shp_km, fill = NA, colour = "black", linewidth = 0.4) +
  geom_point(data = coordinates, aes(x = UTMX, y = UTMY), shape = 21,
             fill = "white", colour = "black", size = 2) +
  scale_fill_viridis_c(name = "Exc. prob.", option = "magma",
    na.value = "transparent", limits=c(0, 1)) +
  coord_sf(xlim = c(309, 529), ylim = c(4875, 5159), expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_blank())


#Combine both plots

# Combine plots side-by-side
combined_plot <- p1 + p2 +
  plot_layout(ncol = 2)

# Save as PNG
ggsave(
  filename = "combined_plot.png",
  plot = combined_plot,
  width = 12,
  height = 6,
  units = "in",
  dpi = 300
)





#The code below will average PM10 over time for each monitoring station
#and produce a plot of the averaged data
PM10_mean <- Piemonte_data %>%
  group_by(Station.ID) %>%
  summarise(PM10_mean = mean(PM10, na.rm = TRUE),
    .groups = "drop")

# Add the spatial coordinates back to the averaged data
PM10_mean <- Piemonte_data %>%
  select(Station.ID, UTMX, UTMY) %>%
  distinct() %>%
  left_join(PM10_mean, by = "Station.ID")

# Convert to sf
PM10_mean_sf <- st_as_sf(PM10_mean, coords = c("UTMX", "UTMY"),
  crs = st_crs(Piemonte_shp_km))

# Plot
point.plot <- ggplot() +
  geom_sf(data = Piemonte_shp_km, fill = NA, colour = "black", linewidth = 0.4) +
  geom_sf(data = PM10_mean_sf, aes(colour = PM10_mean), size = 3) +
  scale_colour_viridis_c(name = "Mean PM10", option = "viridis") +
  coord_sf(xlim = c(309, 529), ylim = c(4875, 5159), expand = FALSE) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.title = element_blank())

point.plot
