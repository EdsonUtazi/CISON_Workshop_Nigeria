#This R script performs K-fold cross-validation and computes some validation statistics which can be averaged 
#after the script has been run

#Load libraries
library(INLA)
library(inlabru)
library(terra); library(maptools)
library(gtools); library(sp); library(spdep)
library(sf)
library(ggplot2)

set.seed(500)

#Set working directory - folder for the current exercise
setwd("path_to_directory")

#loading the data
vaxdata <- read.csv("Data/Practice_data_outcome.csv", header=TRUE)
vaxcov  <- read.csv("Data/Covariates_selected.csv", header=TRUE)

#Align both data sets
data.merge.all <- merge(vaxdata, vaxcov[,1:11], by="ID")

#Delete clusters where <=1 individual was sampled
zero.clust <- which(is.na(data.merge.all$total_sampled)|data.merge.all$total_sampled<=1)
if (length(zero.clust)>0){
  data.merge.all <- data.merge.all[-zero.clust,]
}

#Coordinates
coords.all    <- cbind(data.merge.all$LONGNUM, data.merge.all$LATNUM)

#ID
ID.all <- data.merge$ID


#######################################################
#Start k-fold cross-validation loop
cv  <- 10   #i.e. 10-fold cross validation
lim <- floor(nrow(coords.all)/cv) #No of observations in each fold

#Random sample for random allocation of validation clusters
srand <- sample(1:nrow(coords.all),nrow(coords.all), replace=FALSE)

#Output matrix
val.out <- matrix(0, cv, 5) #4 model-evaluation criteria

#Cross-validation loop
for (kk in 1:cv){
 if (kk < cv) {qq <- (((kk-1)*lim)+1):(lim*kk); samp.c <- srand[qq]}
 if (kk == cv) {qq <- (((kk-1)*lim)+1):nrow(coords.all); samp.c <- srand[qq]}
 
coords.nc 	<- coords.all[samp.c,]
Numvacc.nc	<- data.merge.all$vax_count[samp.c]
weights.nc	<- data.merge.all$total_sampled[samp.c]
vaxcov.nc	<- data.merge.all[samp.c, 7:15]
ID.nc <- ID.all[samp.c]

yp.nc=np.nc=rep(NA, length(Numvacc.nc))

#Use the rest for model estimation
coords 	<- coords.all[-samp.c,]
Numvacc	<- data.merge.all$vax_count[-samp.c]
weights	<- data.merge.all$total_sampled[-samp.c]
vaxcov	<- data.merge.all[-samp.c, 7:15]
ID <- ID.all[-samp.c]

#Data frame for model-fitting
data.mod <- data.frame(
    ID,
    Numvacc,
    weights,
    xp1 = vaxcov$l_Pigs_density,
    xp2 = vaxcov$Dist_to_edge_cult_areas,
    xp3 = vaxcov$l_Travel_time_urban_areas,
    xp4 = vaxcov$l_Distance_to_conflicts,
    xp5 = factor(vaxcov$Urban_rural1, levels = c(0, 1),
                 labels = c("rural", "urban")),
    xp6 = vaxcov$Avg_precipitation / 1000,  # scaling
    xp7 = vaxcov$Avg_EVI / 1000,            # scaling
    xp8 = vaxcov$l_Land_surface_temp,
    xp9 = vaxcov$Proximity_national_borders,
    LONGNUM = coords[,1],
    LATNUM = coords[,2])

#Create a dummy variable for the urban category  
#in xp5
data.mod <- data.mod %>%
  mutate(xp5_urban = as.numeric(xp5 == "urban"))

#Convert data.mod to an sf object
data.mod <- st_as_sf(data.mod, coords = c("LONGNUM", "LATNUM"), 
                       crs = epsg)

#Data frame for validation
data.val <- data.frame(
  ID = ID.nc,
xp1     = vaxcov.nc$l_Pigs_density,
xp2     = vaxcov.nc$Dist_to_edge_cult_areas,
xp3	    = vaxcov.nc$l_Travel_time_urban_areas,
xp4	    = vaxcov.nc$l_Distance_to_conflicts,
xp5	    = vaxcov.nc$Urban_rural1,
xp6     = vaxcov.nc$Avg_precipitation/1000, #Note the scaling
xp7     = vaxcov.nc$Avg_EVI/1000, #Note the scaling
xp8     = vaxcov.nc$l_Land_surface_temp,
xp9     = vaxcov.nc$Proximity_national_borders,
LONGNUM = coords.nc[,1],
LATNUM = coords.nc[,2])

data.val <- data.val %>%
  mutate(xp5_urban = as.numeric(xp5 == "urban"))

data.val <- st_as_sf(data.val, coords = c("LONGNUM", "LATNUM"), 
                     crs = epsg)

#meshfit: fine triangulated mesh
#Read in admin 0 shapefile for Nigeria for creating a mesh
shp_ng  <- st_read("Shapefiles/gadm41_NGA_0.shp", crs = 4326) # WGS84 geographic
c.bnd <- st_coordinates(shp_ng)[,1:2]

meshfit <- inla.mesh.2d(loc=st_coordinates(data.merge[,1:2]), 
                        loc.domain=c.bnd, 
                        max.edge=c(0.3, 0.6), cutoff=0.15)	#mesh	

#For priors
alpha <- 2 #This implies that the smoothness parameter nu=1
r0    <- 0.48 #This is 5% of the extent of Nigeria in the north-south direction 
              #(i.e. 0.05*(ymax-ymin) of any of the rasters above)

#Matern SPDE model object using inla.pcmatern
spde <- inla.spde2.pcmatern(mesh=meshfit, alpha=alpha, prior.range=c(r0, 0.01), 
                            prior.sigma=c(3, 0.01)) 

# Define additional prior distributions 
hyper.prec = list(theta = list(prior="pc.prec", param=c(3,0.01))) #IID variance
control.fixed = list(mean=0, prec=1/1000, mean.intercept=0, prec.intercept=1/1000) #Regression coefficients

#Model fitting using inlabru
cmp <- ~ -1 + Intercept(1) + xp1 + xp2 + xp3 + xp4 + xp5_urban +
  xp6 + xp7 + xp8 + xp9 +
  field(geometry, model = spde) + #spatial random effect
  f.iid(ID, model="iid", hyper = hyper.prec, constr = TRUE) #iid random effect

likelihood <- bru_obs(
  Numvacc ~ .,
  family = "binomial",
  Ntrials = data.mod$weights,
  data = data.mod)

mod.fit <- bru(
  components = cmp,
  likelihood,
  options = list(
    control.compute = list(waic = TRUE, dic = TRUE, cpo = TRUE, config = TRUE),
    control.fixed = control.fixed,
    control.inla = list(
      int.strategy = "eb")) #int.strategy = "eb", strategy = 'simplified.laplace'
)


#Prediction samples
n.samples <- 1000 #Can change to 100 for trial runs

#Generate samples from iid term
iid.sd <- sqrt(1/mod.fit$summary.hyperpar["Precision for f.iid", 1]) 
iid.samp <- matrix(rnorm(nrow(data.val)*n.samples, 0, iid.sd), 
                   nrow=nrow(data.val), ncol=n.samples)

#Generate samples for other model components
s.samp <- generate(mod.fit, data.val,
                   ~  Intercept + xp1 + xp2 + xp3 + xp4 + xp5_urban +
                     xp6 + xp7 + xp8 + xp9 + field, 
                   n.samples = n.samples) 

#Combine both samples and take the inverse logit
inv.linpred <- inv.logit(s.samp + iid.samp)

pred.val    <- data.frame(t(apply(inv.linpred, 1, FUN=function(x){ 
  c(mean(x), sd(x), quantile(x, probs=c(0.025,0.5,0.975)))}))) 
colnames(pred.val) <- c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
fitted.mean.val <- pred.val$mean
fitted.low.val  <- pred.val$`0.025quant`
fitted.up.val   <- pred.val$`0.975quant`

prob.val <- Numvacc.nc/weights.nc #Observed probabilities
corr     <- cor(fitted.mean.val,prob.val) #Correlation
rsq.val  <- (cor(fitted.mean.val, prob.val))^2 #R-squred
RMSE.val <- sqrt(sum((fitted.mean.val-prob.val)^2)/length(prob.val)) #RMSE
MAE      <- sum(abs(fitted.mean.val-prob.val))/length(prob.val) #MAE

#count <- 0
#for(r in 1:length(prob.val)){
#  if ((prob.val[r] >= fitted.low.val[r]) && (prob.val[r] <= fitted.up.val[r])) count <- count + 1
#}
#cov.rate.val <- (count/length(prob.val))*100

perc_bias <- (sum(fitted.mean.val-prob.val)/sum(prob.val))*100

val.out[kk, ] <- c(corr, rsq.val, RMSE.val, MAE, perc_bias)
print(val.out[kk, ])
}

colnames(val.out) <- c("correl", "rsq.val", "RMSE.val", "MAE", "perc_bias")
write.csv(val.out, paste0("Cross_validation_results_random.csv"))
