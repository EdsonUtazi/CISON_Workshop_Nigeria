

######################################################################################################
#Variograms

library(gstat)
library(INLA)
#library(PrevMap)

#Set working directory
setwd("path_to_folder")

#loading the data
dat      <- read.csv("Practice_data_outcome.csv", header=TRUE)
dat.cov  <- read.csv("Covariates_selected.csv", header=TRUE)

data.merge <- merge(dat, dat.cov[,1:11], by="ID")
coords <- data.frame(LONGNUM = data.merge$LONGNUM, LATNUM = data.merge$LATNUM)


#Fit variogram using gstat package

#method 1: Fit a non-spatial binomial regression model to extract the residuals
#Covariate names
covnames <- names(dat.cov)[2:10]

#Formula 
form <- paste("cbind(vax_count, total_sampled-vax_count)","~", paste(covnames, collapse=" + "))

#Fit the model
fit  <- glm(form, data = data.merge, family = binomial(logit))

#Extract the residuals
data.merge$residuals1 <- resid(fit) #Deviance residuals


#Method 2: Fit a Binomial GLM in INLA with an iid random effect, and extract the
#posterior estimaets of the iid random effect and use these to fit the variogram

# Create an observation-level index
data.merge$iid <- 1:nrow(data.merge)
data.merge$intercept <- rep(1, nrow(data.merge))

# Formula
#Note that since Urban_rural1 is binary (0/1), there is no need to create a separate
#dummy variable for this
form2 <- paste(
  "vax_count ~ -1 + intercept +",
  paste(covnames, collapse = " + "),
  "+ f(iid, model = 'iid')"
)

# Fit the model
fit2 <- inla(
  as.formula(form2),
  data = data.merge,
  Ntrials = data.merge$total_sampled,
  family = "binomial",
  control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE)
)

# Extract posterior mean of the iid random effect
data.merge$residuals2 <- fit2$summary.random$iid$mean



###########################################################################
#Calculate empirical variogram using the first residuals from method 1
vgm1 <- variogram(residuals1 ~ 1, data = data.merge, ~LONGNUM+LATNUM)  
plot(vgm1) #Empirical variogram

#Fit a variogram model - Mat, Exp, Sph, etc
#Initial valeus for the parameters of the model can be provided tas well. These parameters
#are estimated using the specified model if not provided
#foo <- fit.variogram(vgm1, vgm("Exp"))
foo <- fit.variogram(vgm1, vgm("Mat", fit.kappa = TRUE))
#foo <- fit.variogram(vgm1, vgm("Sph"))
#foo1 <- fit.variogram(vgm1, vgm(0.4768, "Exp", 0.5627, nugget=0.5))

#View estimated parameters
foo

#Plot fitted variogram model
plot(vgm1, model=foo) #variogram + fitted line

#Plot empirical variogram, variogram fit as well as estimated parameters
plot(vgm1$dist, vgm1$gamma, ylim=c(0,2), col = "cornflowerblue", xlab = "distance", ylab = "semivariance")
lines(variogramLine(vgm(psill = foo$psill[2], "Mat", range = foo$range[2], 
                        kappa = foo$kappa[2], nugget = foo$psill[1]), 10), type = 'l') #psill, model, range, nugget
abline(h = foo$psill[1], col = "blue")
abline(v = foo$range[2], col = "red")
abline(h = foo$psill[1] + foo$psill[2], col = "green")


#Allocate the residuals at random to the locations and plot the empirical variogram
data.merge$residuals.rand <- sample(data.merge$residuals1, nrow(data.merge), replace = FALSE)

#Variogram - observe that the variogram is no longer an increasing function of distance
vgm1 <- variogram(residuals.rand~1, data = data.merge, ~LONGNUM+LATNUM)  #Empirical variogram
plot(vgm1) #Empirical variogram plot
####################################################################################################################



###########################################################################
#Calculate empirical variogram using the first residuals from method 2
vgm1 <- variogram(residuals2 ~ 1, data = data.merge, ~LONGNUM+LATNUM)  
plot(vgm1) #Empirical variogram

#Fit a variogram model - Mat, Exp, Sph, etc
#Initial valeus for the parameters of the model can be provided tas well. These parameters
#are estimated using the specified model if not provided
#foo <- fit.variogram(vgm1, vgm("Exp"))
foo <- fit.variogram(vgm1, vgm("Mat", fit.kappa = TRUE))
#foo <- fit.variogram(vgm1, vgm("Sph"))
#foo1 <- fit.variogram(vgm1, vgm(0.4768, "Exp", 0.5627, nugget=0.5))

#View estimated parameters
foo

#Plot fitted variogram model
plot(vgm1, model=foo) #variogram + fitted line

#Plot empirical variogram, variogram fit as well as estimated parameters
plot(vgm1$dist, vgm1$gamma, ylim=c(0,0.8), col = "cornflowerblue", xlab = "distance", ylab = "semivariance")
lines(variogramLine(vgm(psill = foo$psill[2], "Mat", range = foo$range[2], 
                        kappa = foo$kappa[2], nugget = foo$psill[1]), 10), type = 'l') #psill, model, range, nugget
abline(h = foo$psill[1], col = "blue")
abline(v = foo$range[2], col = "red")
abline(h = foo$psill[1] + foo$psill[2], col = "green")


#Allocate the residuals at random to the locations and plot the empirical variogram
data.merge$residuals.rand <- sample(data.merge$residuals2, nrow(data.merge), replace = FALSE)

#Variogram - observe that the variogram is no longer an increasing function of distance
vgm1 <- variogram(residuals.rand~1, data = data.merge, ~LONGNUM+LATNUM)  #Empirical variogram
plot(vgm1) #Empirical variogram plot
####################################################################################################################











###############################################################################################
# #Variogram using the PrevMap package - includes a permutation test which basically involves assigning
# #the residuals at random to the locations multiple times, and then using the resulting variogram estimates to generate 
# #a variogram envelope based on the quantiles (usually 2.5% and 97.5%) of the estimates
# 
# #Compare this with PrevMap variog.
# 
# # Testing for residual spatial correlation
# # using the variogram
# spat.corr.diagnostic(vax_count ~ l_Pigs_density + Dist_to_edge_cult_areas + l_Travel_time_urban_areas +
#                        l_Distance_to_conflicts + Avg_precipitation + Proximity_national_borders + 
#                        Avg_EVI + l_Land_surface_temp + Urban_rural1,
#                      units.m = ~ total_sampled,
#                      coords = ~LONGNUM+LATNUM,
#                      likelihood = "Binomial", n.sim = 200, kappa=1,
#                      data = data.merge)
# 
# #Delete clusters where TotChild is zero from vaxdata only
# zero.clust <- which(is.na(data.merge$total_2_1)|data.merge$total_2_1<=0)
# if (length(zero.clust)>0){
#   data.merge.2 <- data.merge[-zero.clust,]
# }

###############################################################################################

