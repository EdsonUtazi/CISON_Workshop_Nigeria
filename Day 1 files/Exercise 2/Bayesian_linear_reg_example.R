################################################################################
################Linear regression example using INLA and inlabru################
################################################################################


library(INLA)
library(inlabru)

# Setting working directory
setwd("path_to_folder")

# Read in the csv file
mydata <- read.csv("NMMAPSraw.csv")
head(mydata)
class(mydata)
summary(mydata)

# Plot histogram of pm10
hist(mydata$pm10)

#It looks like the pm10 requires a transformation to encourage normality,
#but we will go ahead and fit the model for illustration.


# Fit the simple linear regression based on Bayesian framework
formula  <-  pm10 ~ 1 + temperature 
#formula <- pm10 ~ temperature #this formula will also work
mod.inla <- inla(formula,family="gaussian",data=mydata)
round(mod.inla$summary.fixed[,1:5],3)

# Plotting the estimated parameters (posteriors)
plot(mod.inla$marginals.fixed[[1]],type="l",main="",ylab="",xlab=expression(beta[0]))
plot(mod.inla$marginals.fixed[[2]],type="l",main="",ylab="",xlab=expression(beta[1]))


# Fit the frequentist linear regression and compare with the Bayesian results
summary(lm(formula,data=mydata))# Very similar results with INLA


# We can potentially change the prior for the regression parameters 
# We can achieve this in R-INLA using the option control.fixed as below:

# Fit the simple linear regression based on Bayesian framework
mod.inla2 <- inla(formula,family="gaussian",data=mydata, 
          control.fixed=list(mean=0, prec=1,
                             mean.intercept=0, prec.intercept=0.0001))

round(mod.inla2$summary.fixed[,1:5],3)

# Plotting the estimated parameters (posteriors)
plot(mod.inla2$marginals.fixed[[1]],type="l",main="",ylab="",
     xlab=expression(beta[0]))
plot(mod.inla2$marginals.fixed[[2]],type="l",main="",ylab="",
     xlab=expression(beta[1]))


# We can also change the prior for the precision as below: 
mod.inla3 <- inla(formula,family="gaussian", data=mydata, 
                  control.fixed=list(mean=0, prec=1,
                    mean.intercept=0, prec.intercept=0.001),
                  control.family=list(hyper=list(prec=list(prior = "pc.prec",
                    param = c(1, 0.05)))))
                                                   
round(mod.inla3$summary.fixed[,1:5],3)

#To view and compare the estimated parameters for all the models
mod.inla$summary.hyperpar
mod.inla2$summary.hyperpar
mod.inla3$summary.hyperpar


#To fit the model using inlabru
#Components
cmp <- ~ Intercept(1) + temperature

#Likelihood
likelihood <- bru_obs(
  formula = pm10 ~ .,
  family  = "gaussian",
  data    = mydata)

#Model
mod.bru <- bru(cmp, likelihood,
               options = list(control.fixed = list(mean=0, prec=1,
                                      mean.intercept=0, prec.intercept=0.001),
                 control.family=list(hyper=list(prec=list(prior = "pc.prec",
                                      param = c(1, 0.05))))))

summary(mod.bru)
################################################################################










