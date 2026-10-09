
library(INLA)
library(gganimate)
library(sf)
library(gifski)


setwd("path_to_folder")

#Years: 1968 to 1988
#Number of counties: 88 

data     <- read.csv("OhioRespMort.csv", header=T)
Ohio_shp <- st_read("tl_2010_39_county00.shp")

names(Ohio_shp)
Ohio_shp$NAME00 #Names of counties

#Create a new ID/numbering of the counties in the data to match 
#the order in the shapefile
data$county.1 <- NA
nn <- unique(data$NAME)
for (i in nn){
  dd <- which(Ohio_shp$NAME00==i)
  data$county.1[data$NAME==i] <- dd
}

#Order data by county.1 and year - optional
data <- data[order(data$county.1, data$year),]
head(data)

#Check for accurary
subset(data, NAME=="Franklin") #Should be #2 in county.1


#Check crs of shapefile
crs(Ohio_shp, describe = TRUE)
Ohio_shp <- st_as_sf(Ohio_shp) #convert to an sf object

#To use the leaflet R package to create maps with the shapefile, we will convert
#it to the geographic CRS
#Ohio_shp <- st_transform(Ohio_shp, crs = 4269) #Transform to lat/lon coords


#Calculate the standardized incidence ratio (SIR)
data$SIR <- data$y / data$E

#Add the data frame to the shapefile
Ohio_shp$NAME <- Ohio_shp$NAME00
Ohio_shp <- Ohio_shp %>%
  left_join(data, by = "NAME")


#Model fitting
#Create a neighbourhood list using poly2nb and the original shapefile
Ohio_shp1 <- st_read("tl_2010_39_county00.shp")
nb.Ohio <- poly2nb(Ohio_shp1) #neighbourhood list
head(nb.Ohio) 

#Convert to a file to be used by INLA
nb2INLA("graph.Ohio.adj", nb.Ohio)
graph.inla <- inla.read.graph(filename = "graph.Ohio.adj")

#Ohio.adj <- inla.read.graph(filename = "Ohio.graph")
#Ohio.adj$nbs[[1]]
#graph.inla$nbs[[2]]

# plot(st_geometry(Ohio_shp1),
#      col = "grey90",
#      border = "grey50",
#      main = "Ohio counties")
# 
# plot(st_geometry(Ohio_shp1[2, ]),
#      col = "red",
#      border = "black",
#      add = TRUE)


#Create ID variables for the random effects to be used in the model
Ohio_shp$ID.area <- Ohio_shp$county.1
Ohio_shp$ID.area.int <- Ohio_shp$ID.area
Ohio_shp$ID.year <- Ohio_shp$year
Ohio_shp$ID.year1 <- Ohio_shp$year
Ohio_shp$ID.year.int <- Ohio_shp$year
Ohio_shp$ID.area.year <- paste(Ohio_shp$ID.area, Ohio_shp$ID.year, sep=".")


#Set the prior distributions
#prior_prec <- list(theta = list(prior="pc.prec", param=c(3, 0.01)))
#control.fixed <- list(mean=0, prec=1/1000, mean.intercept=0, 
#                      prec.intercept=1/1000) #for regression coefficients


#Dynamic trend model
formula.dt <- y ~ f(ID.area, model = "bym", graph=graph.inla) +
  f(ID.area.year, model = "iid") + ID.year1

mod.dt <- inla(formula.dt,
               family = "poisson", data = Ohio_shp, E = E,
               control.predictor = list(compute = TRUE),
               control.fixed = control.fixed,
               control.compute = list(return.marginals.predictor = TRUE, 
                                      waic = TRUE, dic = TRUE))

summary(mod.dt)


#Type I interaction
formula.T1 <- y ~ f(ID.area, model = "bym", graph=graph.inla,
                    scale.model = TRUE, constr = TRUE) +
  f(ID.year, model = "rw2", scale.model = TRUE, constr = TRUE) +
  f(ID.year1, model = "iid", constr = TRUE) +
  f(ID.area.year, model = "iid", constr = TRUE)

mod.T1 <- inla(formula.T1,
               family = "poisson", data = Ohio_shp, E = E,
               control.predictor = list(compute = TRUE),
               control.fixed = control.fixed,
               control.compute = list(return.marginals.predictor = TRUE, 
                                      waic = TRUE, dic = TRUE))

summary(mod.T1)


#Type II interaction
formula.T2 <- y ~ f(ID.area, model = "bym", graph=graph.inla, 
                    scale.model = TRUE, constr = TRUE) +
  f(ID.year, model = "rw2", scale.model = TRUE, constr = TRUE) +
  f(ID.area.int, model = "iid", group = ID.year.int,
    control.group = list(model = "rw2"), constr = TRUE)

mod.T2 <- inla(formula.T2,
               family = "poisson", data = Ohio_shp, E = E,
               control.predictor = list(compute = TRUE),
               control.fixed = control.fixed,
               control.compute = list(return.marginals.predictor = TRUE, 
                                      waic = TRUE, dic = TRUE))

summary(mod.T2)


#Type III interaction
formula.T3 <- y ~ f(ID.area, model = "bym", graph=graph.inla, 
                    scale.model = TRUE, constr = TRUE) +
  f(ID.year, model = "rw2", scale.model = TRUE, constr = TRUE) +
  f(ID.year1, model = "iid") +
  f(ID.year.int, model = "iid", group = ID.area.int, 
    control.group = list(model = "besag", graph=graph.inla), constr = TRUE)

mod.T3 <- inla(formula.T3,
               family = "poisson", data = Ohio_shp, E = E,
               control.predictor = list(compute = TRUE),
               control.fixed = control.fixed,
               control.compute = list(return.marginals.predictor = TRUE, 
                                      waic = TRUE, dic = TRUE))

summary(mod.T3)


#Type IV interaction
formula.T4 <- y ~ f(ID.area, model = "bym", graph=graph.inla, 
                    scale.model = TRUE, constr = TRUE) +
  f(ID.year, model = "rw2", scale.model = TRUE, constr = TRUE) +
  f(ID.year1, model = "iid", constr = TRUE) +
  f(ID.area.int, model = "besag", graph=graph.inla, 
    group = ID.year.int, control.group = list(model = "rw2"), 
    scale.model = TRUE, constr = TRUE)

mod.T4 <- inla(formula.T4,
               family = "poisson", data = Ohio_shp, E = E,
               control.predictor = list(compute = TRUE),
               control.fixed = control.fixed,
               control.compute = list(return.marginals.predictor = TRUE, 
                                      waic = TRUE, dic = TRUE))

summary(mod.T4)



#Choose a model and calculate the relative risks (RRs)
mod.dt$waic$waic
mod.T1$waic$waic
mod.T2$waic$waic
mod.T3$waic$waic
mod.T4$waic$waic


Ohio_shp$RR  <- mod.dt$summary.fitted.values[, "mean"]
Ohio_shp$low <- mod.dt$summary.fitted.values[, "0.025quant"]
Ohio_shp$up  <- mod.dt$summary.fitted.values[, "0.975quant"]

Ohio_shp$year_new <- Ohio_shp$year + 1967

#Plot SIR
pp1 <- ggplot(Ohio_shp) + geom_sf(aes(fill = SIR)) +
  facet_wrap(~year_new, dir = "h", ncol = 7) +
  ggtitle("SIR") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red", lim=c(0,2.6))


#Plot the estimated RRs
pp2 <- ggplot(Ohio_shp) + geom_sf(aes(fill = RR)) +
  facet_wrap(~year_new, dir = "h", ncol = 7) +
  ggtitle("RR") + theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  ) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red", lim=c(0,2.6))


#An animation using gganimate - SIR (modify for RR)
p <- ggplot(Ohio_shp) + geom_sf(aes(fill = SIR)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()) +
  scale_fill_gradient2(
    midpoint = 1, low = "blue", mid = "white", high = "red", lim=c(0,2.6)) +
  transition_time(year_new) +
  labs(title = "Year: {round(frame_time, 0)}")

# Render the animation
anim <- animate(p,
  nframes = length(unique(Ohio_shp$year_new)) * 5,
  fps = 5, width = 800, height = 600,
  renderer = gifski_renderer())

# Save as GIF
anim_save("Ohio_SIR_animation.gif", animation = anim)



