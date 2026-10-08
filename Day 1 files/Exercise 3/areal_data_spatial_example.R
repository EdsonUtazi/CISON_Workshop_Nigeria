
#Spatial areal data analysis example

library(SpatialEpi)
library(sf)
library(dplyr)
library(leaflet)
library(ggplot2)
library(mapview)
library(servr)
library(spdep)
library(INLA)
library(patchwork)


#To see all the data in the SpatialEpi package
data(package = "SpatialEpi")

#Load the Scotland lip cancer data
data(scotland)

#See detailed description of the data
?scotland

#Explore the data attributes/components
names(scotland)

#View data to be modelled
#AFF is the proportion of the population engaged in agriculture, fishing or 
#forestry
View(scotland$data)


#View the map of the counties in Scotland saved in the spatial.polygon component
#of the data file
shp_scot <- scotland$spatial.polygon
plot(shp_scot)

#Check number of areas/counties in the shapefile
nrow(shp_scot)


#Check the CRS (coordinate reference system) of the shapefile
#(contains none)
st_crs(shp_scot)$Name

#Assign the following CRS to the shapefile 
#This is the original CRS for the data
#We use proj4string() because the polygon is an sp object
#Otherwise, we could use crs() for terra objects and st_crs() for sf objects
#Observe that CRS is the projected CRS 

proj4string(shp_scot) <- "+proj=tmerc +lat_0=49 +lon_0=-2
+k=0.9996012717 +x_0=400000 +y_0=-100000 +datum=OSGB36
+units=km +no_defs"


#To use the leaflet R package to create maps with the shapefile, we will convert
#it to the geographic CRS
shp_scot <- st_as_sf(shp_scot) #convert to an sf object
shp_scot <- st_transform(shp_scot, crs = 4326) #Transform to lat/lon coords


#Data frame for modelling
dat <- scotland$data[,c("county.names", "cases", "expected", "AFF")]
names(dat) <- c("county", "Y", "E", "AFF")


#Calculate the standardized incidence ratio (SIR)
dat$SIR <- dat$Y / dat$E

#Add the data frame to the shapefile
shp_scot <- dplyr::bind_cols(shp_scot, dat)
head(shp_scot)

#Make a static map of SIR using ggplot2
#First create up to 4 SIR classes for a discrete legend
shp_scot <- shp_scot %>%
  mutate(
    SIR_class = cut(
      SIR,
      breaks = 4,
      include.lowest = TRUE
    )
  )

# Static map
map.static <- ggplot(shp_scot) +
  geom_sf(aes(fill = SIR_class), colour = "grey50", linewidth = 0.2) +
  scale_fill_viridis_d(
    name = "SIR",
    option = "viridis",
    direction = 1
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "")


map.static

#You can save this map using ggsave


#make an interactive map using the leaflet package
map.interac <- leaflet(shp_scot) %>% addTiles()

pal <- colorNumeric(palette = "YlOrRd", domain = shp_scot$SIR)

map.interac <- map.interac %>%
  addPolygons(
    color = "grey", weight = 1,
    fillColor = ~ pal(SIR), fillOpacity = 0.5
  ) %>%
  addLegend(
    pal = pal, values = ~SIR, opacity = 0.5,
    title = "SIR", position = "bottomright"
  )

#Display the interactive map (find in the Viewer tab)
map.interac


#To share the interactive map or open it in a browser, you need to save it to a folder
#Path to folder
file.path <- "C:/Users/ceu1c14/OneDrive - University of Southampton/Documents/RShiny_practice/Interac_maps/"

#Save it as an html widget to the folder
htmlwidgets::saveWidget(
  map.interac,
  paste0(file.path, "SIR_map.html")
)

#Then run the function below to to generate a local address to view the file in 
#R console
servr::httd(file.path)


#Copy the local address and paste it in your browser, then click the file 
#to open it


#We can use map view as well
#Make sure to select an appropriate basemap 
#You can also export the interactive file using the RStudio Export button 
mapview::mapview(shp_scot, zcol = "SIR") #NOTE: Change labels

#Model fitting
#We will fit and compare different CAR models using the INLA functions provided
#below.

#Create a neighbourhood list using poly2nb
nb.scot <- poly2nb(shp_scot)

#Convert to a file to be used by INLA
nb2INLA("graph.scot.adj", nb.scot)
graph.inla <- inla.read.graph(filename = "graph.scot.adj")

#Create ID variables for the areas in the shapefile
shp_scot$idarea1 <- 1:nrow(shp_scot)
shp_scot$idarea2 <- 1:nrow(shp_scot)

#Set the prior distributions
prior_prec <- list(theta = list(prior="pc.prec", param=c(3, 0.01)))
control.fixed <- list(mean=0, prec=1/1000, mean.intercept=0, 
                      prec.intercept=1/1000) #for regression coefficients

#Model formula for BYM
formula1 <- Y ~ AFF +
  f(idarea1, model = "besag", graph = graph.inla, scale.model = TRUE, 
    hyper = prior_prec) +
  f(idarea2, model = "iid") #, hyper = prior_prec

#Fit BYM model
res1 <- inla(formula1,
            family = "poisson", data = shp_scot, E = E,
            control.predictor = list(compute = TRUE),
            control.fixed = control.fixed,
            control.compute = list(return.marginals.predictor = TRUE, 
                                   waic = TRUE, dic = TRUE))

summary(res1)

#Model formula for BYM2
hyper.bym2 <- list(prec = list(prior = "pc.prec", param = c(3, 0.01)),
  phi = list(prior = "pc", param = c(0.5, 0.9)))

formula2 <- Y ~ AFF +
  f(idarea1, model = "bym2", graph = graph.inla,
    hyper = hyper.bym2)

#Fit BYM2 model
res2 <- inla(formula2,
             family = "poisson", data = shp_scot, E = E,
             control.predictor = list(compute = TRUE),
             control.fixed = control.fixed,
             control.compute = list(return.marginals.predictor = TRUE, 
                                    waic = TRUE, dic = TRUE))

summary(res2)



#Model formula for Leroux model
#A default prior is used for the spatial autocorrelation parameter
hyper.leroux <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))

formula3 <- Y ~ AFF +
  f(idarea1, model = "besagproper2", graph = graph.inla,
    hyper = hyper.leroux)

#Fit Leroux model
res3 <- inla(formula3,
             family = "poisson", data = shp_scot, E = E,
             control.predictor = list(compute = TRUE),
             control.fixed = control.fixed,
             control.compute = list(return.marginals.predictor = TRUE, 
                                    waic = TRUE, dic = TRUE))

summary(res3)


#Compare the WAIC values of all three models
res1$waic$waic
res2$waic$waic
res3$waic$waic


#Compare the dic values of all three models
res1$dic$dic
res2$dic$dic
res3$dic$dic


#Choose a model and use this to calculate the relative risks

shp_scot$RR <- res2$summary.fitted.values[, "mean"]
shp_scot$sd <- res2$summary.fitted.values[, "sd"]
shp_scot$lower <- res2$summary.fitted.values[, "0.025quant"]
shp_scot$upper <- res2$summary.fitted.values[, "0.975quant"]


#Plot the RR estimates and associated uncertainties using ggplot2
map.RR <- ggplot(shp_scot) +
  geom_sf(aes(fill = RR), colour = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "RR",
    option = "viridis",
    limits = c(0, 7),
    direction = 1
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "")

map.sd <- ggplot(shp_scot) +
  geom_sf(aes(fill = sd), colour = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "SD",
    option = "magma",
    direction = 1
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "")

#Continuous scale map for SIR for comparison
map.SIR <- ggplot(shp_scot) +
  geom_sf(aes(fill = SIR), colour = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "SIR",
    option = "viridis",
    limits = c(0, 7),
    direction = 1
  ) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  labs(title = "")

combined_map <- map.SIR + map.RR + map.sd +
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "a")

combined_map

#Save map
ggsave(filename = "combined_maps.png",
  plot = combined_map, width = 15,
  height = 5, units = "in", dpi = 300)

#Q: Between the SIR map and the RR map, which looks smoother?


#You can also calculate and plot the exceedance probabilities p(RR>1)
prob_RR_gt1 <- sapply(
  res2$marginals.fitted.values,
  function(m) {
    1 - inla.pmarginal(q=1, m)
  }
)

#An alternative function
#prob_RR_gt1 <- sapply(
#  res2$marginals.linear.predictor,
#  function(m) {
#    1 - inla.pmarginal(0, m)
#  }
#)

#Add to shapefile
shp_scot$prob_RR_gt1 <- prob_RR_gt1


map.EP <- ggplot(shp_scot) +
  geom_sf(aes(fill = prob_RR_gt1), colour = "grey50", linewidth = 0.2) +
  scale_fill_viridis_c(
    name = "P(RR > 1)",
    option = "magma",
    limits = c(0, 1)
  ) +
  theme_void()


#RR vs EP plot
combined_map2 <- map.RR + map.EP+
  plot_layout(ncol = 3) +
  plot_annotation(tag_levels = "a")

combined_map2



