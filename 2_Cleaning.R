rm(list=ls())
require(dplyr)
require(tidyr)
require(rgbif)
require(magrittr)
require(countrycode)
require(CoordinateCleaner)
require(rnaturalearthdata)
require(rgdal)
require(sp)
require(ggplot2)
require(adehabitatHR)
require(scales)
require(doParallel)
require(sf)
require(raster)



load("../output/RData_Plant_Species_Data.RData") 


### Cleaning the data by clean_coordinates function
dat <- Plant

dat <- dat %>% filter(!is.na(decimallongitude)) %>% filter(!is.na(decimallatitude))

dat <- dat %>% filter(-180 < decimallongitude & decimallongitude <180 &
                                     -90 < decimallatitude & decimallatitude < 90)

flags_dat <- clean_coordinates(x = dat, lon = "decimallongitude", lat = "decimallatitude",
                                  species = "species",
                                  tests = c("capitals", "centroids", "duplicates", "equal", "gbif", "institutions",
                                            "seas", "zeros"),   
                                  seas_scale = 110)

Clean_Data <- flags_dat

Flagged_Data <- Clean_Data %>% filter(.summary == FALSE)

Clean_Data <- Clean_Data %>% filter(.summary == TRUE)

head(Clean_Data)


### Cleaning by count more than 20

Clean_Data <- Clean_Data[ , 1:6] 

#Taxon_lim <- Clean_Data %>% group_by(species, taxonkey) %>% summarise(COUNT=n()) 

Taxon_lim <- Clean_Data %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- Taxon_lim %>% filter(COUNT >20)

Clean_Data <- Clean_Data %>% filter(Clean_Data$species%in% Taxon_lim$species)

Taxon_lim <- as.data.frame(Taxon_lim)

Taxon_lim

CC_Data <- c()
for(i in 1:nrow(Taxon_lim)){
  GBIF <- Clean_Data %>% filter(species == Taxon_lim[i, 1])
  GBIF <- GBIF %>% cc_outl(mltpl = 12)
  CC_Data <- rbind(CC_Data, GBIF)
  print(i)
}

Clean_Data <- CC_Data

save(Clean_Data, file = "../output/RData_Plant_Species_Clean_Data.RData")


### Downloading shape file and cropping by that areas

library(rgdal)
library(raster)
library(ggplot2)
library(rgeos)
library(mapview)
library(leaflet)
library(broom) 
library(sp)
library(sf)
options(stringsAsFactors = FALSE)

# making vector file using sp package
plant_shape <- readOGR("../data/level3/level3.shp")

class(plant_shape)

extent(plant_shape)

crs(plant_shape)

plot(plant_shape, main = "Plant Shape")

plant_shape_simp <- gSimplify(plant_shape, 
                              tol = .1,             ### tol = 3 : more coarse one
                              topologyPreserve = TRUE)

plot(plant_shape_simp,
     main = "Plant Shape with Boundaries Simplified")

ggplot() +
  geom_path(data = plant_shape_simp, aes(x = long, y = lat, group = group)) +
  labs(title = "Plant Shape - using ggplot")

ggplot() +
  geom_path(data = plant_shape, aes(x = long, y = lat, group = group)) + 
  coord_fixed() + 
  labs(title = "My awesome ggplot map of coastlines",
       subtitle = "my awesome subtitle",
       x = "", y = "") 

mapview(plant_shape_simp)

leaflet(plant_shape_simp) %>%
  addTiles() %>% 
  addPolylines(color = "#444444", weight = 1, smoothFactor = 0.5,
               opacity = 1.0)

# making vector file using sf function
plant_shape <- st_read("../data/level3/level3.shp")

st_geometry_type(plant_shape)

st_crs(plant_shape)

st_bbox(plant_shape)

plant_shape

plot(plant_shape[2])

ggplot() +
  geom_sf(data=plant_shape, aes(fill = COUNT))  

# making raster file converted from vector file

ext <- extent(-180, 180, -90, 84)

xy <- abs(apply(as.matrix(bbox(ext)), 1, diff))

n <- 5

r <- raster(ext, ncol=xy[1]*n, nrow=xy[2]*n)

plant_shape_ras <- rasterize(plant_shape, r)

plot(plant_shape_ras)

plant_shape@data
