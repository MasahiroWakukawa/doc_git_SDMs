rm(list= ls())
#require(biomod2)
require(dismo)
require(raster)
require(magrittr)
require(foreach)
#require(doParallel)
require(ff)
require(SSDM)
require(dplyr)
require(ggplot2)
require(rgeos)
require(maptools)
library(knitr)
library(kableExtra)
library(rJava)
library(maps)



load("../output/RData_Plant_Species_Thinned_Data_5km.RData")     

Thinned_Data_5km <- Thinned_Data_5km %>% filter(-175 < decimallongitude & decimallongitude <175 &
                                                  -89 < decimallatitude & decimallatitude < 89)

Taxon_lim <- Thinned_Data_5km %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

Taxon_lim

#Bioclim 1 : Annual mean temperature
#Bioclim 4 : temperature seasonality
#Bioclim 12 : annual precipitation
#Bioclim 15 : Precipitation seasonality

Bioclim_01 <- raster("../data/wc2/wc2.1_2.5m_bio_1.tif")

Bioclim_04 <- raster("../data/wc2/wc2.1_2.5m_bio_4.tif")

Bioclim_12 <- raster("../data/wc2/wc2.1_2.5m_bio_12.tif")

Bioclim_15 <- raster("../data/wc2/wc2.1_2.5m_bio_15.tif")   

Bioclim_Stack <- stack(Bioclim_01,Bioclim_04,Bioclim_12,Bioclim_15)

Bioclim_01_RCP26 <- raster("../data/bc26bi70/wc2.1_2.5m_bio_1.tif")

Bioclim_04_RCP26 <- raster("../data/bc26bi70/wc2.1_2.5m_bio_4.tif")

Bioclim_12_RCP26 <- raster("../data/bc26bi70/wc2.1_2.5m_bio_12.tif")

Bioclim_15_RCP26 <- raster("../data/bc26bi70/wc2.1_2.5m_bio_15.tif")

Bioclim_Stack_RCP26 <- stack(Bioclim_01_RCP26,Bioclim_04_RCP26,Bioclim_12_RCP26,Bioclim_15_RCP26)

Bioclim_01_RCP85 <- raster("../data/bc85bi70/wc2.1_2.5m_bio_1.tif")

Bioclim_04_RCP85 <- raster("../data/bc85bi70/wc2.1_2.5m_bio_4.tif")

Bioclim_12_RCP85 <- raster("../data/bc85bi70/wc2.1_2.5m_bio_12.tif")

Bioclim_15_RCP85 <- raster("../data/bc85bi70/wc2.1_2.5m_bio_15.tif")

Bioclim_Stack_RCP85 <- stack(Bioclim_01_RCP85,Bioclim_04_RCP85,Bioclim_12_RCP85,Bioclim_15_RCP85)

Occurrence_Data_5km <- Thinned_Data_5km[c("species", "decimallongitude", "decimallatitude")]

Occurrence_Data_5km$species <- as.character(Occurrence_Data_5km$species)

## Produce stacked analysis in South America

Taxon_lim

Thinned_Data_5km_under100 <- c()

for (i in c(2, 3, 17, 18, 23, 25, 26, 39, 40, 41, 44, 45, 46, 49, 52)) {
  
  OC_test <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Thinned_Data_5km_under100 <- rbind(Thinned_Data_5km_under100, OC_test)
  
}

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(-90 < decimallongitude & decimallongitude < -30 &
                                                                    -60 < decimallatitude & decimallatitude < 15)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude > -80 | decimallatitude < 5)

wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify()

ggplot()+ coord_fixed()+geom_map(data =wm, map = wm,
                                 
                                 aes(group = group, map_id= region),
                                 
                                 fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
  
  geom_point(data = Thinned_Data_5km_under100, aes(x = decimallongitude, y = decimallatitude),
             
             colour = "darkred", size = 1)

Coords <- Thinned_Data_5km_under100[,c("decimallongitude","decimallatitude")]

Coords <- SpatialPointsDataFrame(coords = Coords, data = Thinned_Data_5km_under100,
                                 proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

buffer_200km <- raster::buffer(Coords, width = 200000)

Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))

Boundary_200km <- stack(mask(Crop_200km, buffer_200km))

SSDM_SouthAmerica_Presence <- stack_modelling(c('MARS', 'CTA', 'RF'), Thinned_Data_5km_under100, Boundary_200km, rep=1, ensemble.thresh = 0,
                                        Xcol = 'decimallongitude', Ycol = 'decimallatitude', Spcol = 'species', method = "pSSDM", verbose = FALSE, tmp = TRUE)   

plot(SSDM_SouthAmerica_Presence@diversity.map, main = 'SSDM\nfor Plant\nwith MARS, CTA and RF algorithms')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_SouthAmerica_Presence.RData", SSDM_SouthAmerica_Presence)

##future prediction RCP2.6

buffer_205km <- raster::buffer(Coords, width = 205000)

Crop_205km_RCP26 <- crop(Bioclim_Stack_RCP26, extent(buffer_205km))

Boundary_205km_RCP26 <- stack(mask(Crop_205km_RCP26, buffer_205km))

SSDM_SouthAmerica_RCP26 <- SSDM::project(SSDM_SouthAmerica_Presence, Boundary_205km_RCP26)

plot(SSDM_SouthAmerica_RCP26@diversity.map, main = 'SSDM\nwith MARS, CTA and RF algorithms +1 temp and -20 prec')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_SouthAmerica_RCP26.RData", SSDM_SouthAmerica_RCP26)

##future prediction RCP8.5

Crop_205km_RCP85 <- crop(Bioclim_Stack_RCP85, extent(buffer_205km))

Boundary_205km_RCP85 <- stack(mask(Crop_205km_RCP85, buffer_205km))

SSDM_SouthAmerica_RCP85 <- SSDM::project(SSDM_SouthAmerica_Presence, Boundary_205km_RCP85)

plot(SSDM_SouthAmerica_RCP85@diversity.map, main = 'SSDM\nfor Plant\nwith MARS, CTA and RF algorithms +4 temp and -55 prec')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_SouthAmerica_RCP85.RData", SSDM_SouthAmerica_RCP85)

###calculate all the SouthAmerica map

load(file = "../output/RData_SouthAmerica_Presence.RData")

load(file = "../output/RData_SouthAmerica_RCP26.RData")

load(file = "../output/RData_SouthAmerica_RCP85.RData")

plot(SSDM_SouthAmerica_Presence@diversity.map, main = 'SouthAmerica Presence')

ggsave(paste("../figs/","Stacked_analysis/SouthAmerica_Presence.png",sep = ""), plot(SSDM_SouthAmerica_Presence@diversity.map, main = 'SouthAmerica Presence'), device = "png", height = 10, width = 10,dpi=300)

plot(SSDM_SouthAmerica_RCP26@diversity.map, main = 'SouthAmerica RCP2.6')

ggsave(paste("../figs/","Stacked_analysis/SouthAmerica_RCP26.png",sep = ""), plot(SSDM_SouthAmerica_RCP26@diversity.map, main = 'SouthAmerica RCP2.6'), device = "png", height = 10, width = 10,dpi=300)

plot(SSDM_SouthAmerica_RCP85@diversity.map, main = 'SouthAmerica RCP8.5')

ggsave(paste("../figs/","Stacked_analysis/SouthAmerica_RCP85.png",sep = ""), plot(SSDM_SouthAmerica_RCP85@diversity.map, main = 'SouthAmerica RCP8.5'), device = "png", height = 10, width = 10,dpi=300)

SSDM_SouthAmerica_Presence@diversity.map <- readAll(SSDM_SouthAmerica_Presence@diversity.map)

SSDM_SouthAmerica_RCP26@diversity.map <- readAll(SSDM_SouthAmerica_RCP26@diversity.map)

SSDM_SouthAmerica_RCP85@diversity.map <- readAll(SSDM_SouthAmerica_RCP85@diversity.map)

Presence_value <- sum(SSDM_SouthAmerica_Presence@diversity.map@data@values, na.rm = TRUE)/length(SSDM_SouthAmerica_Presence@diversity.map@data@values)

RCP26_value <- sum(SSDM_SouthAmerica_RCP26@diversity.map@data@values, na.rm = TRUE)/length(SSDM_SouthAmerica_RCP26@diversity.map@data@values)

RCP85_value <- sum(SSDM_SouthAmerica_RCP85@diversity.map@data@values, na.rm = TRUE)/length(SSDM_SouthAmerica_RCP85@diversity.map@data@values)

SouthAmerica_Scenarios <- cbind(log(RCP26_value/Presence_value), log(RCP85_value/Presence_value))

colnames(SouthAmerica_Scenarios) <- c("log(RCP2.6/Presence)", "log(RCP8.5/Presence)")

write.table(file = paste("../figs/","Stacked_analysis/SouthAmerica_Evaluation",sep=""), SouthAmerica_Scenarios)

