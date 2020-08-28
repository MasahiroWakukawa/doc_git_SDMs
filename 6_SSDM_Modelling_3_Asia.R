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

## Produce stacked analysis in Asia

Taxon_lim

Thinned_Data_5km_under100 <- c()

#5, 26, 30, 40, 41, 44, 46, 52

for (i in c(4, 6, 12, 15, 48, 52, 53)) {
  
  OC_test <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Thinned_Data_5km_under100 <- rbind(Thinned_Data_5km_under100, OC_test)
  
}

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(30 < decimallongitude & decimallongitude <175 &
                                                                    -15 < decimallatitude & decimallatitude < 89)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude > 45 | decimallatitude < 43)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude > 45 | decimallatitude > 20)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude < 125 | decimallatitude > 0)

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

SSDM_Asia_Presence_1 <- stack_modelling(c('MARS', 'CTA', 'RF'), Thinned_Data_5km_under100, Boundary_200km, rep=1, ensemble.thresh = 0,
                                             Xcol = 'decimallongitude', Ycol = 'decimallatitude', Spcol = 'species', method = "pSSDM", verbose = FALSE, tmp = TRUE)   

plot(SSDM_Asia_Presence_1@diversity.map, main = 'SSDM\nfor Plant\nwith MARS, CTA and RF algorithms')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_Asia_Presence_1.RData", SSDM_Asia_Presence_1)

##future prediction RCP2.6

buffer_205km <- raster::buffer(Coords, width = 205000)

Crop_205km_RCP26 <- crop(Bioclim_Stack_RCP26, extent(buffer_205km))

Boundary_205km_RCP26 <- stack(mask(Crop_205km_RCP26, buffer_205km))

SSDM_Asia_RCP26_1 <- SSDM::project(SSDM_Asia_Presence_1, Boundary_205km_RCP26)

plot(SSDM_Asia_RCP26_1@diversity.map, main = 'SSDM\nwith MARS, CTA and RF algorithms +1 temp and -20 prec')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_Asia_RCP26_1.RData", SSDM_Asia_RCP26_1)

##future prediction RCP8.5

Crop_205km_RCP85 <- crop(Bioclim_Stack_RCP85, extent(buffer_205km))

Boundary_205km_RCP85 <- stack(mask(Crop_205km_RCP85, buffer_205km))

SSDM_Asia_RCP85_1 <- SSDM::project(SSDM_Asia_Presence_1, Boundary_205km_RCP85)

plot(SSDM_Asia_RCP85_1@diversity.map, main = 'SSDM\nfor Plant\nwith MARS, CTA and RF algorithms +4 temp and -55 prec')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_Asia_RCP85_1.RData", SSDM_Asia_RCP85_1)


##Asia 2nd

Taxon_lim

Thinned_Data_5km_under100 <- c()

for (i in c(5, 26, 30, 40, 41, 44, 46, 52)) {
  
  OC_test <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Thinned_Data_5km_under100 <- rbind(Thinned_Data_5km_under100, OC_test)
  
}

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(30 < decimallongitude & decimallongitude <175 &
                                                                    -15 < decimallatitude & decimallatitude < 89)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude > 45 | decimallatitude < 43)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude > 60 | decimallatitude > 20)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude < 125 | decimallatitude > 0)

Thinned_Data_5km_under100 <- Thinned_Data_5km_under100 %>% filter(decimallongitude > 33 | decimallatitude > 35)

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

SSDM_Asia_Presence_2 <- stack_modelling(c('MARS', 'CTA', 'RF'), Thinned_Data_5km_under100, Boundary_200km, rep=1, ensemble.thresh = 0,
                                             Xcol = 'decimallongitude', Ycol = 'decimallatitude', Spcol = 'species', method = "pSSDM", verbose = FALSE, tmp = TRUE)   

plot(SSDM_Asia_Presence_2@diversity.map, main = 'SSDM\nfor Plant\nwith MARS, CTA and RF algorithms')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_Asia_Presence_2.RData", SSDM_Asia_Presence_2)

##future prediction RCP2.6

buffer_205km <- raster::buffer(Coords, width = 205000)

Crop_205km_RCP26 <- crop(Bioclim_Stack_RCP26, extent(buffer_205km))

Boundary_205km_RCP26 <- stack(mask(Crop_205km_RCP26, buffer_205km))

SSDM_Asia_RCP26_2 <- SSDM::project(SSDM_Asia_Presence_2, Boundary_205km_RCP26)

plot(SSDM_Asia_RCP26_2@diversity.map, main = 'SSDM\nwith MARS, CTA and RF algorithms +1 temp and -20 prec')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_Asia_RCP26_2.RData", SSDM_Asia_RCP26_2)

##future prediction RCP8.5

Crop_205km_RCP85 <- crop(Bioclim_Stack_RCP85, extent(buffer_205km))

Boundary_205km_RCP85 <- stack(mask(Crop_205km_RCP85, buffer_205km))

SSDM_Asia_RCP85_2 <- SSDM::project(SSDM_Asia_Presence_2, Boundary_205km_RCP85)

plot(SSDM_Asia_RCP85_2@diversity.map, main = 'SSDM\nfor Plant\nwith MARS, CTA and RF algorithms +4 temp and -55 prec')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

save(file = "../output/RData_Asia_RCP85_2.RData", SSDM_Asia_RCP85_2)


###stack all the Asia map

Asia_Presence <- mosaic(SSDM_Asia_Presence_1@diversity.map, SSDM_Asia_Presence_2@diversity.map, fun=sum)

plot(Asia_Presence, main = 'Asia Presence')

map("world", add=TRUE, lwd=0.6, interior = FALSE, col = 1)

ggsave(paste("../figs/","Stacked_analysis/Asia_Presence.png",sep = ""), plot(Asia_Presence, main = 'Asia Presence'), device = "png", height = 10, width = 10,dpi=300)

save(file = "../output/RData_Asia_Presence_all.RData", Asia_Presence)

Asia_RCP26 <- mosaic(SSDM_Asia_RCP26_1@diversity.map, SSDM_Asia_RCP26_2@diversity.map, fun=sum)

plot(Asia_RCP26, main = 'Asia RCP2.6')

ggsave(paste("../figs/","Stacked_analysis/Asia_RCP26.png",sep = ""), plot(Asia_RCP26, main = 'Asia RCP2.6'), device = "png", height = 10, width = 10,dpi=300)

save(file = "../output/RData_Asia_RCP26_all.RData", Asia_RCP26)

Asia_RCP85 <- mosaic(SSDM_Asia_RCP85_1@diversity.map, SSDM_Asia_RCP85_2@diversity.map, fun=sum)

plot(Asia_RCP85, main = 'Asia RCP8.5')

ggsave(paste("../figs/","Stacked_analysis/Asia_RCP85.png",sep = ""), plot(Asia_RCP85, main = 'Asia RCP8.5'), device = "png", height = 10, width = 10,dpi=300)

save(file = "../output/RData_Asia_RCP85_all.RData", Asia_RCP85)


###calculate all the Asia map

load(file = "../output/RData_Asia_Presence_all.RData")

load(file = "../output/RData_Asia_RCP26_all.RData")

load(file = "../output/RData_Asia_RCP85_all.RData")

Asia_Presence <- readAll(Asia_Presence)

Asia_RCP26 <- readAll(Asia_RCP26)

Asia_RCP85 <- readAll(Asia_RCP85)

Presence_value <- sum(Asia_Presence@data@values)/length(Asia_Presence@data@values)

RCP26_value <- sum(Asia_RCP26@data@values)/length(Asia_RCP26@data@values)

RCP85_value <- sum(Asia_RCP85@data@values)/length(Asia_RCP85@data@values)

Asia_Scenarios <- cbind(log(RCP26_value/Presence_value), log(RCP85_value/Presence_value))

colnames(Asia_Scenarios) <- c("log(RCP2.6/Presence)", "log(RCP8.5/Presence)")

write.table(file = paste("../figs/","Stacked_analysis/Asia_Evaluation",sep=""), Asia_Scenarios)

