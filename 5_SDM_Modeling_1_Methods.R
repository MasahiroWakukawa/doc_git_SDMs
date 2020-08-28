rm(list= ls())
#require(biomod2)
require(dismo)
require(raster)
require(magrittr)
require(foreach)
require(doParallel)
require(ff)
require(SSDM)
require(dplyr)
require(ggplot2)
require(rgeos)
require(maptools)
library(knitr)
library(kableExtra)


load("../output/RData_Plant_Species_Thinned_Data_5km.RData")     ### Load Presence and Absence Data

Thinned_Data_5km <- Thinned_Data_5km %>% filter(-175 < decimallongitude & decimallongitude <175 &
                                                  -89 < decimallatitude & decimallatitude < 89)

Taxon_lim <- Thinned_Data_5km %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

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

### SDMs by SSDM package


Taxon_lim


##Using GLM method

for (i in 1:nrow(Taxon_lim)){
  
  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  dir.create(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", Taxon_lim[i,1],sep = ""))
  
  ### Creating GLM models
  
  SDM_GLM_200km<- modelling('GLM', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "GLM_200km", sep = "_"), SDM_GLM_200km)
  
  plot(SDM_GLM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith GLM 200km',sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GLM_200km.png",sep=""),sep = ""), plot = plot(SDM_GLM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith GLM 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  Bioclim_Stack_new <- stack(Bioclim_01+1.0, Bioclim_04, Bioclim_12-20, Bioclim_15)
  
  buffer_208km <- raster::buffer(Coords, width = 208000)
  
  Crop_208km_new <- crop(Bioclim_Stack_new, extent(buffer_208km))
  
  Boundary_208km_new <- stack(mask(Crop_208km_new, buffer_208km))
  
  SDM_projection_208km <- SSDM::project(SDM_GLM_200km, Boundary_208km_new)
  
  plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GLM 208km +1 temp and -20 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/GLM_208km_RCP2_6.png",sep=""),sep = ""), plot = plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GLM 208km +1 temp and -20 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Bioclim_Stack_new_2 <- stack(Bioclim_01+4.0, Bioclim_04, Bioclim_12-55, Bioclim_15)
  
  Crop_208km_new_2 <- crop(Bioclim_Stack_new_2, extent(buffer_208km))
  
  Boundary_208km_new_2 <- stack(mask(Crop_208km_new_2, buffer_208km))
  
  SDM_projection_2_208km <- SSDM::project(SDM_GLM_200km, Boundary_208km_new_2)
  
  plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GLM 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/GLM_208km_RCP8_5.png",sep=""),sep = ""), plot = plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GLM 208km +4 temp and -55 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating evaluation table
  
  xy_presence <- SDM_GLM_200km@data[, 1:3]
  
  xy <- xy_presence[,c("X","Y")]
  
  sp_coord <- SpatialPointsDataFrame(coords = xy, data = xy_presence,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  existence_probability <- raster::extract(SDM_GLM_200km@projection, sp_coord, cellnumbers=TRUE)
  
  accuracy <- mean(abs(xy_presence[,3]-existence_probability[,2]))
  
  discrimination <- SDM_GLM_200km@evaluation[, "AUC"]
  
  calibration <- abs(sum(existence_probability[, "Projection"])/nrow(existence_probability)-sum(xy_presence[, "Presence"])/nrow(xy_presence))
  
  precision <- sum(sqrt(existence_probability[, "Projection"]*(1-existence_probability[, "Projection"])))/nrow(existence_probability)

  evaluate_table <- rbind(c(accuracy, discrimination, calibration, precision))
  
  colnames(evaluate_table) <- c("accuracy", "discrimination", "calibration", "precision")
  
  write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GLM_table",sep=""),sep = ""), evaluate_table)
  
  print(i)
  
}

rm(list = ls(pattern = "GLM*"))


##Using GAM method

for (i in 1:nrow(Taxon_lim)){
  
  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  
  ### Creating GAM models
  
  SDM_GAM_200km<- modelling('GAM', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "GAM_200km", sep = "_"), SDM_GAM_200km)
  
  plot(SDM_GAM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith GAM 200km',sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GAM_200km.png",sep=""),sep = ""), plot = plot(SDM_GAM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith GAM 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  Bioclim_Stack_new <- stack(Bioclim_01+1.0, Bioclim_04, Bioclim_12-20, Bioclim_15)
  
  buffer_208km <- raster::buffer(Coords, width = 208000)
  
  Crop_208km_new <- crop(Bioclim_Stack_new, extent(buffer_208km))
  
  Boundary_208km_new <- stack(mask(Crop_208km_new, buffer_208km))
  
  SDM_projection_208km <- SSDM::project(SDM_GAM_200km, Boundary_208km_new)
  
  plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GAM 208km +1 temp and -20 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/GAM_208km_RCP2_6.png",sep=""),sep = ""), plot = plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GAM 208km +1 temp and -20 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Bioclim_Stack_new_2 <- stack(Bioclim_01+4.0, Bioclim_04, Bioclim_12-55, Bioclim_15)
  
  Crop_208km_new_2 <- crop(Bioclim_Stack_new_2, extent(buffer_208km))
  
  Boundary_208km_new_2 <- stack(mask(Crop_208km_new_2, buffer_208km))
  
  SDM_projection_2_208km <- SSDM::project(SDM_GAM_200km, Boundary_208km_new_2)
  
  plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GAM 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/GAM_208km_RCP8_5.png",sep=""),sep = ""), plot = plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GAM 208km +4 temp and -55 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating evaluation table
  
  xy_presence <- SDM_GAM_200km@data[, 1:3]
  
  xy <- xy_presence[,c("X","Y")]
  
  sp_coord <- SpatialPointsDataFrame(coords = xy, data = xy_presence,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  existence_probability <- raster::extract(SDM_GAM_200km@projection, sp_coord, cellnumbers=TRUE)
  
  accuracy <- mean(abs(xy_presence[,3]-existence_probability[,2]))
  
  discrimination <- SDM_GAM_200km@evaluation[, "AUC"]
  
  calibration <- abs(sum(existence_probability[, "Projection"])/nrow(existence_probability)-sum(xy_presence[, "Presence"])/nrow(xy_presence))
  
  precision <- sum(sqrt(existence_probability[, "Projection"]*(1-existence_probability[, "Projection"])))/nrow(existence_probability)
  
  evaluate_table <- rbind(c(accuracy, discrimination, calibration, precision))
  
  colnames(evaluate_table) <- c("accuracy", "discrimination", "calibration", "precision")
  
  write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GAM_table",sep=""),sep = ""), evaluate_table)
  
  print(i)
}

rm(list = ls(pattern = "GAM*"))


##Using MARS method

for (i in 1:nrow(Taxon_lim)){

  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  dir.create(paste("../figs/", Taxon_lim[i,1],sep = ""))
  
  ### Creating MARS models
  
  SDM_MARS_200km<- modelling('MARS', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "MARS_200km", sep = "_"), SDM_MARS_200km)
  
  plot(SDM_MARS_200km@projection, main = paste(Taxon_lim[i,1], '\nwith MARS 200km',sep = ""))
  
  ggsave(paste("../figs/",paste(Taxon_lim[i,1],"/MARS_200km.png",sep=""),sep = ""), plot = plot(SDM_MARS_200km@projection, main = paste(Taxon_lim[i,1], '\nwith MARS 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  buffer_205km <- raster::buffer(Coords, width = 205000)
  
  Crop_205km_new <- crop(Bioclim_Stack_RCP26, extent(buffer_205km))
  
  Boundary_205km_new <- stack(mask(Crop_205km_new, buffer_205km))
  
  SDM_projection_205km <- SSDM::project(SDM_MARS_200km, Boundary_205km_new)
  
  plot(SDM_projection_205km@projection, main = paste(Taxon_lim[i,1], '\nwith MARS 205km +1 temp and -20 prec', sep = ""))
  
  ggsave(paste("../figs/", paste(Taxon_lim[i,1],"/MARS_205km_RCP26.png",sep=""),sep = ""), plot = plot(SDM_projection_205km@projection, main = paste(Taxon_lim[i,1], '\nwith MARS 205km RCP2.6', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Crop_205km_new_2 <- crop(Bioclim_Stack_RCP85, extent(buffer_205km))
  
  Boundary_205km_new_2 <- stack(mask(Crop_205km_new_2, buffer_205km))
  
  SDM_projection_2_205km <- SSDM::project(SDM_MARS_200km, Boundary_205km_new_2)
  
  plot(SDM_projection_2_205km@projection, main = paste(Taxon_lim[i,1], '\nwith MARS 205km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("../figs/", paste(Taxon_lim[i,1],"/MARS_205km_RCP85.png",sep=""),sep = ""), plot = plot(SDM_projection_2_205km@projection, main = paste(Taxon_lim[i,1], '\nwith MARS 205km RCP8.5', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Compare Presence vs RCP2.6 and RCP8.5
  
  SDM_MARS_200km@projection <- readAll(SDM_MARS_200km@projection)
  
  SDM_projection_205km@projection <- readAll(SDM_projection_205km@projection)
  
  SDM_projection_2_205km@projection <- readAll(SDM_projection_2_205km@projection)
  
  Presence_value <- sum(SDM_MARS_200km@projection@data@values, na.rm = TRUE)/length(SDM_MARS_200km@projection@data@values)
  
  RCP26_value <- sum(SDM_projection_205km@projection@data@values, na.rm = TRUE)/length(SDM_projection_205km@projection@data@values)
  
  RCP85_value <- sum(SDM_projection_2_205km@projection@data@values, na.rm = TRUE)/length(SDM_projection_2_205km@projection@data@values)
  
  Scenarios <- cbind(log(RCP26_value/Presence_value), log(RCP85_value/Presence_value))
  
  colnames(Scenarios) <- c("log(RCP2.6/Presence)", "log(RCP8.5/Presence)")
  
  write.table(file = paste("../figs/", paste(Taxon_lim[i,1],"/Scenarios_MARS",sep=""),sep = ""), Scenarios)
  
  print(i)
}

rm(list = ls(pattern = "MARS*"))


### Creating evaluation table

xy_presence <- SDM_MARS_200km@data[, 1:3]

xy <- xy_presence[,c("X","Y")]

sp_coord <- SpatialPointsDataFrame(coords = xy, data = xy_presence,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

existence_probability <- raster::extract(SDM_MARS_200km@projection, sp_coord, cellnumbers=TRUE)

accuracy <- mean(abs(xy_presence[,3]-existence_probability[,2]))

discrimination <- SDM_MARS_200km@evaluation[, "AUC"]

calibration <- abs(sum(existence_probability[, "Projection"])/nrow(existence_probability)-sum(xy_presence[, "Presence"])/nrow(xy_presence))

precision <- sum(sqrt(existence_probability[, "Projection"]*(1-existence_probability[, "Projection"])))/nrow(existence_probability)

evaluate_table <- rbind(c(accuracy, discrimination, calibration, precision))

colnames(evaluate_table) <- c("accuracy", "discrimination", "calibration", "precision")

write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/MARS_table",sep=""),sep = ""), evaluate_table)


### ESDM modelling

ESDM <- ensemble_modelling(c('CTA', 'MARS'), subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[1, 1]), 
                           ExplVar_10x, rep=1, Xcol = 'decimallongitude', Ycol = 'decimallatitude', ensemble.thresh = 0, verbose = FALSE)

```{r, echo=FALSE}
plot(ESDM@projection, main = 'ESDM\nfor Acianthus elegans\nwith CTA and MARS algorithms')
```



### SSDMs by SSDM package

Plant1_3 <- c()

for (i in c(1,3)) {
  OC_test <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Plant1_3 <- rbind(Plant1_3, OC_test)
  
  }

xmin_10x <- floor(min(Plant1_3$decimallongitude))                             #### Get the extent of the SDMs

xmax_10x <- ceiling(max(Plant1_3$decimallongitude))

ymin_10x <- floor(min(Plant1_3$decimallatitude))

ymax_10x <- ceiling(max(Plant1_3$decimallatitude))

extent_10x <- extent(xmin_10x,xmax_10x,ymin_10x,ymax_10x)

ExplVar_10x<- stack(crop(Bioclim_Stack, extent_10x)) 

Stack_SDM <- stack_modelling(c('CTA', 'SVM'), Plant1_3, ExplVar_10x, rep=1, ensemble.thresh = 0,
                        Xcol = 'decimallongitude', Ycol = 'decimallatitude', Spcol = 'species', method = "pSSDM", verbose = FALSE, tmp = TRUE)   

plot(Stack_SDM@diversity.map, main = 'SSDM\nfor Plant\nwith CTA and SVM algorithms')


knitr::kable(SDM@evaluation)

knitr::kable(ESDM@evaluation)

Bioclim_Stack_new <- stack(Bioclim_01+1.0, Bioclim_04, Bioclim_12-100, Bioclim_15)

SDM_projection <- project(SDM, Bioclim_Stack_new)

plot(SDM_projection@projection, main = 'SDM\nfor Acianthus elegans\nwith GLM algorithm +1 temp and -100 prec')


### SDM for all plant species
xmin_10x <- floor(min(Occurrence_Data_5km$decimallongitude))                             #### Get the extent of the SDMs

xmax_10x <- ceiling(max(Occurrence_Data_5km$decimallongitude))

ymin_10x <- floor(min(Occurrence_Data_5km$decimallatitude))

ymax_10x <- ceiling(max(Occurrence_Data_5km$decimallatitude))

extent_10x <- extent(xmin_10x,xmax_10x,ymin_10x,ymax_10x)

ExplVar_10x<- stack(crop(Bioclim_Stack, extent_10x)) 

Stack_SDM <- stack_modelling('GLM', Occurrence_Data_5km, ExplVar_10x, rep=1, ensemble.thresh = 0,
                             Xcol = 'decimallongitude', Ycol = 'decimallatitude', Spcol = 'species', method = "pSSDM", verbose = FALSE, tmp = TRUE)   


