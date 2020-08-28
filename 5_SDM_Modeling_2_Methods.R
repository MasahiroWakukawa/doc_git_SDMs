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

load("../output/RData_Plant_Species_Thinned_Data_5km.RData")

Thinned_Data_5km <- Thinned_Data_5km %>% filter(-175 < decimallongitude & decimallongitude <175 &
                                                  -89 < decimallatitude & decimallatitude < 89)


Taxon_lim <- Thinned_Data_5km %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

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

Taxon_lim

##Using GBM method

for (i in 45:nrow(Taxon_lim)){
  
  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  
  ### Creating GBM models
  
  SDM_GBM_200km<- modelling('GBM', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "GBM_200km", sep = "_"), SDM_GBM_200km)
  
  plot(SDM_GBM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith GBM 200km',sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GBM_200km.png",sep=""),sep = ""), plot = plot(SDM_GBM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith GBM 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  Bioclim_Stack_new <- stack(Bioclim_01+1.0, Bioclim_04, Bioclim_12-20, Bioclim_15)
  
  buffer_208km <- raster::buffer(Coords, width = 208000)
  
  Crop_208km_new <- crop(Bioclim_Stack_new, extent(buffer_208km))
  
  Boundary_208km_new <- stack(mask(Crop_208km_new, buffer_208km))
  
  SDM_projection_208km <- SSDM::project(SDM_GBM_200km, Boundary_208km_new)
  
  plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GBM 208km +1 temp and -20 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/GBM_208km_RCP2_6.png",sep=""),sep = ""), plot = plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GBM 208km +1 temp and -20 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Bioclim_Stack_new_2 <- stack(Bioclim_01+4.0, Bioclim_04, Bioclim_12-55, Bioclim_15)
  
  Crop_208km_new_2 <- crop(Bioclim_Stack_new_2, extent(buffer_208km))
  
  Boundary_208km_new_2 <- stack(mask(Crop_208km_new_2, buffer_208km))
  
  SDM_projection_2_208km <- SSDM::project(SDM_GBM_200km, Boundary_208km_new_2)
  
  plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GBM 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/GBM_208km_RCP8_5.png",sep=""),sep = ""), plot = plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith GBM 208km +4 temp and -55 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating evaluation table
  
  xy_presence <- SDM_GBM_200km@data[, 1:3]
  
  xy <- xy_presence[,c("X","Y")]
  
  sp_coord <- SpatialPointsDataFrame(coords = xy, data = xy_presence,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  existence_probability <- raster::extract(SDM_GBM_200km@projection, sp_coord, cellnumbers=TRUE)
  
  accuracy <- mean(abs(xy_presence[,3]-existence_probability[,2]))
  
  discrimination <- SDM_GBM_200km@evaluation[, "AUC"]
  
  calibration <- abs(sum(existence_probability[, "Projection"])/nrow(existence_probability)-sum(xy_presence[, "Presence"])/nrow(xy_presence))
  
  precision <- sum(sqrt(existence_probability[, "Projection"]*(1-existence_probability[, "Projection"])))/nrow(existence_probability)
  
  evaluate_table <- rbind(c(accuracy, discrimination, calibration, precision))
  
  colnames(evaluate_table) <- c("accuracy", "discrimination", "calibration", "precision")
  
  write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GBM_table",sep=""),sep = ""), evaluate_table)
  
  print(i)
  
}

rm(list = ls(pattern = "GBM*"))


##Using CTA method

for (i in 1:nrow(Taxon_lim)){
  
  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  
  ### Creating CTA models
  
  SDM_CTA_200km<- modelling('CTA', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "CTA_200km", sep = "_"), SDM_CTA_200km)
  
  plot(SDM_CTA_200km@projection, main = paste(Taxon_lim[i,1], '\nwith CTA 200km',sep = ""))
  
  ggsave(paste("../figs/",paste(Taxon_lim[i,1],"/CTA_200km.png",sep=""),sep = ""), plot = plot(SDM_CTA_200km@projection, main = paste(Taxon_lim[i,1], '\nwith CTA 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  buffer_205km <- raster::buffer(Coords, width = 205000)
  
  Crop_205km_new <- crop(Bioclim_Stack_RCP26, extent(buffer_205km))
  
  Boundary_205km_new <- stack(mask(Crop_205km_new, buffer_205km))
  
  SDM_projection_205km <- SSDM::project(SDM_CTA_200km, Boundary_205km_new)
  
  plot(SDM_projection_205km@projection, main = paste(Taxon_lim[i,1], '\nwith CTA RCP2.6', sep = ""))
  
  ggsave(paste("../figs/", paste(Taxon_lim[i,1],"/CTA_205km_RCP26.png",sep=""),sep = ""), plot = plot(SDM_projection_205km@projection, main = paste(Taxon_lim[i,1], '\nwith CTA 205km RCP2.6', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Crop_205km_new_2 <- crop(Bioclim_Stack_RCP85, extent(buffer_205km))
  
  Boundary_205km_new_2 <- stack(mask(Crop_205km_new_2, buffer_205km))
  
  SDM_projection_2_205km <- SSDM::project(SDM_CTA_200km, Boundary_205km_new_2)
  
  plot(SDM_projection_2_205km@projection, main = paste(Taxon_lim[i,1], '\nwith CTA 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("../figs/", paste(Taxon_lim[i,1],"/CTA_205km_RCP85.png",sep=""),sep = ""), plot = plot(SDM_projection_2_205km@projection, main = paste(Taxon_lim[i,1], '\nwith CTA 205km RCP8.5', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Compare Presence vs RCP2.6 and RCP8.5
  
  SDM_CTA_200km@projection <- readAll(SDM_CTA_200km@projection)
  
  SDM_projection_205km@projection <- readAll(SDM_projection_205km@projection)
  
  SDM_projection_2_205km@projection <- readAll(SDM_projection_2_205km@projection)
  
  Presence_value <- sum(SDM_CTA_200km@projection@data@values, na.rm = TRUE)/length(SDM_CTA_200km@projection@data@values)
  
  RCP26_value <- sum(SDM_projection_205km@projection@data@values, na.rm = TRUE)/length(SDM_projection_205km@projection@data@values)
  
  RCP85_value <- sum(SDM_projection_2_205km@projection@data@values, na.rm = TRUE)/length(SDM_projection_2_205km@projection@data@values)
  
  Scenarios <- cbind(log(RCP26_value/Presence_value), log(RCP85_value/Presence_value))
  
  colnames(Scenarios) <- c("log(RCP2.6/Presence)", "log(RCP8.5/Presence)")
  
  write.table(file = paste("../figs/", paste(Taxon_lim[i,1],"/Scenarios_CTA",sep=""),sep = ""), Scenarios)
  
  print(i)
}

rm(list = ls(pattern = "CTA*"))


##Using RF method

for (i in c(29:43, 45:nrow(Taxon_lim))){

  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  
  ### Creating RF models
  
  SDM_RF_200km<- modelling('RF', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "RF_200km", sep = "_"), SDM_RF_200km)
  
  plot(SDM_RF_200km@projection, main = paste(Taxon_lim[i,1], '\nwith RF 200km',sep = ""))
  
  ggsave(paste("../figs/",paste(Taxon_lim[i,1],"/RF_200km.png",sep=""),sep = ""), plot = plot(SDM_RF_200km@projection, main = paste(Taxon_lim[i,1], '\nwith RF 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  buffer_205km <- raster::buffer(Coords, width = 205000)
  
  Crop_205km_new <- crop(Bioclim_Stack_RCP26, extent(buffer_205km))
  
  Boundary_205km_new <- stack(mask(Crop_205km_new, buffer_205km))
  
  SDM_projection_205km <- SSDM::project(SDM_RF_200km, Boundary_205km_new)
  
  plot(SDM_projection_205km@projection, main = paste(Taxon_lim[i,1], '\nwith RF RCP2.6', sep = ""))
  
  ggsave(paste("../figs/", paste(Taxon_lim[i,1],"/RF_205km_RCP26.png",sep=""),sep = ""), plot = plot(SDM_projection_205km@projection, main = paste(Taxon_lim[i,1], '\nwith RF 205km RCP2.6', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Crop_205km_new_2 <- crop(Bioclim_Stack_RCP85, extent(buffer_205km))
  
  Boundary_205km_new_2 <- stack(mask(Crop_205km_new_2, buffer_205km))
  
  SDM_projection_2_205km <- SSDM::project(SDM_RF_200km, Boundary_205km_new_2)
  
  plot(SDM_projection_2_205km@projection, main = paste(Taxon_lim[i,1], '\nwith RF 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("../figs/", paste(Taxon_lim[i,1],"/RF_205km_RCP85.png",sep=""),sep = ""), plot = plot(SDM_projection_2_205km@projection, main = paste(Taxon_lim[i,1], '\nwith RF 205km RCP8.5', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Compare Presence vs RCP2.6 and RCP8.5
  
  SDM_RF_200km@projection <- readAll(SDM_RF_200km@projection)
  
  SDM_projection_205km@projection <- readAll(SDM_projection_205km@projection)
  
  SDM_projection_2_205km@projection <- readAll(SDM_projection_2_205km@projection)
  
  Presence_value <- sum(SDM_RF_200km@projection@data@values, na.rm = TRUE)/length(SDM_RF_200km@projection@data@values)
  
  RCP26_value <- sum(SDM_projection_205km@projection@data@values, na.rm = TRUE)/length(SDM_projection_205km@projection@data@values)
  
  RCP85_value <- sum(SDM_projection_2_205km@projection@data@values, na.rm = TRUE)/length(SDM_projection_2_205km@projection@data@values)
  
  Scenarios <- cbind(log(RCP26_value/Presence_value), log(RCP85_value/Presence_value))
  
  colnames(Scenarios) <- c("log(RCP2.6/Presence)", "log(RCP8.5/Presence)")
  
  write.table(file = paste("../figs/", paste(Taxon_lim[i,1],"/Scenarios_RF",sep=""),sep = ""), Scenarios)
  
  print(i)
}

rm(list = ls(pattern = "RF*"))

