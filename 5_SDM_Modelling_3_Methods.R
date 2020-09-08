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
library(rJava)


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

Occurrence_Data_5km <- Thinned_Data_5km[c("species", "decimallongitude", "decimallatitude")]

Occurrence_Data_5km$species <- as.character(Occurrence_Data_5km$species)

Taxon_lim


##Using MAXENT method

for (i in 20:nrow(Taxon_lim)){
  
  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  
  ### Creating MAXENT models
  
  SDM_MAXENT_200km<- modelling('MAXENT', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "MAXENT_200km", sep = "_"), SDM_MAXENT_200km)
  
  plot(SDM_MAXENT_200km@projection, main = paste(Taxon_lim[i,1], '\nwith MAXENT 200km',sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/MAXENT_200km.png",sep=""),sep = ""), plot = plot(SDM_MAXENT_200km@projection, main = paste(Taxon_lim[i,1], '\nwith MAXENT 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  Bioclim_Stack_new <- stack(Bioclim_01+1.0, Bioclim_04, Bioclim_12-20, Bioclim_15)
  
  buffer_208km <- raster::buffer(Coords, width = 208000)
  
  Crop_208km_new <- crop(Bioclim_Stack_new, extent(buffer_208km))
  
  Boundary_208km_new <- stack(mask(Crop_208km_new, buffer_208km))
  
  SDM_projection_208km <- SSDM::project(SDM_MAXENT_200km, Boundary_208km_new)
  
  plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith MAXENT 208km +1 temp and -20 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/MAXENT_208km_RCP2_6.png",sep=""),sep = ""), plot = plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith MAXENT 208km +1 temp and -20 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Bioclim_Stack_new_2 <- stack(Bioclim_01+4.0, Bioclim_04, Bioclim_12-55, Bioclim_15)
  
  Crop_208km_new_2 <- crop(Bioclim_Stack_new_2, extent(buffer_208km))
  
  Boundary_208km_new_2 <- stack(mask(Crop_208km_new_2, buffer_208km))
  
  SDM_projection_2_208km <- SSDM::project(SDM_MAXENT_200km, Boundary_208km_new_2)
  
  plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith MAXENT 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/MAXENT_208km_RCP8_5.png",sep=""),sep = ""), plot = plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith MAXENT 208km +4 temp and -55 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating evaluation table
  
  xy_presence <- SDM_MAXENT_200km@data[, 1:3]
  
  xy <- xy_presence[,c("X","Y")]
  
  sp_coord <- SpatialPointsDataFrame(coords = xy, data = xy_presence,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  existence_probability <- raster::extract(SDM_MAXENT_200km@projection, sp_coord, cellnumbers=TRUE)
  
  accuracy <- mean(abs(xy_presence[,3]-existence_probability[,2]))
  
  discrimination <- SDM_MAXENT_200km@evaluation[, "AUC"]
  
  calibration <- abs(sum(existence_probability[, "Projection"])/nrow(existence_probability)-sum(xy_presence[, "Presence"])/nrow(xy_presence))
  
  precision <- sum(sqrt(existence_probability[, "Projection"]*(1-existence_probability[, "Projection"])))/nrow(existence_probability)
  
  evaluate_table <- rbind(c(accuracy, discrimination, calibration, precision))
  
  colnames(evaluate_table) <- c("accuracy", "discrimination", "calibration", "precision")
  
  write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/MAXENT_table",sep=""),sep = ""), evaluate_table)
  
  print(i)
  
}

rm(list = ls(pattern = "MAXENT*"))


##Using ANN method

for (i in 1:nrow(Taxon_lim)){
  
  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  
  ### Creating ANN models
  
  SDM_ANN_200km<- modelling('ANN', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "ANN_200km", sep = "_"), SDM_ANN_200km)
  
  plot(SDM_ANN_200km@projection, main = paste(Taxon_lim[i,1], '\nwith ANN 200km',sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/ANN_200km.png",sep=""),sep = ""), plot = plot(SDM_ANN_200km@projection, main = paste(Taxon_lim[i,1], '\nwith ANN 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  Bioclim_Stack_new <- stack(Bioclim_01+1.0, Bioclim_04, Bioclim_12-20, Bioclim_15)
  
  buffer_208km <- raster::buffer(Coords, width = 208000)
  
  Crop_208km_new <- crop(Bioclim_Stack_new, extent(buffer_208km))
  
  Boundary_208km_new <- stack(mask(Crop_208km_new, buffer_208km))
  
  SDM_projection_208km <- SSDM::project(SDM_ANN_200km, Boundary_208km_new)
  
  plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith ANN 208km +1 temp and -20 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/ANN_208km_RCP2_6.png",sep=""),sep = ""), plot = plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith ANN 208km +1 temp and -20 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Bioclim_Stack_new_2 <- stack(Bioclim_01+4.0, Bioclim_04, Bioclim_12-55, Bioclim_15)
  
  Crop_208km_new_2 <- crop(Bioclim_Stack_new_2, extent(buffer_208km))
  
  Boundary_208km_new_2 <- stack(mask(Crop_208km_new_2, buffer_208km))
  
  SDM_projection_2_208km <- SSDM::project(SDM_ANN_200km, Boundary_208km_new_2)
  
  plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith ANN 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/ANN_208km_RCP8_5.png",sep=""),sep = ""), plot = plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith ANN 208km +4 temp and -55 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating evaluation table
  
  xy_presence <- SDM_ANN_200km@data[, 1:3]
  
  xy <- xy_presence[,c("X","Y")]
  
  sp_coord <- SpatialPointsDataFrame(coords = xy, data = xy_presence,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  existence_probability <- raster::extract(SDM_ANN_200km@projection, sp_coord, cellnumbers=TRUE)
  
  accuracy <- mean(abs(xy_presence[,3]-existence_probability[,2]))
  
  discrimination <- SDM_ANN_200km@evaluation[, "AUC"]
  
  calibration <- abs(sum(existence_probability[, "Projection"])/nrow(existence_probability)-sum(xy_presence[, "Presence"])/nrow(xy_presence))
  
  precision <- sum(sqrt(existence_probability[, "Projection"]*(1-existence_probability[, "Projection"])))/nrow(existence_probability)
  
  evaluate_table <- rbind(c(accuracy, discrimination, calibration, precision))
  
  colnames(evaluate_table) <- c("accuracy", "discrimination", "calibration", "precision")
  
  write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/ANN_table",sep=""),sep = ""), evaluate_table)
  
  print(i)
  
}

rm(list = ls(pattern = "ANN*"))


##Using SVM method

for (i in 1:nrow(Taxon_lim)){
  
  OC <- subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1])
  
  Coords <- OC[,c("decimallongitude","decimallatitude")]
  
  Coords <- SpatialPointsDataFrame(coords = Coords, data = OC,
                                   proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  ### Creating buffer
  
  buffer_200km <- raster::buffer(Coords, width = 200000)
  
  Crop_200km <- crop(Bioclim_Stack, extent(buffer_200km))
  
  Boundary_200km <- stack(mask(Crop_200km, buffer_200km))
  
  
  ### Creating SVM models
  
  SDM_SVM_200km<- modelling('SVM', subset(Occurrence_Data_5km,  Occurrence_Data_5km$species == Taxon_lim[i, 1]), 
                            Boundary_200km, Xcol = 'decimallongitude', Ycol = 'decimallatitude', verbose = FALSE)
  
  assign(paste(Taxon_lim[i,1], "SVM_200km", sep = "_"), SDM_SVM_200km)
  
  plot(SDM_SVM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith SVM 200km',sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/SVM_200km.png",sep=""),sep = ""), plot = plot(SDM_SVM_200km@projection, main = paste(Taxon_lim[i,1], '\nwith SVM 200km',sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP2.6
  
  Bioclim_Stack_new <- stack(Bioclim_01+1.0, Bioclim_04, Bioclim_12-20, Bioclim_15)
  
  buffer_208km <- raster::buffer(Coords, width = 208000)
  
  Crop_208km_new <- crop(Bioclim_Stack_new, extent(buffer_208km))
  
  Boundary_208km_new <- stack(mask(Crop_208km_new, buffer_208km))
  
  SDM_projection_208km <- SSDM::project(SDM_SVM_200km, Boundary_208km_new)
  
  plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith SVM 208km +1 temp and -20 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/SVM_208km_RCP2_6.png",sep=""),sep = ""), plot = plot(SDM_projection_208km@projection, main = paste(Taxon_lim[i,1], '\nwith SVM 208km +1 temp and -20 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating new climatic environment for future prediction RCP8.5
  
  Bioclim_Stack_new_2 <- stack(Bioclim_01+4.0, Bioclim_04, Bioclim_12-55, Bioclim_15)
  
  Crop_208km_new_2 <- crop(Bioclim_Stack_new_2, extent(buffer_208km))
  
  Boundary_208km_new_2 <- stack(mask(Crop_208km_new_2, buffer_208km))
  
  SDM_projection_2_208km <- SSDM::project(SDM_SVM_200km, Boundary_208km_new_2)
  
  plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith SVM 208km +4 temp and -55 prec', sep = ""))
  
  ggsave(paste("~/Documents/Patrick_Masa_code/SDMs/figs/", paste(Taxon_lim[i,1],"/SVM_208km_RCP8_5.png",sep=""),sep = ""), plot = plot(SDM_projection_2_208km@projection, main = paste(Taxon_lim[i,1], '\nwith SVM 208km +4 temp and -55 prec', sep = "")), device = "png", height = 10, width = 10,dpi=300)
  
  ### Creating evaluation table
  
  xy_presence <- SDM_SVM_200km@data[, 1:3]
  
  xy <- xy_presence[,c("X","Y")]
  
  sp_coord <- SpatialPointsDataFrame(coords = xy, data = xy_presence,
                                     proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  existence_probability <- raster::extract(SDM_SVM_200km@projection, sp_coord, cellnumbers=TRUE)
  
  accuracy <- mean(abs(xy_presence[,3]-existence_probability[,2]))
  
  discrimination <- SDM_SVM_200km@evaluation[, "AUC"]
  
  calibration <- abs(sum(existence_probability[, "Projection"])/nrow(existence_probability)-sum(xy_presence[, "Presence"])/nrow(xy_presence))
  
  precision <- sum(sqrt(existence_probability[, "Projection"]*(1-existence_probability[, "Projection"])))/nrow(existence_probability)
  
  evaluate_table <- rbind(c(accuracy, discrimination, calibration, precision))
  
  colnames(evaluate_table) <- c("accuracy", "discrimination", "calibration", "precision")
  
  write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/SVM_table",sep=""),sep = ""), evaluate_table)
  
  print(i)
  
}

rm(list = ls(pattern = "SVM*"))


