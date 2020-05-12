rm(list = ls())
require(dplyr)
require(tidyr)
require(ape)
require(raster)
require(sp)
require(rgdal)
require(spThin)   ####jeffreyhanson branch for githbud download 
require(tiff)
require(doParallel) 
#require(gurobi)
#require(devtools)
#install_github("jeffreyhanson/spThin")

load("../output/RData_European_Bee_Species_Finish_Data.RData")                          ## load finish Data for the thinning process
#load("RData_EUropean_Bee_Species_GBIF_Data.RData")                            ## load GBIF_Data for maps to visualise effects of thinning



###Thinning by thin function with 5km distance
Taxon_lim <- Finish_Data %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

CC_Data <- c()

for(i in 1:65){
  
  GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  Coord <- GBIF[, c("species", "decimallongitude","decimallatitude")]
  
  Thin <- thin(Coord, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )
  
  Thin <- as.data.frame(Thin)
  
  Thin_5km <- GBIF %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)
  
  CC_Data <- rbind(CC_Data, Thin_5km)
  
  print(i)
  
}

## Apis mellifera species (Taxon_lim[66,]) have too many data, so thinning by separated small dataset
GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[66, 2])

GBIF_1 <- GBIF[1:4000,]

Coord_1 <- GBIF_1[, c("species", "decimallongitude","decimallatitude")]

Thin_1 <- thin(Coord_1, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_1 <- as.data.frame(Thin_1)

GBIF_2 <- GBIF[4001:8000,]

Coord_2 <- GBIF_2[, c("species", "decimallongitude","decimallatitude")]

Thin_2 <- thin(Coord_2, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_2 <- as.data.frame(Thin_2)

GBIF_3 <- GBIF[8001:12000,]

Coord_3 <- GBIF_3[, c("species", "decimallongitude","decimallatitude")]

Thin_3 <- thin(Coord_3, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_3 <- as.data.frame(Thin_3)

GBIF_4 <- GBIF[12001:16339,]

Coord_4 <- GBIF_4[, c("species", "decimallongitude","decimallatitude")]

Thin_4 <- thin(Coord_4, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_4 <- as.data.frame(Thin_4)

Thin <- rbind(Thin_1, Thin_2, Thin_3, Thin_4)

Thin_5km <- GBIF %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)

Coord_66 <- Thin_5km[, c("species", "decimallongitude","decimallatitude")]

Thin_66 <- thin(Coord_66, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
                thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
                write.log.file = FALSE, verbose = FALSE
                )

Thin_66 <- as.data.frame(Thin_66)

Thin_5km <- GBIF %>% filter(decimallongitude %in% Thin_66$Longitude & decimallatitude %in% Thin_66$Latitude)

CC_Data <- rbind(CC_Data, Thin_5km)

## Thinning by thin function again with 5km distance
for(i in 67:nrow(Taxon_lim)){
  
  GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  Coord <- GBIF[, c("species", "decimallongitude","decimallatitude")]
  
  Thin <- thin(Coord, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )
  
  Thin <- as.data.frame(Thin)
  
  Thin_5km <- GBIF %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)
  
  CC_Data <- rbind(CC_Data, Thin_5km)
  
  print(i)
  
}

Bee_TaxonKey <- data.frame(CC_Data$species, CC_Data$taxonKey)

Bee_Occurrence_Count <- Bee_TaxonKey %>% group_by_all() %>% summarise(COUNT=n())

Bee_Occurrence_Count <- Bee_Occurrence_Count %>% filter(COUNT > 20)

Thinned_Data_5km <- CC_Data %>% filter(taxonKey %in% Bee_Occurrence_Count$CC_Data.taxonKey)

save(file = "../output/RData_European_Bee_Species_Thinned_Data_5km.RData", Thinned_Data_5km)

?thin

### Thinning by thin function with 10km distance
Taxon_lim <- Finish_Data %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

CC_Data <- c()

for(i in 1:65){
  
  GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  Coord <- GBIF[, c("species", "decimallongitude","decimallatitude")]
  
  Thin <- thin(Coord, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )
  
  Thin <- as.data.frame(Thin)
  
  Thin_10km <- GBIF %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)
  
  CC_Data <- rbind(CC_Data, Thin_10km)
  
  print(i)
  
}

## Apis mellifera species (Taxon_lim[66,]) have too many data, so thinning by separated small dataset
GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[66, 2])

GBIF_1 <- GBIF[1:4000,]

Coord_1 <- GBIF_1[, c("species", "decimallongitude","decimallatitude")]

Thin_1 <- thin(Coord_1, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_1 <- as.data.frame(Thin_1)

GBIF_2 <- GBIF[4001:8000,]

Coord_2 <- GBIF_2[, c("species", "decimallongitude","decimallatitude")]

Thin_2 <- thin(Coord_2, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_2 <- as.data.frame(Thin_2)

GBIF_3 <- GBIF[8001:12000,]

Coord_3 <- GBIF_3[, c("species", "decimallongitude","decimallatitude")]

Thin_3 <- thin(Coord_3, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_3 <- as.data.frame(Thin_3)

GBIF_4 <- GBIF[12001:16339,]

Coord_4 <- GBIF_4[, c("species", "decimallongitude","decimallatitude")]

Thin_4 <- thin(Coord_4, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )

Thin_4 <- as.data.frame(Thin_4)

Thin <- rbind(Thin_1, Thin_2, Thin_3, Thin_4)

Thin_10km <- GBIF %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)

Coord_66 <- Thin_10km[, c("species", "decimallongitude","decimallatitude")]

Thin_66 <- thin(Coord_66, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
                thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
                write.log.file = FALSE, verbose = FALSE
                )

Thin_66 <- as.data.frame(Thin_66)

Thin_10km <- GBIF %>% filter(decimallongitude %in% Thin_66$Longitude & decimallatitude %in% Thin_66$Latitude)

CC_Data <- rbind(CC_Data, Thin_10km)

## Thinning by thin function again with 10km distance
for(i in 67:nrow(Taxon_lim)){
  
  GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  Coord <- GBIF[, c("species", "decimallongitude","decimallatitude")]
  
  Thin <- thin(Coord, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 10, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )
  
  Thin <- as.data.frame(Thin)
  
  Thin_10km <- GBIF %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)
  
  CC_Data <- rbind(CC_Data, Thin_10km)
  
  print(i)
  
}

Bee_TaxonKey <- data.frame(CC_Data$species, CC_Data$taxonKey)

Bee_Occurrence_Count <- Bee_TaxonKey %>% group_by_all() %>% summarise(COUNT=n())

Bee_Occurrence_Count <- Bee_Occurrence_Count %>% filter(COUNT > 20)

Thinned_Data_10km <- CC_Data %>% filter(taxonKey %in% Bee_Occurrence_Count$CC_Data.taxonKey)

save(file = "../output/RData_European_Bee_Species_Thinned_Data_10km.RData", Thinned_Data_10km)



##Patrick code
#registerDoParallel(cores= 4)

#Thinned_Data_5km <- foreach(i = 1:length(Finish_Data),
#                            .packages = c("dplyr","tidyr","ape","raster","sp","rgdal","spThin","tiff","gurobi")) %dopar% {
#                              Coordinates <- Finish_Data[[i]][,c("decimallongitude","decimallatitude")]   ### extract the coordinate columns from Finish Data
#                              Thin <- spThin(Coordinates,                                 ## "spThin" package will perform the thinning, this will randomly remove
#                                             x.col = "decimallongitude",                  ## points that are within a certain distance from each other -- the algorithm 
#                                             y.col = "decimallatitude",                   ## "gurobi will then run to identify the maximum number of points that are at least
#                                             dist = 5000,                                 ## that distance apart to leave an optimum dataset. -- This takes a while to run 
#                                             method = "gurobi",                           ## so be wary.
#                                             great.circle.distance = TRUE
#                              )
#                              t1 <- Thin[[1]]                                    ## Extract the dataframe for the output
#                              Coords_Thin <- data.frame(t1@coords)              ## Extract Coordinates 
#                              Coords_Thin$min_Dist <- Thin@mindist                                                          ## Extract minimum distance points are from each other 
#                              assign(paste(Finish_Data[[i]][["species"]][1],"Thinned_Data_5km",sep = "_"),Coords_Thin)   #### assign dataframe to environment
#                            }

#registerDoSEQ()

#for(i in 1:length(Thinned_Data_5km)){
#  names(Thinned_Data_5km)[i] <- paste(Finish_Data[[i]]["species"][1],"Thinned_Data_5km",sep = "_")
#}

#for(i in 1:length(Thinned_Data_5km)){                                           ## Change coordinate column names to match previous datasets
#  colnames(Thinned_Data_5km[[i]])[1:2]<-c("decimallongitude","decimallatitude")
#  Thinned_Data_5km[[i]]$species <- Finish_Data[[i]][["species"]][1]
#}

#for(i in 1:length(Thinned_Data_5km)){                   #### Remove those cases that have less that 25 points 
#  if(nrow(Thinned_Data_5km[[i]]) > 25){
#    assign(paste(Thinned_Data_5km[[i]][["species"]][1],"Thinned_Data_5km",sep = "_"),Thinned_Data_5km[[i]])
#  }
#}

#rm(Thinned_Data_5km)
#Thinned_Data_5km<-setNames(lapply(ls(pattern ="Thinned_Data_5km*"), function(x) get(x)),(ls(pattern="Thinned_Data_5km*")))




#save(file = "RData_European_Bee_Species_Thinned_Data_5km.RData",Thinned_Data_5km)     ### Save 


#for(i in 1:length(Finish_Data)){                                               ### Repeat for 10km minimum distance apart -- This took 6 hours. 
#  Coordinates <- Finish_Data[[i]][,c("decimallongitude","decimallatitude")]
#  Thin <- spThin(Coordinates,
#                 x.col = "decimallongitude",
#                 y.col = "decimallatitude",
#                 dist = 10000,
#                 method = "gurobi",
#                 great.circle.distance = TRUE
#  )
#  t1 <- Thin[[1]]
#  Coords_Thin <- data.frame(t1@coords)
#  Coords_Thin$min_Dist <- Thin@mindist
#  assign(paste(Finish_Data[[i]][["species"]][1],"Thinned_Data_10km",sep = "_"),Coords_Thin)
#}

#Thinned_Data_10km<-setNames(lapply(ls(pattern ="Thinned_Data_10km*"), function(x) get(x)),(ls(pattern="Thinned_Data_10km*")))


#for(i in 1:length(Thinned_Data_10km)){
#  colnames(Thinned_Data_10km[[i]])[1:2]<-c("decimallongitude","decimallatitude")
#  Thinned_Data_10km[[i]]$species <- Finish_Data[[i]][["species"]][1]
#}

#for(i in 1:length(Thinned_Data_10km)){                   #### Remove those cases that have less that 25 points 
#  if(nrow(Thinned_Data_10km[[i]]) > 25){
#    assign(paste(Thinned_Data_10km[[i]][["species"]][1],"Thinned_Data_10km",sep = "_"),Thinned_Data_10km[[i]])
#  }
#}
#Thinned_Data_10km<-setNames(lapply(ls(pattern ="Thinned_Data_10km*"), function(x) get(x)),(ls(pattern="Thinned_Data_10km*")))


#save(file = "RData_European_Bee_Species_Thinned_Data_10km.RData",Thinned_Data_10km)

##Masa test code
#Taxon_lim <- Finish_Data %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

#Taxon_lim <- as.data.frame(Taxon_lim)

#GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[66, 2])

#Coord <- GBIF[, c("species", "decimallongitude","decimallatitude")]

#Thin <- thin(Coord, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
#             thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
#             write.log.file = FALSE, verbose = FALSE
#             )

#Thin <- as.data.frame(Thin)

#Thin_5km <- GBIF %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)

