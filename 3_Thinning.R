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

load("../output/RData_Plant_Species_Clean_Data.RData")                          


### Thinning by thin function with 5km distance
Taxon_lim <- Clean_Data %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

CC_Data <- c()

for(i in 1:48){
  
  Plant_data <- Clean_Data %>% filter(species == Taxon_lim[i, 1])
  
  Coord <- Plant_data[, c("species", "decimallongitude","decimallatitude")]
  
  Thin <- thin(Coord, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
               thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
               write.log.file = FALSE, verbose = FALSE
               )
  
  Thin <- as.data.frame(Thin)
  
  Thin_5km <- Plant_data %>% filter(decimallongitude %in% Thin$Longitude & decimallatitude %in% Thin$Latitude)
  
  CC_Data <- rbind(CC_Data, Thin_5km)
  
  print(i)
  
}

## To devide Hordeum murinum species dataset, which has too many occurrence records

Plant_data <- Clean_Data %>% filter(species == Taxon_lim[49, 1])

Plant_data_1 <- Plant_data[1:10000, ]

Plant_data_2 <- Plant_data[10001:20000, ]

Plant_data_3 <- Plant_data[20001:30000, ]

Plant_data_4 <- Plant_data[30001:40000, ]

Plant_data_5 <- Plant_data[40001:51754, ]

Coord1 <- Plant_data_1[, c("species", "decimallongitude","decimallatitude")]

Coord2 <- Plant_data_2[, c("species", "decimallongitude","decimallatitude")]

Coord3 <- Plant_data_3[, c("species", "decimallongitude","decimallatitude")]

Coord4 <- Plant_data_4[, c("species", "decimallongitude","decimallatitude")]

Coord5 <- Plant_data_5[, c("species", "decimallongitude","decimallatitude")]

Thin1 <- thin(Coord1, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
             thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
             write.log.file = FALSE, verbose = FALSE
)

Thin1 <- as.data.frame(Thin1)

Thin2 <- thin(Coord2, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
              thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
              write.log.file = FALSE, verbose = FALSE
)

Thin2 <- as.data.frame(Thin2)

Thin3 <- thin(Coord3, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
              thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
              write.log.file = FALSE, verbose = FALSE
)

Thin3 <- as.data.frame(Thin3)

Thin4 <- thin(Coord4, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
              thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
              write.log.file = FALSE, verbose = FALSE
)

Thin4 <- as.data.frame(Thin4)

Thin5 <- thin(Coord5, spec.col = "species", long.col = "decimallongitude", lat.col = "decimallatitude", 
              thin.par = 5, reps = 1, locs.thinned.list.return = TRUE, write.files = FALSE,
              write.log.file = FALSE, verbose = FALSE
)

Thin5 <- as.data.frame(Thin5)

Thin_5km_1 <- Plant_data %>% filter(decimallongitude %in% Thin1$Longitude & decimallatitude %in% Thin1$Latitude)

Thin_5km_2 <- Plant_data %>% filter(decimallongitude %in% Thin2$Longitude & decimallatitude %in% Thin2$Latitude)

Thin_5km_3 <- Plant_data %>% filter(decimallongitude %in% Thin3$Longitude & decimallatitude %in% Thin3$Latitude)

Thin_5km_4 <- Plant_data %>% filter(decimallongitude %in% Thin4$Longitude & decimallatitude %in% Thin4$Latitude)

Thin_5km_5 <- Plant_data %>% filter(decimallongitude %in% Thin5$Longitude & decimallatitude %in% Thin5$Latitude)

CC_Data <- rbind(CC_Data, Thin_5km_1, Thin_5km_2, Thin_5km_3, Thin_5km_4, Thin_5km_5)

Plant_TaxonKey <- data.frame(CC_Data$species)

###Save the thinned data with more than 20 records per species

Plant_Occurrence_Count <- Plant_TaxonKey %>% group_by(CC_Data.species) %>% summarise(COUNT=n())

Plant_Occurrence_Count <- as.data.frame(Plant_Occurrence_Count)

Plant_Occurrence_Count <- Plant_Occurrence_Count %>% filter(COUNT > 20)

Thinned_Data_5km <- CC_Data %>% filter(species %in% Plant_Occurrence_Count$CC_Data.species)

Taxon_lim <- Plant_Occurrence_Count

save(file = "../output/RData_Plant_Species_Thinned_Data_5km.RData", Thinned_Data_5km)

