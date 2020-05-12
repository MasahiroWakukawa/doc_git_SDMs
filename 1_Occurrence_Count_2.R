rm(list=ls())
require(dplyr)
require(tidyr)
require(rgbif)
require(magrittr)


### Reading real data of bee species
PREDICTS_Andrenidae <- read.delim("../data/Andrenidae.csv", header = TRUE, stringsAsFactors = FALSE, 
                                  colClasses = c(rep("factor", 1), rep("NULL", 6), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 3),
                                                 rep("factor", 1), rep("NULL", 3), rep("factor", 1), rep("NULL", 1), rep("factor", 3), rep("NULL", 8), rep("factor", 2), 
                                                 rep("NULL", 1), rep("factor", 2), rep("NULL", 13)), quote="")

PREDICTS_Apidae <- read.delim("../data/Apidae.csv", header = TRUE, stringsAsFactors = FALSE,
                              colClasses = c(rep("factor", 1), rep("NULL", 6), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 3),
                                             rep("factor", 1), rep("NULL", 3), rep("factor", 1), rep("NULL", 1), rep("factor", 3), rep("NULL", 8), rep("factor", 2), 
                                             rep("NULL", 1), rep("factor", 2), rep("NULL", 13)), quote="")

PREDICTS_Colletidae <- read.delim("../data/Colletidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                  colClasses = c(rep("factor", 1), rep("NULL", 6), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 3),
                                                 rep("factor", 1), rep("NULL", 3), rep("factor", 1), rep("NULL", 1), rep("factor", 3), rep("NULL", 8), rep("factor", 2), 
                                                 rep("NULL", 1), rep("factor", 2), rep("NULL", 13)), quote="")

PREDICTS_Halictidae <- read.delim("../data/Halictidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                  colClasses = c(rep("factor", 1), rep("NULL", 6), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 3),
                                                 rep("factor", 1), rep("NULL", 3), rep("factor", 1), rep("NULL", 1), rep("factor", 3), rep("NULL", 8), rep("factor", 2), 
                                                 rep("NULL", 1), rep("factor", 2), rep("NULL", 13)), quote="")

PREDICTS_Megachilidae <- read.delim("../data/Megachilidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                    colClasses = c(rep("factor", 1), rep("NULL", 6), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 3),
                                                   rep("factor", 1), rep("NULL", 3), rep("factor", 1), rep("NULL", 1), rep("factor", 3), rep("NULL", 8), rep("factor", 2), 
                                                   rep("NULL", 1), rep("factor", 2), rep("NULL", 13)), quote="")

PREDICTS_Melittidae <- read.delim("../data/Melittidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                  colClasses = c(rep("factor", 1), rep("NULL", 6), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 3),
                                                 rep("factor", 1), rep("NULL", 3), rep("factor", 1), rep("NULL", 1), rep("factor", 3), rep("NULL", 8), rep("factor", 2), 
                                                 rep("NULL", 1), rep("factor", 2), rep("NULL", 13)), quote="")

PREDICTS_Stenotritidae <- read.delim("../data/Stenotritidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                     colClasses = c(rep("factor", 1), rep("NULL", 6), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 1), rep("factor", 1), rep("NULL", 3),
                                                    rep("factor", 1), rep("NULL", 3), rep("factor", 1), rep("NULL", 1), rep("factor", 3), rep("NULL", 8), rep("factor", 2), 
                                                    rep("NULL", 1), rep("factor", 2), rep("NULL", 13)), quote="")

PREDICTS <- rbind(PREDICTS_Andrenidae, PREDICTS_Apidae, PREDICTS_Colletidae, PREDICTS_Halictidae,
                  PREDICTS_Megachilidae, PREDICTS_Melittidae, PREDICTS_Stenotritidae)

load("../output/RData_European_Bee_Species.RData")

Species<-European_Bee_Species

### Cleaning the data by name, taxonkey, latitude and longitude
PREDICTS <- PREDICTS %>% filter(species %in% Species$Scientific_Name)

PREDICTS <- PREDICTS %>% filter(taxonKey %in% Species$GBIF_TaxonKey)

PREDICTS <- PREDICTS %>% filter(species!="" & decimalLatitude!="" & decimalLongitude!="")

colnames(PREDICTS)[7:8]<-c("decimallongitude","decimallatitude")

GBIF_Data <- PREDICTS

save(file = "../output/RData_European_Bee_Species_GBIF_Data.RData", GBIF_Data)


##Patrick code
## Extracting the decimal coordinates of the occurrences of each european Species from GBIF.
## have also extended the code to extract extra information from the meta-data that may infrom the data cleaning process
#for(i in 1:nrow(Species)){
#  Data <- occ_search(taxonKey = Species[i,2], hasCoordinate = TRUE,hasGeospatialIssue = FALSE,limit = 150000, return = "data")
#  GBIF <- Data %>%
#    dplyr::select(one_of("species", "decimalLongitude", "decimalLatitude", "countryCode", "individualCount",
#                         "gbifID", "family", "taxonRank", "coordinateUncertaintyInMeters", "year",
#                         "basisOfRecord", "institutionCode", "datasetName", "geodeticDatum",
#                         "dataGeneralizations","verbatimSRS","verbatimCoordinateSystem","georeferenceRemarks"))
#  colnames(GBIF)[2:3]<-c("decimallongitude","decimallatitude")
#  assign(paste(Species[i,1],"GBIF_Data",sep="_"),GBIF)
#}
#nrow(Species)
#for(i in 1:nrow(Species)){
#  Data <- occ_search(taxonKey = Species[i,2], hasCoordinate = TRUE,hasGeospatialIssue = FALSE,limit = 150000, return = "data")
#  GBIF <- Data %>%
#    dplyr::select(one_of("species", "decimalLongitude", "decimalLatitude", "countryCode", "individualCount",
#                         "gbifID", "family", "taxonRank", "coordinateUncertaintyInMeters", "year",
#                         "basisOfRecord", "institutionCode", "datasetName", "geodeticDatum",
#                         "dataGeneralizations","verbatimSRS","verbatimCoordinateSystem","georeferenceRemarks"))
#  colnames(GBIF)[2:3]<-c("decimallongitude","decimallatitude")
#  assign(paste(Species[i,1],"GBIF_Data",sep="_"),GBIF)
#  print(i)
#}

#PREDICTS <- PREDICTS %>% filter(species %in% Species$Scientific_Name)
#PREDICTS <- PREDICTS %>% filter(taxonKey %in% Species$GBIF_TaxonKey)
#PREDICTS <- PREDICTS %>% filter(species!="" & decimalLatitude!="" & decimalLongitude!="")

### In the future, this part is necessary to generate each species' data sets###
#for(i in 1:nrow(Species)){
#  GBIF <- PREDICTS %>% filter(taxonKey == Species[i, 2])
#  colnames(GBIF)[7:8]<-c("decimallongitude","decimallatitude")
#  assign(paste(Species[i,1],"GBIF_Data",sep="_"),GBIF)
#  print(i)
#}
#GBIF_Data<-setNames(lapply(ls(pattern ="GBIF_Data*"), function(x) get(x)),(ls(pattern="GBIF_Data*")))

#save(file = "RData_European_Bee_Species_GBIF_Data.RData", GBIF_Data)

#GBIF_Data <- as.data.frame(GBIF_Data)
#GBIF_Data$Amegilla.albigena_GBIF_Data.datasetName

#GBIF_Species <- c()                                             #### Identify mis-matches between the scientific name and name given by GBIF
#for(i in 1:length(GBIF_Data)){                                   ### Could be indicative of species synonyms or outdated species names 
#GBIF_Species[i] <- GBIF_Data[[i]][["species"]][1]
#}
#?read.delim
#Species_1 <- cbind(Species,GBIF_Species)
#Species_1
#Species
#Species_1$Scientific_Name <- as.character(Species$Scientific_Name)
#Species_1$GBIF_Species <- as.character(Species$GBIF_Species)
#Species_1$GBIF_Species <- as.character(Species_1$GBIF_Species)

#rm_Species <- Species_1 %>%
#  filter(Scientific_Name != GBIF_Species)

#keep_Species <- Species_1 %>%
#  filter(Scientific_Name == GBIF_Species)

#### reviewing the species that have been removed due to mismatches only one -- sephcodes gibbus can be kept - name is a synonym but is still a distinct bee species that can be included.

#GBIF_Data[["Sphecodes gibbus_GBIF_Data"]][["species"]] = "Sphecodes gibbus"

#GBIF_Species_2 <- c()                                             #### Identify mis-matches between the scientific name and name given by GBIF
#for(i in 1:length(GBIF_Data)){                                   ### Could be indicative of species synonyms or outdated species names 
#  GBIF_Species_2[i] <- GBIF_Data[[i]][["species"]][1]
#}

#Species_2 <- cbind(Species,GBIF_Species_2)
#Species_2
#Species
#Species$Scientific_Name <- as.character(Species_2$Scientific_Name)
#Species_2$Scientific_Name <- as.character(Species_2$Scientific_Name)
#Species$GBIF_Species <- as.character(Species_2$GBIF_Species)
#Species_2$GBIF_Species <- as.character(Species_2$GBIF_Species)
#flags <- Species_2$Scientific_Name == Species_2$GBIF_Species_2
#flags
#GBIF_Data <- GBIF_Data[flags]
#GBIF_Data
#save(file = "RData_European_Bee_Species_GBIF_Data.RData", GBIF_Data)
