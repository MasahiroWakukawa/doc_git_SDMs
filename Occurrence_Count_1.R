rm(list=ls())
require(dplyr)
require(tidyr)
require(rgbif)
require(taxize)


## Initial setting
#library(usethis)
#usethis::use_git()

## Create a user-defined function in R: 
#square_number <- function(base){
#  square <- base*base 
#  return(square)
#}

## Calling a user-defined function in R 
#square_number(5)

## Additional argument to R function: 
#exp_number <- function(base, power) {
#  exp <- base^power 
#  return(exp)
#}
#sessionInfo()

## Reading data test
#PREDICTS2 <- readRDS("../taxa-2020-02-20-03-32-15-rds/taxa-2020-02-20-03-32-15.rds") ### Masa data
#PREDICTS2 <- read.delim("../taxa-2020-02-20-03-32-15-rds/file_sample.csv", header = TRUE, stringsAsFactors = FALSE, 
#                        colClasses = c(rep("NULL", 7), rep("character", 3), rep("NULL", 23), rep("integer", 1), rep("NULL", 16)))
#sample <- read.delim("../taxa-2020-02-20-03-32-15-rds/file_sample.csv",header = TRUE)
#sample
#names(sample)


### Reading real data of bee species
PREDICTS_Andrenidae <- read.delim("../data/Andrenidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                  colClasses = c(rep("NULL", 7), rep("factor", 3), rep("NULL", 5), rep("factor", 1), rep("NULL", 5),
                                                 rep("factor", 2), rep("NULL", 10), rep("factor", 1), rep("NULL", 16)), quote="")

PREDICTS_Apidae <- read.delim("../data/Apidae.csv", header = TRUE, stringsAsFactors = FALSE,
                              colClasses = c(rep("NULL", 7), rep("factor", 3), rep("NULL", 5), rep("factor", 1), rep("NULL", 5),
                                             rep("factor", 2), rep("NULL", 10), rep("factor", 1), rep("NULL", 16)), quote="")

PREDICTS_Colletidae <- read.delim("../data/Colletidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                  colClasses = c(rep("NULL", 7), rep("factor", 3), rep("NULL", 5), rep("factor", 1), rep("NULL", 5),
                                                 rep("factor", 2), rep("NULL", 10), rep("factor", 1), rep("NULL", 16)), quote="")

PREDICTS_Halictidae <- read.delim("../data/Halictidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                  colClasses = c(rep("NULL", 7), rep("factor", 3), rep("NULL", 5), rep("factor", 1), rep("NULL", 5),
                                                 rep("factor", 2), rep("NULL", 10), rep("factor", 1), rep("NULL", 16)), quote="")

PREDICTS_Megachilidae <- read.delim("../data/Megachilidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                    colClasses = c(rep("NULL", 7), rep("factor", 3), rep("NULL", 5), rep("factor", 1), rep("NULL", 5),
                                                   rep("factor", 2), rep("NULL", 10), rep("factor", 1), rep("NULL", 16)), quote="")

PREDICTS_Melittidae <- read.delim("../data/Melittidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                  colClasses = c(rep("NULL", 7), rep("factor", 3), rep("NULL", 5), rep("factor", 1), rep("NULL", 5),
                                                 rep("factor", 2), rep("NULL", 10), rep("factor", 1), rep("NULL", 16)), quote="")

PREDICTS_Stenotritidae <- read.delim("../data/Stenotritidae.csv", header = TRUE, stringsAsFactors = FALSE,
                                     colClasses = c(rep("NULL", 7), rep("factor", 3), rep("NULL", 5), rep("factor", 1), rep("NULL", 5),
                                                    rep("factor", 2), rep("NULL", 10), rep("factor", 1), rep("NULL", 16)), quote="")

PREDICTS <- rbind(PREDICTS_Andrenidae, PREDICTS_Apidae, PREDICTS_Colletidae, PREDICTS_Halictidae,
                  PREDICTS_Megachilidae, PREDICTS_Melittidae, PREDICTS_Stenotritidae)

## Check basic functions
#names(PREDICTS) 
#nrow(PREDICTS)
#head(PREDICTS)
#class(PREDICTS)
#?scan
#?read.delim

##Patrick code
#PREDICTS <- readRDS("Raw_Data/Taxa/taxa-2019-05-09-02-34-47.rds")       #### load PREDICTS Taxa list
#PREDICTS$Scientific_Name <- paste(PREDICTS$genus,PREDICTS$species)  #### create column with scientific name
#PREDICTS<- PREDICTS %>%
#  drop_na(species)                                                      ### rm columns with NA
#PREDICTS <- PREDICTS %>% 
#  filter(family %in% Bee_Families)                ## Subset taxa list so that it only contains species in the bee families 
#PREDICTS_Species <- data.frame(unique(PREDICTS$Scientific_Name)) ### create data frame with just the species list from PREDICTS
#PREDICTS_Species
#PREDICTS2_Species=filter(PREDICTS2_Species, unique.PREDICTS.Scientific_Name.!=" ")
#head(PREDICTS_Species)


### Filtering the data
PREDICTS<- PREDICTS %>% filter(species!="" & decimalLatitude!="" & decimalLongitude!="")

### Set some initial definitions
Bee_Families <- c("Andrenidae","Apidae","Colletidae","Halictidae","Megachilidae","Melittidae","Stenotritidae")  ## vector with the names of the seven bee families

Bee_Species <- data.frame(unique(PREDICTS$species))

European_Countries<-c("AT","AD","AL","AX","BA","BY","BE","BG","CY","GG","GI","JE","HR","CH","CY","CZ","DK","EE","FI","FR","DE",
                      "GR","HU","HR","IE","IT","LI","LV","LT","LU","ME","MC","MT","MK","MT","MD","NL","NO","PL","PT","VA","RO",
                      "RU","SK","SI","SM","ES","SE","UA","GB")  ### Vector of European CountryCodes  

save(file = "../output/RData_All_Bee_Speices.RData", Bee_Species)

##Patrick code
## get the family IDs for the three databases we are going to query to get list of species names
#TaxonKey <- taxize::get_ids(names = Bee_Families, db = c("itis","gbif"))   ###the data of "col" is not accessible
#TaxonKey <- taxize::get_ids(names = Bee_Families, db = "gbif")
#?get_ids
#remotes::install_github("ropensci/colpluz")   ###install colpluz package
#library("colpluz")
#?colpluz
#TaxonKey[[1]]
#TaxonKey[[2]]
#TaxonKey_Bee <- data.frame(Bee_Families, TaxonKey[[1]], TaxonKey[[2]])
#TaxonKey_Bee <- data.frame(Bee_Families, TaxonKey[[1]])
#TaxonKey_Bee
#TaxonKey_Bee <- TaxonKey_Bee[,c(1:2,8)]
#TaxonKey_Bee <- TaxonKey_Bee[,c(1:2)]
#colnames(TaxonKey_Bee)[2:3] <- c("itis", "gbif")
#colnames(TaxonKey_Bee)[2] <- "gbif"
#TaxonKey_Bee

##Patrick code
## query each database in turn 
#for(i in 1:length(Bee_Families)){
#  itis_bee_down <- itis_downstream(tsns = TaxonKey_Bee[i,2], downto = "species", intermediate = FALSE)         ## Catelouge of Life
#  gbif_bee_down <- gbif_downstream(key = TaxonKey_Bee[i,3], downto = "species", intermediate = FALSE, start = 1, limit = 100000000) ## Gbif database
#  Bee_Species <- rbind(itis_bee_down[[1]][2],gbif_bee_down$name)
#  Bee_Species <- Bee_Species %>% distinct(childtaxa_name)
#  assign(paste(Bee_Families[i],"Species", sep = "_"),Bee_Species)
#}
#for(i in 1:length(Bee_Families)){        
#  gbif_bee_down <- gbif_downstream(key = TaxonKey_Bee[i,2], downto = "species", intermediate = FALSE, start = 1, limit = 100000000) ## Gbif database
#  Bee_Species <- data.frame(gbif_bee_down$name)
#  assign(paste(Bee_Families[i],"Species", sep = "_"),Bee_Species)
#}
#Bee_Species
#gbif_bee_down
#?distinct
#?itis_downstream
#?gbif_downstream
#colnames(PREDICTS_Species)="gbif_bee_down.name"  
#Bee_Species <- rbind(Andrenidae_Species,         ### combine all the species downstream of the family and the species in the PREDICTS Database
#                     Apidae_Species,
#                     Colletidae_Species,
#                     Halictidae_Species,
#                     Megachilidae_Species,
#                     Melittidae_Species,
#                     Stenotritidae_Species,
#                     PREDICTS_Species)
#Bee_Species <- Apidae_Species
#Bee_Species
#colnames(Bee_Species)[1] <- c("Scientific_Name")           ### Rename column
#Bee_Species <- Bee_Species %>%              ### Keep only the unique species -- PREDICTS Database adds a single species -- Halictus gemmeus
#  distinct(Scientific_Name)
#save(file = "RData_All_Bee_Speices.RData", Bee_Species)

##Patrick code
##generated the taxon key for each species in the list so that I can query GBIF
#UsageKey<-c()
#for(i in 1:nrow(Bee_Species)){
#  Taxon <- name_backbone(name=(Bee_Species[i,1]),rank="Species",kingdom="Animalia")   ###look up names in the GBIF backbone taxonomy
#  TaxonKey <- Taxon$usageKey
#  UsageKey[i]<-TaxonKey
#}
#?name_backbone
#TaxonKey
#?unique
#Bee_TaxonKey <- data.frame(Bee_Species$Scientific_Name, UsageKey)


###Creat Bee_Occurrence data
Bee_TaxonKey <- data.frame(PREDICTS$species, PREDICTS$taxonKey)

Bee_Occurrence_Count <- Bee_TaxonKey %>% group_by_all() %>% summarise(COUNT=n())

colnames(Bee_Occurrence_Count)[1:3]<-c("Scientific_Name","GBIF_TaxonKey","GBIF_Occurrence_Count")

Bee_Occurrence_Count <- Bee_Occurrence_Count %>% filter(GBIF_Occurrence_Count > 20)

##Patrick code
#Bee_TaxonKey <- distinct(Bee_TaxonKey) 
#Occurrence<-c()
#for(i in 1:nrow(Bee_TaxonKey)){
#  Occurrence[i]<-occ_count(taxonKey = (Bee_TaxonKey[i,2]), georeferenced = TRUE)
#}
#Bee_Occurrence_Count<-data.frame(Bee_TaxonKey,Occurrence)
#Bee_Occurrence_Count


##Patrick code
##Queried GBIF to find the occurrence data for each species
##Then going to query again but now only those species that occur in europe --- Here I think it becomes a bit patchy with GBIF as not all occurrences in the GBIF database
#have been classified a continent so will be missing a lot of data. So I have got the country codes for all european countries and then summed the occurrence data from 
#each country.Hopefully identifying those species which have sufficient occurrence data in Europe.

#Europe_Occ<-c()
#Country_Count<-c()
#Full_Count<-c()
#for(i in 1:nrow(Bee_Occurrence_Count)){
#  Europe=occ_search(taxonKey = Bee_Occurrence_Count[i,2], hasCoordinate = TRUE, hasGeospatialIssue = FALSE, country = c(European_Countries), return = "meta")
#  print(i)
#for(j in 1:50){
#    Europe_Continent=Europe[[j]][["count"]]
#    Country_Count[j]=sum(Europe_Continent)
#  }
#  Full_Count=sum(Country_Count)
#  Europe_Occ[i]<-sum(Full_Count)
#  print(i)
#}


### Creat Europe_Bee_Occurrence data
Bee_TaxonKey_Europe <- data.frame(PREDICTS$species, PREDICTS$taxonKey, PREDICTS$countryCode)

Bee_TaxonKey_Europe <- Bee_TaxonKey_Europe %>% filter(PREDICTS.countryCode %in% European_Countries)

Europe_Occ <- Bee_TaxonKey_Europe %>% group_by(PREDICTS.species, PREDICTS.taxonKey) %>% summarise(COUNT=n())

Europe_Occ <- Europe_Occ %>% filter(COUNT > 20)

colnames(Europe_Occ)[3]<-c("GBIF_Europe_Occurrence_Count")

Bee_Occurrence_Count <- Bee_Occurrence_Count %>% filter(GBIF_TaxonKey %in% Europe_Occ$PREDICTS.taxonKey)

Bee_Occurrence_Count <- as.data.frame(Bee_Occurrence_Count)

Europe_Occurrence_Count <- cbind(Bee_Occurrence_Count, Europe_Occ[3])

European_Bee_Species <- Europe_Occurrence_Count %>% arrange(Scientific_Name)

save(file ="../output/RData_European_Bee_Species.RData", European_Bee_Species)

##Patrick code
#Europe_Occurrence_Count<-cbind(Bee_Occurrence_Count,Europe_Occ)
#colnames(Europe_Occurrence_Count)[4]<-c("GBIF_Europe_Occurrence_Count")
#head(Europe_Occurrence_Count)

#European_Bee_Species <- Europe_Occurrence_Count %>%
#  filter(GBIF_Europe_Occurrence_Count > 20) %>%
#  arrange(Scientific_Name)

