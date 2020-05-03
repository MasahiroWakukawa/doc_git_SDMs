rm(list=ls())
require(dplyr)
require(tidyr)
require(rgbif)
require(magrittr)
require(countrycode)
require(CoordinateCleaner)
require(rnaturalearthdata)

load("../output/RData_European_Bee_Species_GBIF_Data.RData")

## Check basic functions
#class(GBIF_Data)
#length(GBIF_Data)
#?countrycode
#?codelist
#?clean_coordinates
#?setNames
#head(GBIF_Data)

## Patrick code
## Initally flag those coordinates that are a bit dubious.
#for(i in 1:length(GBIF_Data)){
#  dat<-GBIF_Data[[i]]%>%
#    filter(!is.na(decimallongitude)) %>%
#   filter(!is.na(decimallatitude))
#  dat$countryCode <-  countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')
#  dat_NO_NOR <- dat %>%filter(countryCode != "NOR")
#  dat_NOR <- dat %>% filter(countryCode == "NOR")
#  flags_NO_NOR <- clean_coordinates(x = dat_NO_NOR, lon = "decimallongitude", lat = "decimallatitude",
#                             countries = "countryCode", 
#                             species = "species",
#                             tests = c("capitals", "centroids", "equal","gbif", "institutions",
#                                       "zeros", "seas", "duplicates"),   ### remove "countries"
#                             seas_scale = 110)
#  if(NROW(dat_NOR)>0){
#  flags_NOR <- clean_coordinates(x = dat_NOR, lon = "decimallongitude", lat = "decimallatitude",
#                                    countries = "countryCode", 
#                                    species = "species",
#                                    tests = c("capitals", "centroids", "equal","gbif", "institutions",
#                                              "zeros","seas", "duplicates"),
#                                    seas_scale = 110)
#  flags_NOR$.con<- NA
#  } else {
#                                      flags_NOR <-cbind(dat_NOR, data.frame(.val=character(),
#                                                                            .equ=character(),
#                                                                            .zer=character(),
#                                                                            .cap=character(),
#                                                                            .cen=character(),
#                                                                            .sea=character(),
#                                                                            .con=character(),
#                                                                            .gbf=character(),
#                                                                            .inst=character(),
#                                                                            .dpl=character(),
#                                                                            .summary=character(),stringsAsFactors=FALSE))
#                                    }
#                                      
#  flags<-rbind(flags_NO_NOR,flags_NOR)
#  assign(paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_"),flags)
#}
#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*")))


### Cleaning the data by clean_coordinates function, with experimenting the data in Norway
dat <- GBIF_Data

dat$decimallongitude <- as.numeric(levels(dat$decimallongitude))[dat$decimallongitude]

dat$decimallatitude<- as.numeric(levels(dat$decimallatitude))[dat$decimallatitude]

dat <- dat %>% filter(!is.na(decimallongitude)) %>% filter(!is.na(decimallatitude))

dat$countryCode <- countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')

dat_NO_NOR <- dat %>% filter(countryCode != "NOR")

dat_NOR <- dat %>% filter(countryCode == "NOR")

dat_NO_NOR <-dat_NO_NOR %>% filter(-180 < decimallongitude & decimallongitude <180 &
                                     -90 < decimallatitude & decimallatitude < 90)

flags_NO_NOR <- clean_coordinates(x = dat_NO_NOR, lon = "decimallongitude", lat = "decimallatitude",
                                  countries = "countryCode", 
                                  species = "species",
                                  tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                            "zeros", "seas", "duplicates"),   
                                  seas_scale = 110)

dat_NOR <-dat_NOR %>% filter(-180 < decimallongitude & decimallongitude < 180 &
                               -90 < decimallatitude & decimallatitude < 90)

flags_NOR <- clean_coordinates(x = dat_NOR, lon = "decimallongitude", lat = "decimallatitude",
                               countries = "countryCode", 
                               species = "species",
                               tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                         "zeros","seas", "duplicates"),
                               seas_scale = 110)

Clean_Data <- rbind(flags_NO_NOR,flags_NOR)

Flagged_Data <- Clean_Data %>% filter(.summary == FALSE)

Clean_Data <- Clean_Data %>% filter(.summary == TRUE)

##Patrick code
##Just showing an alternative method for automated removal of records.

#for(i in 1){
#dat<-GBIF_Data[[paste(Bee_Species[i,1],"GBIF_Data", sep = "_")]]%>%
#  filter(!is.na(decimallongitude))%>%
#  filter(!is.na(decimallatitude))
#dat$countryCode <-  countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')
#  clean_dat<- data.frame(dat) %>%
#  cc_val() %>%                      
#  cc_equ() %>%                      
#  cc_cap() %>%                      
#  cc_cen() %>%                      
#  cc_coun(iso3 = "countryCode") %>% 
#  cc_gbif() %>%                     
#  cc_inst() %>%                      
#  cc_sea() %>%                      
#  cc_zero() %>%                     
#  cc_dupl()
#assign(paste(Bee_Species[i,1],"clean_Data",sep = "_"),clean_dat)

#for(i in 1:length(GBIF_Data)){
#Flagged_Data<-Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]] %>%
#  filter(.summary==FALSE)
#Data_cleaned<-Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]] %>%
#  filter(.summary==TRUE)
#assign(paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_"),Data_cleaned)
#assign(paste(GBIF_Data[[i]][["species"]][1],"flagged_Data",sep = "_"),Flagged_Data)
#}



#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*")))
#Flagged_Data<-setNames(lapply(ls(pattern ="flagged_Data*"), function(x) get(x)),(ls(pattern="flagged_Data*")))


### Filter out the records that have low precision

#for(i in 1:length(GBIF_Data)){
#if("coordinateUncertaintyInMeters" %in% colnames(Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]])){
#  clean<-Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]] %>%
#  filter(coordinateUncertaintyInMeters/1000 <= 5 | is.na(coordinateUncertaintyInMeters))
#  assign(paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_"),clean)
#}} 

#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*")))


### Cleaning by metadata
dat_NO_COOR <- Clean_Data %>% filter(coordinateUncertaintyInMeters == "")

dat_COOR <- Clean_Data %>% filter(coordinateUncertaintyInMeters != "")

dat_COOR$coordinateUncertaintyInMeters <- as.numeric(levels(dat_COOR$coordinateUncertaintyInMeters))[dat_COOR$coordinateUncertaintyInMeters]

dat_COOR <- dat_COOR %>% filter(coordinateUncertaintyInMeters/1000 <= 5 | is.na(coordinateUncertaintyInMeters))

dat_COOR$coordinateUncertaintyInMeters <- as.factor(dat_COOR$coordinateUncertaintyInMeters)

Clean_Data <- rbind(dat_NO_COOR, dat_COOR)

dat_NO_BAS  <- Clean_Data %>% filter(basisOfRecord == "") 

dat_BAS  <- Clean_Data %>% filter(basisOfRecord != "")

Clean_Data <- dat_BAS %>% filter( basisOfRecord == "HUMAN_OBSERVATION" |
                                    basisOfRecord == "OBSERVATION" |
                                    basisOfRecord == "PRESERVED_SPECIMEN")

dat_NO_IND <- Clean_Data %>% filter(individualCount =="")

dat_IND <- Clean_Data %>% filter(individualCount !="")

dat_IND$individualCount <- as.numeric(levels(dat_IND$individualCount))[dat_IND$individualCount]

dat_IND <- dat_IND %>% filter(individualCount > 0 & individualCount < 99)

dat_IND$individualCount <- as.factor(dat_IND$individualCount)

Clean_Data <- rbind(dat_NO_IND, dat_IND)

dat_NO_YEAR <- Clean_Data %>% filter(year == "")

dat_YEAR <- Clean_Data %>% filter(year != "")

dat_YEAR$year <- as.numeric(levels(dat_YEAR$year))[dat_YEAR$year]

dat_YEAR <- dat_YEAR %>% filter(year > 1945)

dat_YEAR$year <- as.factor(dat_YEAR$year)

Clean_Data <- dat_YEAR

Clean_Data <- Clean_Data[ , 1:13] 

Taxon_lim <- Clean_Data %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

Taxon_lim <- Taxon_lim %>% filter(COUNT >25)

Clean_Data <- Clean_Data %>% filter(Clean_Data$taxonKey %in% Taxon_lim$taxonKey)

Taxon_lim <- as.data.frame(Taxon_lim)

Clean_Data <- Clean_Data %>% filter(decimallatitude > -89)

CC_Data <- c()
for(i in 1:nrow(Taxon_lim)){
  GBIF <- Clean_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  GBIF <- GBIF %>% cc_outl(mltpl = 12)
  CC_Data <- rbind(CC_Data, GBIF)
  print(i)
}

Clean_Data <- CC_Data

save(Clean_Data, file = "../output/RData_European_Bee_Species_Clean_Data.RData") ## save clean data


##Patrick code
## Then a loop to remove records that are from unsuitable data sources e.g fossils 

#for(i in 1:length(GBIF_Data)){
#if("basisOfRecord" %in% colnames(Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]])){
#  clean<-Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]] %>%
#    filter( basisOfRecord == "HUMAN_OBSERVATION" |
#              basisOfRecord == "OBSERVATION" |
#              basisOfRecord == "PRESERVED_SPECIMEN")
#  assign(paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_"),clean)
#}}

#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*")))


## filter out records with absence records and high individual count records also

#for(i in 1:length(GBIF_Data)){
#if("individualCount" %in% colnames(Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]])){
#  clean<-Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]] %>%
#    filter(individualCount > 0 | is.na(individualCount)) %>%
#    filter(individualCount < 99 | is.na(individualCount))
#  assign(paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_"),clean)
#}} 

#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*")))

## filter out the old collection records

#for(i in 1:length(GBIF_Data)){
#if("year" %in% colnames(Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]])){
#  clean<-Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]] %>%
#    filter(year > 1945)
#  assign(paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_"),clean)
#}}

#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*")))  

## Now most of the records have been cleaned we can re-run an outlier analysis to remove some final spurious records -- species with less than 7 unique 
## records cannot be tested

#for(i in 1:length(GBIF_Data)){
#if((NROW(Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]])) >= 7){
#clean<- Clean_Data[[paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_")]] %>%
#      cc_outl(mltpl = 12)
#    assign(paste(GBIF_Data[[i]][["species"]][1],"clean_Data",sep = "_"),clean)
#}}

#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*")))

#rm(list=ls())                                                    ## Clear Environment

#load("RData_European_Bee_Species_Clean_Data.RData")              ## Reload data after automated cleaning


# finally remove species that have less than 25 occurrences remaining after inital cleaning


#for(i in 1:length(Clean_Data)){
#  if ((NROW(Clean_Data[[paste(Clean_Data[[i]][["species"]][1],"clean_Data",sep = "_")]])) >= 20){ 
#    assign(paste(Clean_Data[[i]][["species"]][1],"clean_Data",sep = "_"),Clean_Data[[paste(Clean_Data[[i]][["species"]][1],"clean_Data",sep = "_")]])
#  }}

#Clean_Data<-setNames(lapply(ls(pattern ="clean_Data*"), function(x) get(x)),(ls(pattern="clean_Data*"))) ### 405 species remain with over 25 occurences after inital cleaning 
#save(Clean_Data, file = "RData_European_Bee_Species_Clean_Data.RData") ## Save the clean data RData
