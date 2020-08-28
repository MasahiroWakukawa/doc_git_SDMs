rm(list=ls())
require(dplyr)
require(tidyr)
require(rgbif)
require(taxize)
require(readxl)
require(countrycode)
require(ggplot2)
require(rgeos)
require(maptools)
library(maps)



### Reading real data of plant species

Plant <- read_excel("../data/monocot_specimens.xlsx", range = "A1:AR61721")

Plant <- cbind(Plant$family, Plant$species, Plant$decimallatitude, Plant$decimallongitude, Plant$taxonkey, Plant$countrycode)

Plant <- as.data.frame(Plant)

Plant2 <- read_excel("../data/monocot_specimens2.xlsx", range = "A1:AC89999")

Plant2 <- cbind(Plant2$family, Plant2$species, Plant2$decimallatitude, Plant2$decimallongitude, Plant2$taxonkey, Plant2$countrycode)

Plant2 <- as.data.frame(Plant2)

Plant3 <- read_excel("../data/monocot_specimens3.xlsx", range = "A1:AC90002")

Plant3 <- cbind(Plant3$family, Plant3$species, Plant3$decimallatitude, Plant3$decimallongitude, Plant3$taxonkey, Plant3$countrycode)

Plant3 <- as.data.frame(Plant3)

Plant4 <- read_excel("../data/monocot_specimens4.xlsx", range = "A1:AC51659")

Plant4 <- cbind(Plant4$family, Plant4$species, Plant4$decimallatitude, Plant4$decimallongitude, Plant4$taxonkey, Plant4$countrycode)

Plant4 <- as.data.frame(Plant4)

Plant_all <- rbind(Plant, Plant2, Plant3, Plant4)

colnames(Plant_all)[1:6] <- c("family", "species", "decimallatitude", "decimallongitude", "taxonkey", "countrycode")

count_spcies <- Plant_all %>% group_by(species) %>% summarise(COUNT=n())

wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify()

ggplot()+ coord_fixed()+geom_map(data =wm, map = wm,
                                 
                                 aes(group = group, map_id= region),
                                 
                                 fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
  
  geom_point(data = Plant_all, aes(x = decimallongitude, y = decimallatitude),
             
             colour = "darkred", size = 1)

Plant$countrycode <- countrycode(Plant$countrycode, origin =  'iso2c', destination = 'iso3c')

Plant<- Plant %>% filter(species!="NA" & decimallatitude!="NA" & decimallongitude!="NA")

Plant_TaxonKey <- data.frame(Plant$species, Plant$taxonkey)

Plant_Occurrence_Count <- Plant_TaxonKey %>% group_by_all() %>% summarise(COUNT=n())

head(Plant_Occurrence_Count)

Plant_Occurrence_Count <- Plant_Occurrence_Count %>% filter(COUNT > 20)

Plant <- Plant %>% filter(taxonkey %in% Plant_Occurrence_Count$Plant.taxonkey)

Plant$decimallatitude <- as.numeric(levels(Plant$decimallatitude))[Plant$decimallatitude]

Plant$decimallongitude <- as.numeric(levels(Plant$decimallongitude))[Plant$decimallongitude]

Plant$taxonkey <- as.numeric(levels(Plant$taxonkey))[Plant$taxonkey]

save(file = "../output/RData_Plant_Species_Data.RData", Plant)


