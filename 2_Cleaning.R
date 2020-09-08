rm(list=ls())
require(dplyr)
require(tidyr)
require(rgbif)
require(magrittr)
require(countrycode)
require(CoordinateCleaner)
require(rnaturalearthdata)
require(rgdal)
require(sp)
require(ggplot2)
require(adehabitatHR)
require(scales)
require(doParallel)
require(sf)
require(raster)



load("../output/RData_Plant_Species_Data.RData") 


### Cleaning the data by clean_coordinates function
dat <- Plant

dat <- dat %>% filter(!is.na(decimallongitude)) %>% filter(!is.na(decimallatitude))

dat <- dat %>% filter(-180 < decimallongitude & decimallongitude <180 &
                                     -90 < decimallatitude & decimallatitude < 90)

flags_dat <- clean_coordinates(x = dat, lon = "decimallongitude", lat = "decimallatitude",
                                  species = "species",
                                  tests = c("capitals", "centroids", "duplicates", "equal", "gbif", "institutions",
                                            "seas", "zeros"),   
                                  seas_scale = 110)

Clean_Data <- flags_dat

Flagged_Data <- Clean_Data %>% filter(.summary == FALSE)

Clean_Data <- Clean_Data %>% filter(.summary == TRUE)

head(Clean_Data)


### Cleaning the data with less than 20 records

Clean_Data <- Clean_Data[ , 1:6] 

Taxon_lim <- Clean_Data %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- Taxon_lim %>% filter(COUNT >20)

Clean_Data <- Clean_Data %>% filter(Clean_Data$species%in% Taxon_lim$species)

Taxon_lim <- as.data.frame(Taxon_lim)

Taxon_lim

### Cleaning data with cc_outl function
CC_Data <- c()
for(i in 1:nrow(Taxon_lim)){
  GBIF <- Clean_Data %>% filter(species == Taxon_lim[i, 1])
  GBIF <- GBIF %>% cc_outl(mltpl = 12)
  CC_Data <- rbind(CC_Data, GBIF)
  print(i)
}

Clean_Data <- CC_Data

save(Clean_Data, file = "../output/RData_Plant_Species_Clean_Data.RData")

