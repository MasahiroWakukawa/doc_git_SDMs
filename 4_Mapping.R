rm(list=ls())
require(dplyr)
require(tidyr)
require(ggplot2)
library(ggmap)
library(ggthemes)
library(viridis) # devtools::install_github("sjmgarnier/viridis)
library(scales)
library(grid)
library(gridExtra)


load("../output/RData_Plant_Species_Data.RData")

load("../output/RData_Plant_Species_Clean_Data.RData")

load("../output/RData_Plant_Species_Thinned_Data_5km.RData")

wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify()

Taxon_lim <- Thinned_Data_5km %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

Species_count <- Clean_Data %>% group_by(species) %>% summarise(COUNT=n()) 

### Comparing maps before and after cleaning

for(i in 1:nrow(Taxon_lim)){
  GBIF <- Plant %>% filter(species == Taxon_lim[i, 1])
  
  Clean <- Clean_Data %>% filter(species == Taxon_lim[i, 1])
  
  before<-ggplot()+ coord_fixed()+
    
    geom_map(data =wm, map = wm,
             
             aes(group = group, map_id= region),
             
             fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
    
    geom_point(data = GBIF, aes(x = decimallongitude, y = decimallatitude),
               
               colour = "darkred", size = 1)+
    
    lims(x = c(floor(min(GBIF$decimallongitude)-10),     ### limit the extent of the map 
               
               ceiling(max(GBIF$decimallongitude)+10)),  ### to capture the points
         
         y = c(floor(min(GBIF$decimallatitude)-10),
               
               ceiling(max(GBIF$decimallatitude)+10)))+
    
    theme_bw()
  
  after<-ggplot()+ coord_fixed()+
    
    geom_map(data =wm, map = wm,
             
             aes(group = group, map_id= region),
             
             fill = "darkgrey")+
    
    geom_point(data = Clean, aes(x = decimallongitude, y = decimallatitude),
               
               colour = "blue", size = 1)+
    
    lims(x = c(floor(min(GBIF$decimallongitude)-10),
               
               ceiling(max(GBIF$decimallongitude)+10)), 
         
         y = c(floor(min(GBIF$decimallatitude)-10),
               
               ceiling(max(GBIF$decimallatitude)+10)))+
    
    theme_bw()
  
  ggsave(paste(GBIF$species,"_1_before.png",sep = ""), plot = before, device = "png", path="~/Documents/Patrick_Masa_code/SDMs/figs/clean_maps_plant/", height = 10, width = 10, dpi=300) ## save the before, after, and polish maps
  
  ggsave(paste(Clean$species,"_2_after.png",sep = ""), plot = after, device = "png", path="~/Documents/Patrick_Masa_code/SDMs/figs/clean_maps_plant/", height = 10, width = 10, dpi=300)
  
  print(i)
  
}


### Creating maps for thinned data at 5km
for(i in 1:nrow(Taxon_lim)){
  
  Thinned_5km <- Thinned_Data_5km %>% filter(species == Taxon_lim[i, 1])
  
  thinned<-ggplot()+ coord_fixed()+
    
    geom_map(data =wm, map = wm,
             
             aes(group = group, map_id= region),
             
             fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
    
    geom_point(data = Thinned_5km, aes(x = decimallongitude, y = decimallatitude),
               
               colour = "darkred", size = 1)+
    
    lims(x = c(floor(min(Thinned_5km$decimallongitude)-10),     ### limit the extent of the map 
               
               ceiling(max(Thinned_5km$decimallongitude)+10)),  ### to capture the points
         
         y = c(floor(min(Thinned_5km$decimallatitude)-10),
               
               ceiling(max(Thinned_5km$decimallatitude)+10)))+
    
    theme_bw()
  
  ggsave(paste(Thinned_5km$species,"_4_Thinned.png",sep = ""), plot = thinned, device = "png", path="~/Documents/Patrick_Masa_code/SDMs/figs/clean_maps_plant/", height = 10, width = 10, dpi=300) ## save the before, after, and polish maps
  
  print(i)
  
}

ggplot()+ coord_fixed()+
  
  geom_map(data =wm, map = wm,
           
           aes(group = group, map_id= region),
           
           fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
  
  geom_point(data = Thinned_Data_5km, aes(x = decimallongitude, y = decimallatitude),
             
             colour = "darkred", size = 1)+
  
  lims(x = c(floor(min(Thinned_Data_5km$decimallongitude)-10),     ### limit the extent of the map 
             
             ceiling(max(Thinned_Data_5km$decimallongitude)+10)),  ### to capture the points
       
       y = c(floor(min(Thinned_Data_5km$decimallatitude)-10),
             
             ceiling(max(Thinned_Data_5km$decimallatitude)+10)))+
  
  theme_bw()


### Create a heatmap of thinned data (Not looking good and improvement is necessary)

wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify()

ggplot()+ coord_fixed()+geom_map(data =wm, map = wm,
                                 
                                 aes(group = group, map_id= region),
                                 
                                 fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
  
  geom_point(data = Thinned_Data_5km, aes(x = decimallongitude, y = decimallatitude),
             
             colour = "darkred", size = 1)

##1st option
qmplot(data = Thinned_Data_5km, x = decimallongitude, y = decimallatitude, maptype = "toner-lite", geom = "density2d", color = I("red"))+
  
  stat_density2d(aes(fill = log(..level..)), geom = "polygon", alpha = 0.3, color=NA)+

  scale_fill_gradient2("Density", low = "white", mid = "yellow", high = "red", midpoint = -20)



##2nd option
ggplot()+ coord_fixed()+geom_map(data =wm, map = wm,
                                 
                                 aes(group = group, map_id= region),
                                 
                                 fill = "darkgrey")+
  
  geom_density2d(data = Thinned_Data_5km, mapping = aes(x = decimallongitude, y = decimallatitude), size = 0.3)+ 
  
  stat_density2d(data = Thinned_Data_5km,
                 
                 aes(x = decimallongitude, y = decimallatitude, fill = log(..level..), alpha = ..level..), size = 0.01,
                 
                 bins = 16, geom = "polygon") +
  
  scale_fill_gradient(low = "green", high = "red")+
  
  scale_alpha(range = c(0, 0.3), guide = FALSE) 

