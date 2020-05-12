rm(list=ls())
require(dplyr)
require(tidyr)
require(ggplot2)
require(doParallel)

load("../output/RData_European_Bee_Species_GBIF_Data.RData")
load("../output/RData_European_Bee_Species_Clean_Data.RData")
load("../output/RData_European_Bee_Species_Finish_Data.RData")
load("../output/RData_European_Bee_Species_Thinned_Data_5km.RData")

wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify()

Taxon_lim <- Thinned_Data_5km %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)


### Creating maps before and after cleaning
for(i in 1:3){
  GBIF <- GBIF_Data %>% filter(taxonKey == Taxon_lim[i, 2])
            
            GBIF$decimallongitude <- as.numeric(levels(GBIF$decimallongitude))[GBIF$decimallongitude]
            
            GBIF$decimallatitude <- as.numeric(levels(GBIF$decimallatitude))[GBIF$decimallatitude]
            
            Clean <- Clean_Data %>% filter(taxonKey == Taxon_lim[i, 2])
            
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
            
            ggsave(paste(GBIF$species,"_1_before.png",sep = ""), plot = before, device = "png", path="~/Documents/Patrick_Masa_code/SDMs/figs/clean_maps/", height = 10, width = 10, dpi=300) ## save the before, after, and polish maps
            
            ggsave(paste(Clean$species,"_2_after.png",sep = ""), plot = after, device = "png", path="~/Documents/Patrick_Masa_code/SDMs/figs/clean_maps/", height = 10, width = 10, dpi=300)
            
            print(i)
            
            }


### Creating maps for finished data
for(i in 1:3){
  GBIF <- GBIF_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  GBIF$decimallongitude <- as.numeric(levels(GBIF$decimallongitude))[GBIF$decimallongitude]
  
  GBIF$decimallatitude <- as.numeric(levels(GBIF$decimallatitude))[GBIF$decimallatitude]
  
  Finish <- Finish_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  polish<-ggplot()+ coord_fixed()+
    
    geom_map(data =wm, map = wm,
             
             aes(group = group, map_id= region),
             
             fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
    
    geom_point(data = Finish, aes(x = decimallongitude, y = decimallatitude),
               
               colour = "darkred", size = 1)+
    
    lims(x = c(floor(min(GBIF$decimallongitude)-10),     ### limit the extent of the map 
               
               ceiling(max(GBIF$decimallongitude)+10)),  ### to capture the points
         
         y = c(floor(min(GBIF$decimallatitude)-10),
               
               ceiling(max(GBIF$decimallatitude)+10)))+
    
    theme_bw()
  
  ggsave(paste(Finish$species,"_3_polish.png",sep = ""), plot = polish, device = "png", path="~/Documents/Patrick_Masa_code/SDMs/figs/clean_maps/", height = 10, width = 10, dpi=300) ## save the before, after, and polish maps
  
  print(i)
  
  }


### Creating maps for thinned data at 5km
for(i in 1:3){
  GBIF <- GBIF_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  GBIF$decimallongitude <- as.numeric(levels(GBIF$decimallongitude))[GBIF$decimallongitude]
  
  GBIF$decimallatitude <- as.numeric(levels(GBIF$decimallatitude))[GBIF$decimallatitude]
  
  Thinned_5km <- Thinned_Data_5km %>% filter(taxonKey == Taxon_lim[i, 2])
  
  thinned<-ggplot()+ coord_fixed()+
    
    geom_map(data =wm, map = wm,
             
             aes(group = group, map_id= region),
             
             fill = "darkgrey")+     ## if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
    
    geom_point(data = Thinned_5km, aes(x = decimallongitude, y = decimallatitude),
               
               colour = "darkred", size = 1)+
    
    lims(x = c(floor(min(GBIF$decimallongitude)-10),     ### limit the extent of the map 
               
               ceiling(max(GBIF$decimallongitude)+10)),  ### to capture the points
         
         y = c(floor(min(GBIF$decimallatitude)-10),
               
               ceiling(max(GBIF$decimallatitude)+10)))+
    
    theme_bw()
  
  ggsave(paste(Thinned_5km$species,"_4_Thinned.png",sep = ""), plot = thinned, device = "png", path="~/Documents/Patrick_Masa_code/SDMs/figs/clean_maps/", height = 10, width = 10, dpi=300) ## save the before, after, and polish maps
  
  print(i)
  
}

##Patrick code
## Can see if we can create some before and after comparison maps.
#wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify() # This gets the polygon data for the world map


#registerDoParallel( cores = 4 )


#foreach(i = 1:length(Clean_Data),
#          .packages = c("dplyr","ggplot2","maps","tidyr")) %dopar% {
#  before<-ggplot()+ coord_fixed()+
#    geom_map(data =wm, map = wm,
#             aes(group = group, map_id= region),
#             fill = "darkgrey")+                        ### if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
#    geom_point(data = GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]], aes(x = decimallongitude, y = decimallatitude),
#               colour = "darkred", size = 1)+
#    lims(x = c(floor(min(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]])-10),     ### limit the extent of the map 
#               ceiling(max(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]]))+10),  ### to capture the points
#         y = c(floor(min(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]])-10),
#               ceiling(max(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]]))+10))+
#    theme_bw()
#  after<-ggplot()+ coord_fixed()+
#    geom_map(data =wm, map = wm,
#             aes(group = group, map_id= region),
#             fill = "darkgrey")+
#    geom_point(data = Clean_Data[[i]], aes(x = decimallongitude, y = decimallatitude),
#               colour = "blue", size = 1)+
#    lims(x = c(floor(min(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]])-10),
#               ceiling(max(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]]))+10), 
#         y = c(floor(min(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]])-10),
#               ceiling(max(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]]))+10))+
#    theme_bw()
#  ggsave(paste((paste("Cleaning_6_Cleaning_Maps/",Clean_Data[[i]][["species"]][1],sep="")),"/1_before.png",sep = ""),before,device = "png", height = 10, width = 10,dpi=300) ## save the before, after, and polish maps 
#  ggsave(paste((paste("Cleaning_6_Cleaning_Maps/",Clean_Data[[i]][["species"]][1],sep = "")),"/2_after.png",sep = ""),after,device = "png", height = 10, width = 10,dpi=300)
#}

#registerDoSEQ()

## Now for the Finished cleaning map

#registerDoParallel( cores = 4 )

#foreach(i = 1:length(Finish_Data),
#        .packages = c("dplyr","ggplot2","maps","tidyr")) %dopar% {
#polish<-ggplot()+ coord_fixed()+
#  geom_map(data =wm, map = wm,
#           aes(group = group, map_id= region),
#           fill = "darkgrey")+                        ### if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
#  geom_point(data = Finish_Data[[i]], aes(decimallongitude, decimallatitude),
#             colour = "green", size = 1)+
#  lims(x = c(floor(min(GBIF_Data[[paste(Finish_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]])-10),     ### limit the extent of the map 
#             ceiling(max(GBIF_Data[[paste(Finish_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]]))+10),  ### to capture the points
#       y = c(floor(min(GBIF_Data[[paste(Finish_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]])-10),
#             ceiling(max(GBIF_Data[[paste(Finish_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]]))+10))+
#  theme_bw()
#ggsave(paste((paste("Cleaning_6_Cleaning_Maps/",Finish_Data[[i]][["species"]][1],sep="")),"/4_polish.png",sep = ""),polish,device = "png", height = 10, width = 10,dpi=300)
#}

#registerDoSEQ()
## final Map showing the Thinned data at 5km 

#for(i in 1:length(Thinned_Data_5km)){
#  thinned<-ggplot()+ coord_fixed()+
#    geom_map(data =wm, map = wm,
#             aes(group = group, map_id= region),
#             fill = "darkgrey")+                        ### if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
#    geom_point(data = Thinned_Data_5km[[i]], aes(decimallongitude, decimallatitude),
#               colour = "orange", size = 1)+
#    lims(x = c(floor(min(GBIF_Data[[paste(Thinned_Data_5km[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]])-10),     ### limit the extent of the map 
#               ceiling(max(GBIF_Data[[paste(Thinned_Data_5km[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]]))+10),  ### to capture the points
#         y = c(floor(min(GBIF_Data[[paste(Thinned_Data_5km[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]])-10),
#               ceiling(max(GBIF_Data[[paste(Thinned_Data_5km[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]]))+10))+
#    theme_bw()
#  ggsave(paste((paste("Cleaning_6_Cleaning_Maps/",Thinned_Data_5km[[i]][["species"]][1],sep="")),"/8_Thinned.png",sep = ""),thinned,device = "png", height = 10, width = 10,dpi=300)
#}

