rm(list=ls())
require(dplyr)
require(tidyr)
require(rgdal)
require(sp)
require(ggplot2)
require(magrittr)
require(adehabitatHR)
require(scales)
require(CoordinateCleaner)
require(doParallel)

load("../output/RData_European_Bee_Species_Clean_Data.RData")

### Polishing data by cc_iucn function
Taxon_lim <- Clean_Data %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

CC_Data <- c()

for(i in 1:nrow(Taxon_lim)){
  
  GBIF <- Clean_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  GBIF$species <- droplevels(GBIF$species)
  
  Coordinates<-GBIF[,c("decimallongitude","decimallatitude")]
  
  SpatialCoord<-SpatialPointsDataFrame(coords = Coordinates, data = GBIF,
                                       proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  
  SpatialCoord <- SpatialCoord[ , c("species")]
  
  names(SpatialCoord)[1]<-"id"
  
  Bee_mcp<-mcp(SpatialCoord, percent = 95)
  
  names(Bee_mcp)[1] <- "species"
  
  Final_clean <- cc_iucn(x=GBIF, range=Bee_mcp, value="flagged")
  
  Finish_clean<-GBIF[Final_clean,]
  
  CC_Data <- rbind(CC_Data, Finish_clean)
  
  print(i)
  
}

Finish_Data <- CC_Data

save(Finish_Data, file = "../output/RData_European_Bee_Species_Finish_Data.RData")

##Patrick code
#?dir.create
#for(i in 1:length(Clean_Data)){                                         ### Create a file pertaining to each species in a Maps folder
#  dir.create(paste("Cleaning_6_Cleaning_Maps/",Clean_Data[[i]][["species"]][1],sep = ""))
#}

#registerDoParallel(cores = 4)

#?registerDoParallel
#?.packages
#?doParallel
#?mcp.area
#?SpatialPointsDataFrame

#foreach(i = 1:length(Clean_Data),
#        .packages = c("sp","dplyr","ggplot2","adehabitatHR", "magrittr","rgdal","scales")) %dopar% {
#  Coordinates<-Clean_Data[[i]][,c("decimallongitude","decimallatitude")]
#SpatialCoord<-SpatialPointsDataFrame(coords = Coordinates, data = Clean_Data[[paste(Clean_Data[[i]][["species"]][1],"clean_Data",sep="_")]],   ##Convert regular data frame into
#                                           proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))         ## an SpatialPointsDataFrame    
#SpatialCoord<-SpatialCoord[,c("species")]  ## Coordinates now embedded within SPDF so for congruence with the package 
#names(SpatialCoord)[1]<-"id"               ## rename "species" column "id"
#assign(paste(Clean_Data[[i]][["species"]][1], "Spatial_Coordinates", sep = "_"),SpatialCoord)   ## Create a object in the environment of the Spatail
#jpeg(paste("Cleaning_6_Cleaning_Maps/",paste(Clean_Data[[i]][["species"]][1],"/5_Polygon_Area.jpeg",sep=""),sep = ""))## pointsDataFrame for each species <- this will be needed
#poly_area<-mcp.area(SpatialCoord, percent = seq(50,100, by = 0.1))                                       ## later
#dev.off()                                                                                                ## save the area plots in the same file with the maps
#}

#registerDoSEQ()

#Spatial_Data<-setNames(lapply(ls(pattern ="Spatial_Coordinates*"), function(x) get(x)),(ls(pattern="Spatial_Coordinates*"))) ## collate data in a list
#save(file = "RData_European_Bee_Species_Spatial_Data.RData", Spatial_Data)                                               ## Save the SpatialpointsDataframe

## having manually looked through the graphs and the point maps to determine what percent of points the mcp should use. 
#Scientific_Name <- c()
#for(i in 1:length(Clean_Data)){
#  Scientific_Name <- rbind(Scientific_Name,Clean_Data[[i]][["species"]][1])
#}

#write.csv(file = "European_Bee_Species_MCP_Area_Percent.csv", x =  data.frame(Scientific_Name), row.names = FALSE)  ##Raw_Data/Taxa/European_Bee_Species_MCP_Area_Percent.csv

#Bee_Percent<-read.csv("European_Bee_Species_MCP_Area_Percent.csv")  ### Load manual evaluation of the area plots    ##Raw_Data/Taxa/European_Bee_Species_MCP_Area_Percent.csv

#wm<-map_data("world") %>% filter(region != "Antartica") %>% fortify() # This gets the polygon data for the world map

#registerDoParallel(cores = 4)
#Clean_Data[[147]][["species"]][1]
#?foreach
#?mcp

#foreach(i = 1:length(Clean_Data),
#        .packages = c("sp","dplyr","ggplot2","adehabitatHR", "magrittr","rgdal","scales"),
#        .combine = "list" ) %dopar% {
#          if(NROW(Clean_Data[[i]]) > 5){
#Bee_mcp<-mcp(Spatial_Data[[paste(Clean_Data[[i]][["species"]][1],"Spatial_Coordinates",sep = "_")]],percent = Bee_Percent[i,2])
#Bee_points<-data.frame(Spatial_Data[[paste(Clean_Data[[i]][["species"]][1],"Spatial_Coordinates",sep = "_")]]@coords,
#                       id= Spatial_Data[[paste(Clean_Data[[i]][["species"]][1],"Spatial_Coordinates",sep = "_")]]@data$id) 
#assign(paste(Clean_Data[[i]][["species"]][1],"Minimum_Convex_Polygon", sep = "_"),Bee_mcp)}
#Polygon<-ggplot()+ coord_fixed()+
#  geom_map(data =wm, map = wm,
#           aes(group = group, map_id= region),
#           fill = "darkgrey")+                        ### if you want to add country borders  (colour= "#7f7f7f", size = 0.5)
#  geom_point(data = Bee_points, aes(decimallongitude, decimallatitude),
#             colour = "blue", size = 1)+
#  lims(x = c(floor(min(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]])-10),     ### limit the extent of the map 
#             ceiling(max(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallongitude"]]))+10),  ### to capture the points
#       y = c(floor(min(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]])-10),
#             ceiling(max(GBIF_Data[[paste(Clean_Data[[i]][["species"]][1],"GBIF_Data",sep = "_")]][["decimallatitude"]]))+10))+
#  geom_polygon(data = fortify(Bee_mcp),                                ### visulise the polygon around the points
#               aes(long,lat),alpha= 0.3)+
#theme_bw()
#ggsave(paste("Cleaning_6_Cleaning_Maps/",paste(Clean_Data[[i]][["species"]][1],"/3_Polygon.png",sep=""),sep = ""),Polygon,device = "png", height = 10, width = 10,dpi=300)
#}

#Bee_mcp<-mcp(Spatial_Data[[paste(Clean_Data[[1]][["species"]][1],"Spatial_Coordinates",sep = "_")]],percent = Bee_Percent[1,2])
#?mcp
#Clean_Data[[1]][["species"]][1]

#Polygon_Data<-setNames(lapply(ls(pattern ="Minimum_Convex_Polygon*"), function(x) get(x)),(ls(pattern="Minimum_Convex_Polygon*"))) ## collate all the spatial polygons together
#save(file = "RData_European_Bee_Species_Polygon_Data.RData", Polygon_Data)

### Finally overlay the polygon back onto the map, points that are not within the boundary of the polygon will be removed as the final 
### act of cleaning using the function cc_iucn.

#?cc_iucn

#for(i in 1:length(Polygon_Data)){
#names(Polygon_Data[[i]])[1]<-"species"  ## rename the id column back to species so that it is compatible with the package
#Final_clean<-cc_iucn(x=Clean_Data[[paste(as.character(Polygon_Data[[i]][["species"]][1]),"clean_Data",sep = "_")]],               ### this will remove points that are outside of the polygon
#                     range = Polygon_Data[[i]],
#                     value = "flagged")
#Finish_clean<-Clean_Data[[paste(as.character(Polygon_Data[[i]][["species"]][1]),"clean_Data",sep = "_")]][Final_clean,]
#assign(paste(Clean_Data[[i]][["species"]][1],"Finish_Clean", sep = "_"),Finish_clean)
#}                                                                                              


#Finish_Data<-setNames(lapply(ls(pattern ="Finish_Clean*"), function(x) get(x)),(ls(pattern="Finish_Clean*")))
#save(Finish_Data, file = "RData_European_Bee_Species_Finish_Data.RData") 

##Masa test code
#Taxon_lim <- Clean_Data %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

#Taxon_lim <- as.data.frame(Taxon_lim)

#GBIF <- Clean_Data %>% filter(taxonKey == Taxon_lim[1, 2])

#GBIF$species <- droplevels(GBIF$species)

#Coordinates<-GBIF[,c("decimallongitude","decimallatitude")]

#SpatialCoord<-SpatialPointsDataFrame(coords = Coordinates, data = GBIF,
#                                      proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

#SpatialCoord <- SpatialCoord[ , c("species")]

#names(SpatialCoord)[1]<-"id"

#Bee_mcp<-mcp(SpatialCoord, percent = 95)

#names(Bee_mcp)[1] <- "species"

#Final_clean <- cc_iucn(x=GBIF, range=Bee_mcp, value="flagged")

#Finish_clean<-GBIF[Final_clean,]
