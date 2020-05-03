rm(list = ls())
require(dplyr)
require(tidyr)
require(ape)
require(raster)
require(sp)
require(rgdal)
require(spThin)   ####jeffreyhanson branch for githbud download 
require(tiff)
require(gurobi)
require(doParallel) 


load("RData_European_Bee_Species_Finish_Data.RData")                          ## load finish Data for the thinning process
load("RData_EUropean_Bee_Species_GBIF_Data.RData")                            ## load GBIF_Data for maps to visualise effects of thinning


##Patrick code
registerDoParallel(cores= 4)

Thinned_Data_5km <- foreach(i = 1:length(Finish_Data),
                            .packages = c("dplyr","tidyr","ape","raster","sp","rgdal","spThin","tiff","gurobi")) %dopar% {
                              Coordinates <- Finish_Data[[i]][,c("decimallongitude","decimallatitude")]   ### extract the coordinate columns from Finish Data
                              Thin <- spThin(Coordinates,                                 ## "spThin" package will perform the thinning, this will randomly remove
                                             x.col = "decimallongitude",                  ## points that are within a certain distance from each other -- the algorithm 
                                             y.col = "decimallatitude",                   ## "gurobi will then run to identify the maximum number of points that are at least
                                             dist = 5000,                                 ## that distance apart to leave an optimum dataset. -- This takes a while to run 
                                             method = "gurobi",                           ## so be wary.
                                             great.circle.distance = TRUE
                              )
                              t1 <- Thin[[1]]                                    ## Extract the dataframe for the output
                              Coords_Thin <- data.frame(t1@coords)              ## Extract Coordinates 
                              Coords_Thin$min_Dist <- Thin@mindist                                                          ## Extract minimum distance points are from each other 
                              assign(paste(Finish_Data[[i]][["species"]][1],"Thinned_Data_5km",sep = "_"),Coords_Thin)   #### assign dataframe to environment
                            }

registerDoSEQ()

for(i in 1:length(Thinned_Data_5km)){
  names(Thinned_Data_5km)[i] <- paste(Finish_Data[[i]]["species"][1],"Thinned_Data_5km",sep = "_")
}

for(i in 1:length(Thinned_Data_5km)){                                           ## Change coordinate column names to match previous datasets
  colnames(Thinned_Data_5km[[i]])[1:2]<-c("decimallongitude","decimallatitude")
  Thinned_Data_5km[[i]]$species <- Finish_Data[[i]][["species"]][1]
}

for(i in 1:length(Thinned_Data_5km)){                   #### Remove those cases that have less that 25 points 
  if(nrow(Thinned_Data_5km[[i]]) > 25){
    assign(paste(Thinned_Data_5km[[i]][["species"]][1],"Thinned_Data_5km",sep = "_"),Thinned_Data_5km[[i]])
  }
}

rm(Thinned_Data_5km)
Thinned_Data_5km<-setNames(lapply(ls(pattern ="Thinned_Data_5km*"), function(x) get(x)),(ls(pattern="Thinned_Data_5km*")))




save(file = "RData_European_Bee_Species_Thinned_Data_5km.RData",Thinned_Data_5km)     ### Save 


for(i in 1:length(Finish_Data)){                                               ### Repeat for 10km minimum distance apart -- This took 6 hours. 
  Coordinates <- Finish_Data[[i]][,c("decimallongitude","decimallatitude")]
  Thin <- spThin(Coordinates,
                 x.col = "decimallongitude",
                 y.col = "decimallatitude",
                 dist = 10000,
                 method = "gurobi",
                 great.circle.distance = TRUE
  )
  t1 <- Thin[[1]]
  Coords_Thin <- data.frame(t1@coords)
  Coords_Thin$min_Dist <- Thin@mindist
  assign(paste(Finish_Data[[i]][["species"]][1],"Thinned_Data_10km",sep = "_"),Coords_Thin)
}

Thinned_Data_10km<-setNames(lapply(ls(pattern ="Thinned_Data_10km*"), function(x) get(x)),(ls(pattern="Thinned_Data_10km*")))


for(i in 1:length(Thinned_Data_10km)){
  colnames(Thinned_Data_10km[[i]])[1:2]<-c("decimallongitude","decimallatitude")
  Thinned_Data_10km[[i]]$species <- Finish_Data[[i]][["species"]][1]
}

for(i in 1:length(Thinned_Data_10km)){                   #### Remove those cases that have less that 25 points 
  if(nrow(Thinned_Data_10km[[i]]) > 25){
    assign(paste(Thinned_Data_10km[[i]][["species"]][1],"Thinned_Data_10km",sep = "_"),Thinned_Data_10km[[i]])
  }
}
Thinned_Data_10km<-setNames(lapply(ls(pattern ="Thinned_Data_10km*"), function(x) get(x)),(ls(pattern="Thinned_Data_10km*")))


save(file = "RData_European_Bee_Species_Thinned_Data_10km.RData",Thinned_Data_10km)

Taxon_lim <- Finish_Data %>% group_by(species, taxonKey) %>% summarise(COUNT=n()) 

##Masa code test
Taxon_lim <- as.data.frame(Taxon_lim)

GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[1, 2])

Coord <- GBIF[, c("species", "decimallongitude","decimallatitude")]

Thin <- thin(Coord, lat.col = "decimallatitude", long.col = "decimallongitude", 
             spec.col = "species", thin.par = 5, reps = 10)

?thin

##Masa future code 
CC_Data <- c()

for(i in 1:nrow(Taxon_lim)){
  
  GBIF <- Finish_Data %>% filter(taxonKey == Taxon_lim[i, 2])
  
  Coordinates<-GBIF[,c("decimallongitude","decimallatitude")]
  
  Thin <- spThin(Coord,                                 
                 x.col = "decimallongitude",                 
                 y.col = "decimallatitude",                   
                 dist = 5000,                                 
                 method = "gurobi",                           
                 great.circle.distance = TRUE
  )
  
  t1 <- Thin[[1]]
  
  Coords_Thin <- data.frame(t1@coords)              
  
  Coords_Thin$min_Dist <- Thin@mindist   
  
  Final_clean <- cc_iucn(x=GBIF, range=Bee_mcp, value="flagged")
  
  Finish_clean<-GBIF[Final_clean,]
  
  CC_Data <- rbind(CC_Data, Finish_clean)
  
  print(i)
  
}

