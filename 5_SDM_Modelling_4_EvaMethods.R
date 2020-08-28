rm(list= ls())
#require(biomod2)
require(dismo)
require(raster)
require(magrittr)
require(foreach)
require(doParallel)
require(ff)
require(SSDM)
require(dplyr)
require(ggplot2)
require(rgeos)
require(maptools)
library(knitr)
library(kableExtra)
library(rJava)
require(readxl)

load("../output/RData_Plant_Species_Thinned_Data_5km.RData")     

Thinned_Data_5km <- Thinned_Data_5km %>% filter(-175 < decimallongitude & decimallongitude <175 &
                                                  -89 < decimallatitude & decimallatitude < 89)

Taxon_lim <- Thinned_Data_5km %>% group_by(species) %>% summarise(COUNT=n()) 

Taxon_lim <- as.data.frame(Taxon_lim)

Taxon_lim

#Bioclim 1 : Annual mean temperature
#Bioclim 4 : temperature seasonality
#Bioclim 12 : annual precipitation
#Bioclim 15 : Precipitation seasonality

Bioclim_01 <- raster("../data/wc2/wc2.1_2.5m_bio_1.tif")

Bioclim_04 <- raster("../data/wc2/wc2.1_2.5m_bio_4.tif")

Bioclim_12 <- raster("../data/wc2/wc2.1_2.5m_bio_12.tif")

Bioclim_15 <- raster("../data/wc2/wc2.1_2.5m_bio_15.tif")   

Bioclim_Stack <- stack(Bioclim_01,Bioclim_04,Bioclim_12,Bioclim_15)

Occurrence_Data_5km <- Thinned_Data_5km[c("species", "decimallongitude", "decimallatitude")]

Occurrence_Data_5km$species <- as.character(Occurrence_Data_5km$species)

##Producing GLM table

GLM_table <- c()

for(i in 1:nrow(Taxon_lim)){

GLM_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GLM_table",sep=""),sep = ""))

GLM_table <- rbind(GLM_table, GLM_table_ind)

GLM_table <- c(mean(GLM_table$accuracy), mean(GLM_table$discrimination), mean(GLM_table$calibration), mean(GLM_table$precision))

}

##Producing GAM table

GAM_table <- c()

for(i in 1:nrow(Taxon_lim)){
  
  GAM_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GAM_table",sep=""),sep = ""))
  
  GAM_table <- rbind(GAM_table, GAM_table_ind)
  
  GAM_table <- c(mean(GAM_table$accuracy), mean(GAM_table$discrimination), mean(GAM_table$calibration), mean(GAM_table$precision))
  
}

##Producing MARS table

MARS_table <- c()

for(i in 1:nrow(Taxon_lim)){
  
  MARS_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/MARS_table",sep=""),sep = ""))
  
  MARS_table <- rbind(MARS_table, MARS_table_ind)
  
  MARS_table <- c(mean(MARS_table$accuracy), mean(MARS_table$discrimination), mean(MARS_table$calibration), mean(MARS_table$precision))
  
}

##Producing GBM table

GBM_table <- c()

for(i in 1:nrow(Taxon_lim)){
  
  GBM_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/GBM_table",sep=""),sep = ""))
  
  GBM_table <- rbind(GBM_table, GBM_table_ind)
  
  GBM_table <- c(mean(GBM_table$accuracy), mean(GBM_table$discrimination), mean(GBM_table$calibration), mean(GBM_table$precision))
  
}

##Producing CTA table

CTA_table <- c()

for(i in 1:nrow(Taxon_lim)){
  
  CTA_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/CTA_table",sep=""),sep = ""))
  
  CTA_table <- rbind(CTA_table, CTA_table_ind)
  
  CTA_table <- c(mean(CTA_table$accuracy), mean(CTA_table$discrimination), mean(CTA_table$calibration), mean(CTA_table$precision))
  
}

##Producing RF table

RF_table <- c()

for(i in c(1:43, 45:nrow(Taxon_lim))){
  
  RF_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/RF_table",sep=""),sep = ""))
  
  RF_table <- rbind(RF_table, RF_table_ind)
  
  RF_table <- c(mean(RF_table$accuracy), mean(RF_table$discrimination), mean(RF_table$calibration), mean(RF_table$precision))
  
}

##Producing MAXENT table

MAXENT_table <- c()

for(i in 1:nrow(Taxon_lim)){
  
  MAXENT_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/MAXENT_table",sep=""),sep = ""))
  
  MAXENT_table <- rbind(MAXENT_table, MAXENT_table_ind)
  
  MAXENT_table <- c(mean(MAXENT_table$accuracy), mean(MAXENT_table$discrimination), mean(MAXENT_table$calibration), mean(MAXENT_table$precision))
  
}

##Producing ANN table

ANN_table <- c()

for(i in 1:nrow(Taxon_lim)){
  
  ANN_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/ANN_table",sep=""),sep = ""))
  
  ANN_table <- rbind(ANN_table, ANN_table_ind)
  
  ANN_table <- c(mean(ANN_table$accuracy), mean(ANN_table$discrimination), mean(ANN_table$calibration), mean(ANN_table$precision))
  
}

##Producing SVM table

SVM_table <- c()

for(i in 1:nrow(Taxon_lim)){
  
  SVM_table_ind <- read.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/",paste(Taxon_lim[i,1],"/SVM_table",sep=""),sep = ""))
  
  SVM_table <- rbind(SVM_table, SVM_table_ind)
  
  SVM_table <- c(mean(SVM_table$accuracy), mean(SVM_table$discrimination), mean(SVM_table$calibration), mean(SVM_table$precision))
  
}

## Create final table

Methods_evaluate <- as.data.frame(rbind(GLM_table, GAM_table, MARS_table, GBM_table, CTA_table, RF_table, MAXENT_table, ANN_table, SVM_table))  

colnames(Methods_evaluate) <- c("accuracy", "discrimination", "calibration", "precision")

Average <- c(mean(Methods_evaluate$accuracy), mean(Methods_evaluate$discrimination), mean(Methods_evaluate$calibration), mean(Methods_evaluate$precision[c(1:7,9)]))

Methods_evaluate <- rbind(Methods_evaluate, Average)

rownames(Methods_evaluate) <- c("GLM", "GAM", "MARS", "GBM", "CTA", "RF", "MAXENT", "ANN", "SVM", "Average")

write.table(file = paste("~/Documents/Patrick_Masa_code/SDMs/figs/","Methods_evaluation",sep=""), Methods_evaluate)

Methods_evaluation <- read.table(file = paste("../figs/","Methods_evaluation",sep=""))

Methods_evaluation[1,1] <- round(Methods_evaluation[1,1], digits = 2)

Methods_evaluation[1:10,1:4] <- round(Methods_evaluation[1:10,1:4], digits = 2)

Methods_evaluation[1:10,1:4]

write.table(file = paste("../figs/","Methods_evaluation",sep=""), Methods_evaluation)



```{r echo=FALSE, results='hold'}
library(knitr)
kable(read.table(file = paste("../figs/","Methods_evaluation",sep="")))
```



##Each species response
Resp <- read_excel("../data/each_species_response.xlsx", range = "A1:C169") 

Resp <- as.data.frame(Resp)

Resp1 <- Resp[1:84, ]

Resp1$X <- as.factor(Resp1$X)

Resp1$Y <- as.factor(Resp1$Y)

name1 <- as.factor(Taxon_lim$species[1:28])

ggplot(Resp1, aes(x=X, y=Y, color=Type, shape=Type)) + 
  geom_point(size=3) +
  scale_x_discrete(breaks=c(1.0,2.0,3.0), labels=c("Negative", "Not Clear", "Positive")) +
  scale_y_discrete(breaks=c(1:28), labels = name1) +
  xlab(NULL) +
  ylab(NULL)

Resp2 <- Resp[85:168, ]

Resp2$X <- as.factor(Resp2$X)

Resp2$Y <- as.factor(Resp2$Y)

name2 <- as.factor(Taxon_lim$species[29:56])

ggplot(Resp2, aes(x=X, y=Y, color=Type, shape=Type)) + 
  geom_point(size=3) +
  scale_x_discrete(breaks=c(1.0,2.0,3.0), labels=c("Negative", "Not Clear", "Positive")) +
  scale_y_discrete(breaks=c(29:56), labels = name2) +
  xlab(NULL) +
  ylab(NULL)


##Each species response2

MARS <- c()

CTA <- c()

RF <- c()

for(i in 45:nrow(Taxon_lim)){

a <- read.table(file = paste("../figs/",paste(Taxon_lim[i,1],"/Scenarios_MARS",sep=""),sep = ""))

b <- read.table(file = paste("../figs/",paste(Taxon_lim[i,1],"/Scenarios_CTA",sep=""),sep = ""))

c <- read.table(file = paste("../figs/",paste(Taxon_lim[i,1],"/Scenarios_RF",sep=""),sep = ""))

MARS <- rbind(MARS, a)

CTA <-rbind(CTA, b)

RF <- rbind(RF, c)

}

p<-c()
for(i in 1:55){
test1<-rbind(MARS[i,2],CTA[i,2],RF[i,2])
test2<-rep(i, times=3)
test3<-c("MARS", "CTA", "RF")
test<-cbind(test1,test2,test3)
p<-rbind(p,test)
}
p<-as.data.frame(p)
colnames(p)<-c("X","Y","Type")
p
p$X<-as.numeric(levels(p$X))[p$X]
name <- as.factor(Taxon_lim$species[1:28])
ggplot(p[1:84,], aes(x=X, y=Y, color=Type, shape=Type)) + 
  geom_point(size=3) +
#  scale_x_discrete(breaks=c(-15,-5,5,15)) +
  scale_y_discrete(breaks=c(1:28), labels = name) +
  xlab(NULL) +
  ylab(NULL)

name <- as.factor(Taxon_lim$species[c(29:43,45:56)])
ggplot(p[85:165,], aes(x=X, y=Y, color=Type, shape=Type)) + 
  geom_point(size=3) +
  #  scale_x_discrete(breaks=c(-15,-5,5,15)) +
  scale_y_discrete(breaks=c(29:55), labels = name) +
  xlab(NULL) +
  ylab(NULL)







