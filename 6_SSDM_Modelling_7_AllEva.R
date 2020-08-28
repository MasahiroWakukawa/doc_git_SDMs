require(dismo)
require(raster)
require(magrittr)
require(foreach)
#require(doParallel)
require(ff)
require(SSDM)
require(dplyr)
require(ggplot2)
require(rgeos)
require(maptools)
library(knitr)
library(kableExtra)
library(rJava)
library(maps)


Australia_Scenarios <- read.table(file = paste("../figs/","Stacked_analysis/Australia_Evaluation",sep=""))

Africa_Scenarios <- read.table(file = paste("../figs/","Stacked_analysis/Africa_Evaluation",sep=""))

Asia_Scenarios <- read.table(file = paste("../figs/","Stacked_analysis/Asia_Evaluation",sep=""))

Europe_Scenarios <- read.table(file = paste("../figs/","Stacked_analysis/Europe_Evaluation",sep=""))

SouthAmerica_Scenarios <- read.table(file = paste("../figs/","Stacked_analysis/SouthAmerica_Evaluation",sep=""))

NorthAmerica_Scenarios <- read.table(file = paste("../figs/","Stacked_analysis/NorthAmerica_Evaluation",sep=""))

Scenarios <- rbind(Australia_Scenarios, Africa_Scenarios, Asia_Scenarios, Europe_Scenarios, SouthAmerica_Scenarios, NorthAmerica_Scenarios)

rownames(Scenarios) <- c("Australia", "Africa", "Asia", "Europe", "SouthAmerica", "NorthAmerica")

colnames(Scenarios) <- c("log(RCP2.6/Presence)", "log(RCP8.5/Presence)")

Scenarios[1:6,1:2] <- round(Scenarios[1:6,1:2], digits = 3)

write.table(file = paste("../figs/","Stacked_analysis/AllContinents_Evaluation",sep=""), Scenarios)


```{r echo=FALSE, results='hold'}
library(knitr)
kable(read.table(file = paste("../figs/","Stacked_analysis/AllContinents_Evaluation",sep="")))




load(file = "../output/RData_Australia_Presence_all.RData")

load(file = "../output/RData_Australia_RCP26_all.RData")

load(file = "../output/RData_Australia_RCP85_all.RData")

Australia_difference <- mosaic(Australia_Presence, Australia_RCP85, fun=diff)

plot(Australia_difference, xlab="longitude", ylab="latitude", main="Australia")


load(file = "../output/RData_Africa_Presence.RData")

load(file = "../output/RData_Africa_RCP26.RData")

load(file = "../output/RData_Africa_RCP85.RData")

Africa_difference <- mosaic(SSDM_Africa_Presence@diversity.map, SSDM_Africa_RCP85@diversity.map, fun=diff)

plot(Africa_difference, xlab="longitude", ylab="latitude", main="Africa")


load(file = "../output/RData_Asia_Presence_all.RData")

load(file = "../output/RData_Asia_RCP26_all.RData")

load(file = "../output/RData_Asia_RCP85_all.RData")

Asia_difference <- mosaic(Asia_Presence, Asia_RCP85, fun=diff)

plot(Asia_difference, xlab="longitude", ylab="latitude", main="Asia")


load(file = "../output/RData_Europe_Presence_all.RData")

load(file = "../output/RData_Europe_RCP26_all.RData")

load(file = "../output/RData_Europe_RCP85_all.RData")

Europe_difference <- mosaic(Europe_Presence, Europe_RCP85, fun=diff)

plot(Europe_difference, xlab="longitude", ylab="latitude", main="Europe")


load(file = "../output/RData_SouthAmerica_Presence.RData")

load(file = "../output/RData_SouthAmerica_RCP26.RData")

load(file = "../output/RData_SouthAmerica_RCP85.RData")

SouthAmerica_difference <- mosaic(SSDM_SouthAmerica_Presence@diversity.map, SSDM_SouthAmerica_RCP85@diversity.map, fun=diff)

plot(SouthAmerica_difference, xlab="longitude", ylab="latitude", main="SouthAmerica")


load(file = "../output/RData_NorthAmerica_Presence_all.RData")

load(file = "../output/RData_NorthAmerica_RCP26_all.RData")

load(file = "../output/RData_NorthAmerica_RCP85_all.RData")

NorthAmerica_difference <- mosaic(NorthAmerica_Presence, NorthAmerica_RCP85, fun=diff)

plot(NorthAmerica_difference, xlab="longitude", ylab="latitude", main="NorthAmerica")







```      