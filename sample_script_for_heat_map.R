# load packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ContourFunctions)
library (lattice)
#read data and specify NAN as NA

rawtimeseries <- read.csv (file="Saanich_TimeSeries_Chemical_DATA.csv", header=TRUE, stringsAsFactors=FALSE, na.strings=c("NA", "NAN"))
rawtimeseries<-mutate(rawtimeseries, Date=as.numeric(raw_dat$Date)) 
#making NO3vsdepth


NO3dat<-
  select(rawtimeseries, "Cruise", "Date", "Depth","NO3")

NO3dat<- na_if(NO3dat, "NAN")
NO3dat<- na_if(NO3dat, "ND")
NO3dat<- filter(NO3dat, !is.na(NO3))

x <- NO3dat$Cruise
y <- NO3dat$Depth
z <- as.numeric(as.character(NO3dat$NO3))

x.scale <- list(at=seq(0,100,12))

#lattice plot 
df = data.frame(x,y,z)
levelplot(z ~ x*y, data=df  , xlab="Months",ylab="Depth" ,ylim=c(200,0), col.regions = heat.colors(100)[length(heat.colors(100)):1]   , main="",interpolate=TRUE,scales=list(x=x.scale,y=y))

#just NO3 plot
ggplot() +   geom_point(aes(x=x, y=y, color=z)) +    scale_y_reverse(limits=c(200, 0)) 


