# load packages
library(tidyverse)
library(ggplot2)
library(dplyr)
library(ContourFunctions)
library (lattice)
#read data and specify NAN as NA

rawtimeseries <- read.csv (file="Saanich_TimeSeries_Chemical_DATA.csv", header=TRUE, stringsAsFactors=FALSE, na.strings=c("NA", "NAN"))
#making NO3vsdepth


NO3dat<-
  select(rawtimeseries, "Cruise", "Date", "Depth","NO3")



NO3dat<- na_if(NO3dat, "NAN")
NO3dat<- na_if(NO3dat, "ND")
NO3dat<- filter(NO3dat, !is.na(NO3))
x <- NO3dat$Date
x <- NO3dat$Cruise
y <- NO3dat$Depth
z <- as.numeric(as.character(NO3dat$NO3))

#lattice plot 
df = data.frame(x,y,z)
levelplot(z ~ x*y, data=df  , xlab="X" , col.regions = heat.colors(100)[length(heat.colors(100)):1]   , main="",interpolate=TRUE)

#just NO3 plot
ggplot() +   geom_point(aes(x=x, y=y, color=z)) +    scale_y_reverse(limits=c(200, 0))

#selecting PO4Data
PO4dat<-
    select(filter(rawtimeseries, !is.na(PO4)), "Cruise", "Date", "Depth","PO4")

x_PO4 <- PO4dat$Date
y_PO4 <- PO4dat$Depth
z_PO4 <- as.numeric(as.character(PO4dat$PO4))


# PO4 lattice plot 
df_PO4= data.frame(x_PO4, y_PO4, z_PO4)
levelplot(z_PO4 ~ x_PO4*y_PO4  , xlab="X" , col.regions = heat.colors(100)[length(heat.colors(100)):1]   , main="",interpolate=TRUE)

#gg PO4 plot
ggplot() +   geom_point(aes(x=x_PO4, y=y_PO4, color=z_PO4)) +    scale_y_reverse(limits=c(200, 0))
