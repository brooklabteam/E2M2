rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(maptools) 
library(rgdal)
library(cartography)
library(raster)
library(colorRamps)
library(rmapshaper)
library(cowplot)


#read in gps data
gps.dat <- read.csv(file = "GPS_dat.csv", header = T, stringsAsFactors = F)
head(gps.dat)

#remove faulty data 
gps.dat$error = 0
gps.dat$error[gps.dat$GroupNum==1 & gps.dat$Date=="2018-09-08" |gps.dat$GroupNum==1 & gps.dat$Date=="2018-09-04"] <- 1
gps.dat$error[gps.dat$GroupNum==1 & gps.dat$Longitude<48.4] <- 1
gps.dat = subset(gps.dat, error !=1)


#separate the different lemurs as a factor
gps.dat$GroupNum = as.factor(gps.dat$GroupNum)

#rewrite date in format known to R
gps.dat$Date = as.Date(gps.dat$Date, format = "%m/%d/%y")
#plot all their points, coloring by ID
p1 <- ggplot(data=gps.dat) + geom_point(aes(x=Latitude, y=Longitude, color=GroupNum))
print(p1)

#add a line to connect dots by site
p2 <- ggplot(data=gps.dat) + geom_point(aes(x=Latitude, y=Longitude, color=GroupNum)) +
      geom_line(aes(x=Latitude, y=Longitude, color=GroupNum, group=Date))
print(p2)


#plot all three lemur movements by day
gps.dat$Date = as.factor(gps.dat$Date)

#now plot
p3 <- ggplot(data=gps.dat) + geom_point(aes(x=Latitude, y=Longitude, color=Date)) + facet_grid(~GroupNum) +
      geom_line(aes(x=Latitude, y=Longitude, group=Date, color=Date)) + theme_bw()
print(p3) 

#pdf file
ggsave(file = "side_by_side.pdf",
       units="mm",  
       width=80, 
       height=50, 
       scale=3, 
       dpi=300)

#as a jpeg
ggsave(file = "side_by_side.jpeg",
       units="mm",  
       width=80, 
       height=50, 
       scale=3, 
       dpi=300)

#could view each individually
lem1 <- subset(gps.dat, GroupNum==1)
lem2 <- subset(gps.dat, GroupNum==2)
lem3 <- subset(gps.dat, GroupNum==3)


plem1 <- ggplot(data=lem1) + geom_point(aes(x=Latitude, y=Longitude, color=Date)) + facet_grid(~GroupNum) +
  geom_line(aes(x=Latitude, y=Longitude, group=Date, color=Date))
print(plem1) 
plem2 <- ggplot(data=lem2) + geom_point(aes(x=Latitude, y=Longitude, color=Date)) + facet_grid(~GroupNum) +
  geom_line(aes(x=Latitude, y=Longitude, group=Date, color=Date))
print(plem2) 
plem3 <- ggplot(data=lem3) + geom_point(aes(x=Latitude, y=Longitude, color=Date)) + facet_grid(~GroupNum) +
  geom_line(aes(x=Latitude, y=Longitude, group=Date, color=Date))
print(plem3)


#use cowplot to show together
comb.plot <- cowplot::plot_grid(plem1,plem2,plem3, nrow=3, ncol=1, align = "v")
print(comb.plot)

ggsave(file = "narrow_range.pdf",
       units="mm",  
       width=80, 
       height=70, 
       scale=3, 
       dpi=300)

#as a jpeg
ggsave(file = "narrow_range.jpeg",
       units="mm",  
       width=80, 
       height=70, 
       scale=3, 
       dpi=300)


#add Madagascar - ideaqlly, you want a shape file for Moramanga and Andasibe

mdg_district<-readOGR("MDG_adm3.shp")
mdg_district<-ms_simplify(mdg_district, keep = 0.15, keep_shapes = T) #the function ms_simplify the shapeile to make it smaller and easier to load

plot(mdg_district)
#add your points to the plot
points(gps.dat$Latitude,gps.dat$Longitude,pch=19,col="red", cex = 5)
# too small to see on all of madages

## change XX and YY to different values 
image(mdg_birth, zlim =c(5,5.5), col = blue2red(10))
plot(mdg_district, add = TRUE)


