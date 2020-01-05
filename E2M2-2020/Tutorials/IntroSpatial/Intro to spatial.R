
######## INTRODUCTION TO SPATIAL ANALYSIS IN R
####this is a tutorial that will help you make your first maps and run spatial analysis using R statistical software

#load required packages that we will use for this exercise
library(maptools) 
library(rgdal)
library(cartography)
library(raster)
library(colorRamps)
library(rmapshaper)

###############   Reading and mapping from shapefile data
######
### we will first read a shapefile representing the different fokontany of Antananarivo renivohitra and the different markets of the city
getwd()

fokontany<-readOGR("cua_fkt_1.shp") #the function readOGR from the package rgdal allows you to import vector maps into R,
#remember the shapefile (.shp) must be in the same folder as the other
plot(fokontany, main="map of the fokontany of Antananarivo")
fokontany
head(fokontany@data)

tsena<-read.csv("marketsGPS.csv")
points(tsena$x,tsena$y,pch=16,col="red")

#what if we want to plot only fokontany in which there are no people registered in 2003?
#just like a dataframe, the data embedded in the shapefile can be indexed
nopeople<-fokontany[(fokontany$Pop_2003 ==0),]
plot(fokontany,main="fokontany of Antananarivo")
plot(nopeople,col="green",add=T)
legend("bottomleft","no people registered 
       in 2003",pch=20, col="green",bty=F)
###add=TRUE adds the most recent plot to the previous

#####Let's make a thematic map : we are interested in the population size of the different fokontany of Antananarivo between 1993-2003
####
library(cartography) # package cartography helps make thematic maps (e.g. proportional symbols or choropleth representation)
par(mfrow=c(1,1)) #par(mfrow=c(y,x)) allows you to put several graphs on the same sheet


plot(fokontany, col = "grey60",border = "white", lwd=0.4,main="Population of Antananarivo 1993")
fokontany@data$Pop_1993<-as.numeric(fokontany@data$Pop_1993)

cols <- carto.pal(pal1 = "multi.pal",n1=20)
choroLayer(spdf = fokontany, 
           df = fokontany@data,
           var = "Pop_1993",
           col = cols, 
           border = "grey40", 
           lwd = 0.5, 
           legend.pos = "bottomleft",
           legend.title.txt = "Population",
           add = T)

plot(fokontany, col = "grey60",border = "white", lwd=0.4,main="Population of Antananarivo 2003")
fokontany@data$Pop_2003<-as.numeric(fokontany@data$Pop_2003)

cols <- carto.pal(pal1 = "multi.pal",n1=20)
choroLayer(spdf = fokontany, 
           df = fokontany@data,
           var = "Pop_2003",
           col = cols, 
           border = "grey40", 
           lwd = 0.5, 
           legend.pos = "bottomleft",
           legend.title.txt = "Population",
           add = T)


## ------------------------------------------------------------------------
# Reading Raster Data

# Next we will practice reading in raster data

library(raster)
library(colorRamps)
mdg_pop<-readOGR("MDG_pop/")
## load the file from raster worldpop of estimated births per raster grid cell from Worldpop
mdg_birth<-raster('MDG_Births/MDG_births_pp_v2_2015.tif')
# After reading in these data, we will then plot these raster images using various color schemes 
# and highlighting particular values. 
par(mfrow=c(1,2))
image(mdg_birth, col = blue2red(10)) #the blue2red function from colorRamps 
image(log(mdg_birth+1), col = blue2red(10))
# Now we will also plot these raster images with the Madagascar shp files
par(mfrow=c(1,1))
library(rmapshaper)
mdg_district<-readOGR("MDG_adm3.shp")
mdg_district<-ms_simplify(mdg_district, keep = 0.15, keep_shapes = T) #the function ms_simplify the shapeile to make it smaller and easier to load


image(mdg_birth, zlim=c(0,5), col = blue2red(10), main="birth rates in Madagascar") #here we use the argument zlim to limit the values of birth rate between xx and yy
plot(mdg_district, add = TRUE)

## change XX and YY to different values 
image(mdg_birth, zlim =c(5,5.5), col = blue2red(10))
plot(mdg_district, add = TRUE)

## Where are places in Madagascar with a birth rate above 5?


#############################################################
#### JOHN SNOW EXAMPLE, THE COST OF LIVING NEAR A PUMP
#########################################################
par(mfrow=c(1,1))
sohoSG<- readGDAL("sohoSG.tif")

names(sohoSG) <- c("snowcost_broad", "snowcost_not_broad")
buildings <- readOGR("buildings.shp", "buildings", integer64="allow.loss")
proj4string(sohoSG) <- CRS(proj4string(buildings))

deaths <- readOGR("deaths.shp", "deaths", integer64="allow.loss")
names(deaths) <- c("cat", "long", "lat", "Num_Cases", "snowcost_broad",
                   "snowcost_not_broad", "b_nearer")
head(deaths@data)
head(sohoSG@data)
####overlay deaths and the streetmap
o <- over(deaths, sohoSG)
deaths <- spCbind(deaths, o) #cbind is just a function to join two spatial datasets
deaths$b_nearer <- deaths$snowcost_broad < deaths$snowcost_not_broad
by(deaths$Num_Cases, deaths$b_nearer, sum)

boxplot(snowcost_broad ~ b_nearer, deaths,main="Where do the deceased people live?", ylim=c(0,450),
        ylab="distance", xlab="Broad Street", col=grey.colors(1, 0.8, 0.8, 2.2)) #this makes a boxplot that shows that
nb_pump <- readOGR("nb_pump.shp", "nb_pump")
b_pump <- readOGR("b_pump.shp", "b_pump")
oopar <- par(mar=c(1,1,1,1)+0.1)
library(RColorBrewer)
image(sohoSG, "snowcost_broad", breaks=seq(0,750,50),
      col=colorRampPalette(brewer.pal(7, "Reds"))(15))
plot(buildings, col="white", add=TRUE)
plot(buildings, angle=45, density=10, col="grey70", add=TRUE)

symbols(coordinates(deaths), circles=4*sqrt(deaths$Num_Cases),
        inches=FALSE, add=TRUE, bg=c("brown2","grey40")[deaths$b_nearer+1]) #this function draws symbols on a plot 

plot(nb_pump, add=TRUE, pch=8, cex=1.3, lwd=2)
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=8, col="white")
plot(b_pump, add=TRUE, pch=4, cex=1.5, lwd=6)

rect(528900, 181250, 529300, 181385, border=NA, col="white")
legend(c(528910, 529100), c(181350, 181380),
       legend=c("Broad Street pump","other pumps"), pch=c(4,8), bty="n",
       cex=0.6, y.inter=0.7)

legend(c(528910, 529100), c(181275, 181325),
       legend=c("nearer Broad Street pump","nearer other pump"),
       fill=c("grey40","brown2"), bty="n", cex=0.6, y.inter=0.7)
box()

