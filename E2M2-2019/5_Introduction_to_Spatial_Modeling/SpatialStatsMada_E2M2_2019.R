### Spatial Modeling, Data, and Statistics in R
## author: Amy Wesolowski, Ben Rice

### First, you need to set your current working directory to the location where the spatial data 
### files are located. This will be different for each person.

getwd()
setwd("/Users/benjaminrice/Dropbox/Lab Projects/3 E2M2 2019/Tasks and Lectures/3 Introduction to Spatial Modeling (Thursday Jan 17 pm)/E2M2_Introduction_to_Spatial_Modeling")

## ------------------------------------------------------------------------
# Reading shp files

# For the first part of the lecture, we will go through different spatial data sets and how you can 
# load them into R. There are multiple packages that will enable you to read in .shp files in R, 
# but here we will use maptools. It is important that all of the corresponding files (.dbf, .prj, .shx) 
# are in the same folder as the .shp file. These other files include additional information about 
# the spatial data frame that are necessarily to read a .shp file. We will use two functions, one 
# to read in a point shp file and one to read in a polygon shp file. If an error message comes up 
# that says use rgdal::readOGR or sf::st_read, these messages are just to indicate that other 
# functions to read shp files may be more efficient or newer. First we will read in and load libraries
# for: 

library(maptools)
library(rgdal)
library(rgeos)

### all admin 2s for Madagascar
mdg_admin2_shp<-readShapePoly('MDG_Shp/MDG_adm2.shp', proj4string = CRS('+proj=longlat'))

### what if we just want to plot a single location
subset_shp_file<-mdg_admin2_shp[which(mdg_admin2_shp$NAME_2 == 'Sava'),]

mdg_admin2_shp<-gSimplify(mdg_admin2_shp, tol = 0.01)

subset_shp_file<-gSimplify(subset_shp_file, tol = 0.01)

### plot the madagascar shp file
plot(mdg_admin2_shp)

### now we can plot these images
plot(subset_shp_file, col = 'blue', add = TRUE)

## ------------------------------------------------------------------------
# Reading Raster Data

# Next we will practice reading in raster data and extracting these data (and loading the libraries)

library(raster)
library(rgdal)
library(colorRamps)

## estimated births per raster grid cell from Worldpop
mdg_birth<-raster('MDG_Births/MDG_births_pp_v2_2015.tif')

# After reading in these data, we will then plot these raster images using various color schemes 
# and highlighting particular values. 
par(mfrow=c(1,3))
image(mdg_birth, col = blue2red(100))
image(log(mdg_birth+1), col = blue2red(10))


## TO DO: How does the image change if you change the raster values to be 
## between 0 and 5, or 5 and above? 

# Now we will also plot these raster images with the Madagascar shp files
par(mfrow=c(1,3))

image(mdg_birth, zlim = c(0,5), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)

## change XX to 5
image(mdg_birth, zlim = c(0,XX), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)
## This should give you the same plot

## change XX to 4
image(mdg_birth, zlim = c(0,XX), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)

## change XX to 3
image(mdg_birth, zlim = c(0,XX), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)

## change XX to 1
image(mdg_birth, zlim = c(0,XX), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)

## change XX to 0.5
image(mdg_birth, zlim = c(0,XX), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)

## change XX to 5 again
image(mdg_birth, zlim = c(0,XX), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)

## change YY to 2, XX to 5
image(mdg_birth, zlim = c(YY,XX), col = blue2red(10))
plot(mdg_admin2_shp, add = TRUE)

## Where are places in Madagascar with a birth rate above 5?


## Does this make sense??





## ------------------------------------------------------------------------
# Reading in point data sets and analyzing the distance between locations

# We will work with the Snow (cholera) data again and calculate some additional summary values. 
# The pumps and deaths datapoints were digitized by Rusty Dodson at the National Center for Geographic
# Information & Analysis (NCGIA) at UC Santa Barbara. First we will read in each data set. 

#library(spatstat) # do not need, but has lots of great functions for spatial point process data

### reading in the different files 
deaths_points_file<-read.csv('Snow_deaths.csv', header = TRUE)

pumps_points_file<-read.csv('Snow_pumps.csv', header = TRUE)

street_points_file<-read.csv('Snow_streets.csv', header = TRUE)


# The death and pump data set are x,y coordinates (not longitude/latitude) that can be read 
# and plotted directly. 

# The street data set is read in as endpoints of lines that will need to be plotted in a different 
# manner. First we will remake the plot that we made during the lecture.

# First, let's make the same plot we made in lecture

par(mfrow=c(1,1))
plot(NA, NA, xlim = range(street_points_file$x), ylim = range(street_points_file$y), xlab = '', ylab = '', bty = 'n', xaxt = 'n', yaxt = 'n')

for(ii in 1:max(street_points_file$street)){
  sub_dat<-street_points_file[which(street_points_file$street == ii),]
  lines(c(sub_dat$x[1], sub_dat$x[2]), c(sub_dat$y[1], sub_dat$y[2]))
}

points(deaths_points_file$x, deaths_points_file$y, col = 'red')

points(pumps_points_file$x, pumps_points_file$y, col = 'blue', pch = 16)

legend('bottomleft', legend = c('pumps', 'deaths'), col = c('blue', 'red'), pch = 16)

# What do we observe? Where are the pumps in relation to the deaths?
# What pumps are suspicious?

## We will also write two additional functions that we will use to plot both the basic streets 
##  as well as the pump, death and street data. 

plot_streets<-function(){
  plot(NA, NA, xlim = range(street_points_file$x), ylim = range(street_points_file$y), xlab = '', ylab = '', bty = 'n', xaxt = 'n', yaxt = 'n', add = TRUE)
  
  for(ii in 1:max(street_points_file$street)){
    sub_dat<-street_points_file[which(street_points_file$street == ii),]
    lines(c(sub_dat$x[1], sub_dat$x[2]), c(sub_dat$y[1], sub_dat$y[2]))
  }
}

# ^Note that the inside of this function is exactly the same as the code we wrote above to produce
# the street map

# Functions make life easier

plot_base_snow_plot<-function(){
  par(mfrow=c(1,1))
  plot(NA, NA, xlim = range(street_points_file$x), ylim = range(street_points_file$y), xlab = '', ylab = '', bty = 'n', xaxt = 'n', yaxt = 'n')
  
  for(ii in 1:max(street_points_file$street)){
    sub_dat<-street_points_file[which(street_points_file$street == ii),]
    lines(c(sub_dat$x[1], sub_dat$x[2]), c(sub_dat$y[1], sub_dat$y[2]))
  }
  
  points(deaths_points_file$x, deaths_points_file$y, col = 'red')
  
  points(pumps_points_file$x, pumps_points_file$y, col = 'blue', pch = 16)
  
  legend('bottomleft', legend = c('pumps', 'deaths'), col = c('blue', 'red'), pch = 16)
}

# ^Note that the inside of this function is the exact same as the code we wrote to plot the deaths
# and pumps above

## First we will calculate the average x, y coordinates of all the deaths
## (Where is the center of the cholera deaths?)
## And then we plot a circle around this point 
## (size of the circle is proportional to the standard deviation)

mean_x = mean(deaths_points_file$x)
mean_y = mean(deaths_points_file$y)

sd = sqrt(sum((deaths_points_file$x - mean_x)^2 + (deaths_points_file$y - mean_y)^2)/length(deaths_points_file$x))

plot_base_snow_plot()
bearing = 1:360*pi/180
cx = mean_x + sd * cos(bearing)
cy = mean_y + sd * sin(bearing)
circle <- cbind(cx, cy)
lines(circle, col = 'green', lwd = 2)

## *How would we make the line red and make it much thicker?
## *What information does this circle potentially give us

## Question: What is the pump that is on average closest to the most deaths?

## In the lecture, we calculated the distribution of distances from two sample pumps to all 
## of the cases. Now we will calculate the distances from all other points and identify if 
## the cases really were closer to the Broad Street Pump relative to all other pumps. We will 
## do this by comparing the mean distance. 

numb_pumps<-nrow(pumps_points_file)
mean_distance<-rep(NA, numb_pumps)
for(jj in 1:numb_pumps){
  single_pump_coord<-pumps_points_file[jj,]
  euc_dist_from_single_pump<-sqrt((deaths_points_file$x - single_pump_coord$x)^2 + (deaths_points_file$y - single_pump_coord$y)^2)
  mean_distance[jj] = mean(euc_dist_from_single_pump)
}

## This function above calculated the mean distance. Now let's look at the results:
mean_distance_w_pumps<-data.frame(pumps_points_file$label, mean_distance)
View(mean_distance_w_pumps)

## We can plot this as well:
plot(mean_distance_w_pumps, las = 2, cex.axis = 0.65, xlab = '', ylab = 'Mean Distance')

## How does the mean distance from Broad Street pump compares to the other pumps? 

## TO DO: Would you make the same inference if you compared the distributions of the distances 
## (not just the means)? How would you compare these?

## Hint: calculate the minimum distance or maximum distance reusing the above for loop 



######################################################################################################################## 
# Analyzing vaccination coverage in Madagascar
######################################################################################################################## 

## Now, we will read in data from the 2008-2009 Madagascar Demographic Health Survey 
## (in the file: Madagascar2008-2009.csv). These large scale, cross-sectional surveys 
## are conducted globally, freely available, and may contain relevant covariates. We will 
## first bring the data:

library(mgcv)
dhs<-read.csv('Madagascar2008-2009.csv', header = TRUE)

## We will focus on the variables 'age', 'measles.y', and the geographic location 
## ('long' and 'lat'). measles.y is coded as a 1 if the mother can report on whether 
## the child was vaccinated, otherwise it is coded as a 0. Some of the coordinates are 
## incorrect, so let's make a new variable identifying which rows in the data set have 
## true longitude and latitude coordinates. 

good<-which(dhs$long != 0 & dhs$lat != 0, arr.ind = TRUE)

## With this data, we can fit a very simple non-linear statistical model to vaccination coverage, 
## which essentially ‘smooths’ across age, and ‘smooths’ across space. 
## This is clearly very simplistic! (All models are ____ but some are _____)
## But is presented here as a starting point from which further analyses could proceed. 
## Here, we use the package mgcv, which fits ‘generalized additive models’ 
## (or gams, see Wood et al. 2015), following a syntax much like the regression syntax in R, 
## but where ‘s’ indicates ‘smoothed’ terms. Type ?gam into your console if you want to learn
## more about this function. We fit this model keeping only the ‘good’ values defined above 
## in the data-set.

fit<-gam(measles.y~s(age.in.months)+s(long,lat), family = 'binomial', data = dhs[good,])

## and we can see if these covariates significantly explain the patterns:
summary(fit)

## The stars tell us that both of these covariates significantly improve model fit. 
## We can therefore plot the projected patterns over age and space inferred by this model, 
## first looking at the predicted pattern over age:
plot(fit, select = 1, trans = function(x)exp(x)/(1+exp(x)), xlim = c(0,60), 
     ylim = c(0,1.5), xlab = 'Age in months', ylab = 'Proportion vaccinated')

## Note that since we are using binomial (0,1) data, we’ve effectively fitted a logistic transform, 
## so the function ‘trans’ is used to bring our results back to the 0,1 scale. 
## The pattern predicted by the model broadly makes sense - 
## most vaccination is delivered in very young kids (aged < 12 months) and then vaccination rates 
## tail off. 

## We can also plot patterns across space:

library(maps) 
vis.gam(fit,view=c("long","lat"),plot.type="contour",too.far=0.1,type="response",color="cm",ylim=c(-30,-10),xlim=c(35,55), xlab="", ylab="", main="") 
points(dhs$long[good], dhs$lat[good], pch=19,cex=0.2)
map(add=TRUE)

## Does this make sense?

## This again, broadly, makes sense - in Madagascar, 
## it is reasonable to expect that the highest coverage (shown here in purple) 
## will be around the capital, Antananarivo, which is towards the center of the country; 
## and lower coverage (blue) in the south. 


## To figure out how many susceptible children this distribution of vaccination coverage 
## will result in (which then is of relevance for the SIR models, as we can project incidence 
## through time, but also more formally defines the population at risk, and lets us known
## where most vulnerable children are to be found) we could bring in layers including the 
## details of population density of children aged <5, or birth rates across space from 
## worldpop.org.uk. See for example Takahashi et al. (2015).

# Issues to consider with Generalized Additive Models
# We use gams here as a descriptive tool, and do not go into the details of issues associated 
# with over-fitting, out-of-data prediction, choice of smooth terms, etc, 
# but these should all be considered for more serious use of these methods!

## TO DO: Do you see the same patterns with the first polio vaccination? 
## Why would you (or wouldn't you) expect the same relationships?

## Hint: 
## change fit<-gam(measles.y~s(age.in.months)+s(long,lat), family = 'binomial', data = dhs[good,]) 
## using polio0.y 


















# First, peak at the dhs data frame to remind us what our column names are:
head()

# Subset:
# In other words, filter out the bad latitude and longitude data:
good<-which(dhs$long != 0 & dhs$lat != 0, arr.ind = TRUE)


# Re-fit the model but use POLIO instead of MEASLES vaccination:

## MEASLES WAS:
# fit<-gam(measles.y~s(age.in.months)+s(long,lat), family = 'binomial', data = dhs[good,])
## POLIO: fit_polio


## Check the summary of the fit:
summary()


## Maybe we want to compare the change in vaccination rates by age in months for POLIO too

## MEASLES WAS:
# fit_measles<-gam(measles.y~s(age.in.months)+s(long,lat), family = 'binomial', data = dhs[good,])
# plot(fit_measles, select = 1, trans = function(x)exp(x)/(1+exp(x)), xlim = c(0,60), 
# ylim = c(0,1.5), xlab = 'Age in months', ylab = 'Proportion vaccinated')
## POLIO:


## To remind ourselves and to compare let's plot MEASLES again:

# fit_measles<-gam(measles.y~s(age.in.months)+s(long,lat), family = 'binomial', data = dhs[good,])
# plot(fit_measles, select = 1, trans = function(x)exp(x)/(1+exp(x)), xlim = c(0,60), 
# ylim = c(0,1.5), xlab = 'Age in months', ylab = 'Proportion vaccinated')



## Now let's plot POLIO on a map (add POLIO to the X axis label to make things easier)

# MEASLES WAS:
# vis.gam(fit,view=c("long","lat"),plot.type="contour",too.far=0.1,type="response",color="cm",ylim=c(-30,-10),xlim=c(35,55), xlab="", ylab="", main="") 
# points(dhs$long[good], dhs$lat[good], pch=19,cex=0.2)
# map(add=TRUE)

# POLIO:



## And compare to MEASLES (add MEASLES to the X axis label to make things easier)

vis.gam(fit_measles,view=c("long","lat"),plot.type="contour",too.far=0.1,type="response",
color="cm",ylim=c(-30,-10),xlim=c(35,55), xlab="MEASLES", ylab="", main="") 
points(dhs$long[good], dhs$lat[good], pch=19,cex=0.2)
map(add=TRUE)


## What do we observe is different between POLIO and MEASLES?






