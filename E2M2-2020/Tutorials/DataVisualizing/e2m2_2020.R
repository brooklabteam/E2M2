#############################################################################
##   Exploring and visualizing data with R 2016-2020
#############################################################################
##   Ecological and Epidemiological Modeling in Madagascar (E2M2)
##   03-14 of January 2020
##   Institut Pasteur de Madagascar (IPM), Ambohitrakely, Antananrivo, Madagascar
##   ValBio center, Ranomafana, Fianarantsoa, Madaascar
#############################################################################
##   Hafaliana Christian Ranaivoson, 2016-2020
##   Virology Unit, Institut Pasteur de Madagascar
##   Mention Zoology and Animal Biodiversity, Univesity of Antananarivo
#############################################################################
##   The goal of this tutorial is to introduce participants to data 
##   management and data visualisation using the R statistical software
##   To reach this goal, the tutorial will use a simulated dataset created under
##   R software.
##   The data sets was created from real observed parameters from bat capture
##   It also borrows general charateristics of some hemoparasite infection pattern
##   Though simulated, the data point out some aspect of some hemoprasite infection 
##   to help cope how R software can be used to manage raw data and visualize
##   some pattern of biological events.
#############################################################################
##   Data story: 
##   Small fruit bats "fb" were captured from one locality, with thousands of bats.
##   100 bats were net captured every month of the year, for one year.
##   For each captured individual, a unique Id was given "fb_number".
##   The length of the forearm (Forearm, in mm) was measured
##   and the each bat was weighted (Weight, in gramm), sexed (Sex, f/m), and aged (year).
##   Prevalence (Prev, 0/1) of hemoparasite was obtained by observation under microscopy
##   and recorded as 0 when absent and 1 when present.
##   Parasitemia (Para) was recorded as the number of  microscopy positive field
##   out of 200 field observation under high magnification.
##   All data was recorded on excel software as csv file, "e2m2_FB.csv".
#############################################################################
##   Note that the symbol "#" is used to put note within your code. Anything 
##   placed after this symbol in a row will not be execute as code.
#############################################################################
############            Check your Working Directory      ###########
##  First, check you are working in the right directory
##  get the default directory use by R.

getwd()

##### Manually set your own directory (optional)
##  An alternative to creating R project file is to manually set your own working directory .
##  Do not run the next two lines of codes if you have set an R project file.
##


# setwd("D:\\Chris\\E2M2\\Lectures\\Data")

# setwd("D:/Chris/E2M2/Lectures/Data")

##  File path need to be changed to meet your own file path.
##  Note that both the next two line codes work the same way
## \\ or / are used to set the to your files, do not forget the quote "".

######################################################################
############             Impoting data              ##################
######################################################################
##  Import the dataset named "e2m2_FB.csv" to a new data name.

e2m2 <- read.csv("e2m2_FB.csv",stringsAsFactors = F,header = T)

## In the code above, we tell R not to change strings to a factor by setting stringsAsFactor
## to False "F". If we did not set this, R will change all strings entry to factor. 
## Your data should be listed in the R data enviroment list under R studio.

## You can check with R help details about functions.
## Just add an exlamation mark "?" infront.

?read.csv


## Now, check your imported data with the View() function. You can use this fucntion
## to view any data in the R data environment. Another way is just typing the data name
## directly, but it will be shown in the R consol an won't be as nice as View().

View(e2m2)

## and try

e2m2

#################################################################
################        Exploring data      #####################
#################################################################
## Use the dim() and names() function to get an overview of you data

dim(e2m2)

names(e2m2)

## How many columns and rows are there in the data?
## What variables are we interested in?
## Refer to the data story to have more information.

## Now that you know the variable names, try to print out the content of some of them.
## From the lecture, we have seen the data structure and how to access its content.
## Remember the combination of dataset name and variable name by $ symbol to do so.

e2m2$Forearm

## Now filter the content with [...] to show specific data, for example the Id of fruit bats
## which have a forearm length smaller than 56 mm.

e2m2$Id[e2m2$Forearm < 56]


## You can get the count of bat which have a forearm smaller than 56 mm by using
## length(...) function. Which give 130 bats.

length(e2m2$Id[e2m2$Forearm < 56])

## You are told to print out bat Id that have forearm smaller than 40 mm and their number.
## Try it
## What happened?
## No bats have a forearm shorter than 40mm
## Now use range() function to check it.

range(e2m2$Forearm)

## The range() function is usefull when you want to know the two limits of your variable range.
## Instead, you can also use min(), or max() to get it separately. Try it.
## Note that these fuction only work for numerical and date variables.

## You are asked to tell how many bats have a forearm longer than 70 mm
## how many?

## You can add more logic filter within [...] by using "and" or/and "or" logical conjunction.
## "and" is represented by the "&" symbol. It means both filter is true
## "or" is represented by the "|" symbol. It means one or the other is true
## For example, the next code print out the count of bats which have a forearm longer than
## 70 mm and smaller than 75 mm. which is 238 bats

length(e2m2$Id[e2m2$Forearm > 70 & e2m2$Forearm < 75])

## Now try to get the number of bats that have a forearm length within 60 to 70 mm.

########################     Data structure     ########################
## You are now able to access any data within the dataset. It is time to see what type
## of variable do we have. This is important since the operation we can use on them
## will also depend on their type.

## The function str(...) which means structure can be used to do so.

str(e2m2)


## What can you tell after executing the code above?
## Since we set stringAsFactors option to false, when we imported data, we have only two types
## of variables, numeric or charactere.
## The advantage from this is that it will easier to manipulate the data.
## We are going to set by ourselves the variable type we need and check errors.

########### Knowing our data structure and variable types  ##############
## The way we record data on the field may depends on the goal we want to acheive and
## the way we get them. For example, some people may record age data as age group 
## and use character instead of numbers.
## Sometimes, after entering a thousands records, we may do some typing error or
## we just mixed type of entry.

## Try to identify error in the result from the str(...) function we have applied to 
## our data.
##
## When we told R to just read our data as is (see import data above) it just imported
## the data as numeric or as character. Though numeric can be represented with integer
## (int) for non decimal numbers or numeric if decimal numbers.
##
## Most of the case str(..) function will not suffice to figure out erros in the data, just
## because R is only showing the first five rows. You need to set the variables to the
## needed format first and search for errors.
## 
## You will decide what type you want R to treat your variable. The R software work mainly
## on the following type of variables.
##
## - Factor: is a grouped identity. R group all same values in your variable as categories
##   ex: Sex will be grouped as male or female. as.factor(...) is the funtcion 
##   to factor variables. Levels are the categories.
##
## - Numeric: is numeric variable. as.numeric(...) is the function.
##
## - Date: is time variable. as.Date(...) is the function. For this one, you need to tell R
##   how you recorded your date. It means the format.
##   The default format R is reading Date is "%Y-%m-%d", Year-month-day (ex: "1900-12-30")
##   If the way you record your date is not the same, you need to point it out to R.
##   Here, our date is recorded as "%m/%d/%Y", month/day/Year.
##
##   Then we write: as.Date(e2m2$Date,format="%m/%d/%Y)

## - NA: Not Available. They are missing values. There are function to handle Na values.
##
## Let's set our Variable types:
## Do not assign variable, we do not want R to save them in Rdata yet, just try.
##
## In the folowing code, we are trying to factor the Sex Variable:

as.factor(e2m2$Sex)

## What do we get?
## R printed out the values, and at the bottom, the levels. Levels are the number of category in the
## Variable. Levels can be shown with levels() function. Try:

levels(as.factor(e2m2$Sex))

## Yes, you can even use combination of function in R!
## Here we combined levels() and as.factor(), nested functions.
## How many levels did we get? Does it seem right?
## Of course not, Sex should have only two categories, Male or Female. What was wrong?
## You can see, we typed f in three different ways. f F f. Male "m" seems to be OK.
## R is case sensitive, you need to type exactly the same way.
##
## Let's move to another variables:
## We will deal with Age variable. We have already seen something was wrong with it from the str() function we've used.
## We are going to set age variable as numeric.

as.numeric(e2m2$Age)

## R printed out values, and at the bottom, a warining message. It says NA value are introduced.
## It means that R does not know how to change some values to numeric.

## Let's move on to another variable for now, the Date variable:

as.Date(e2m2$Date,format="%m/%d/%Y")

## R ptinted out the values and no warning. We are fine with the date variable then. Note that R
## is showing the Date as it's default format "%Y-%m-%d", this is a good sign.

## Now, we know which variable has error (Sex and age variables), it is time to correct them. Let's move to..

###################           Correcting Data    #######################

########            Finding the wrong values     #############################"
## Now, we know that Age and Sex variables have values errors. We need to find them.
## A simple way is to go back to the excel data. But...
## The problem is that excel is not realy case sensitive
## If you filter the variable Sex in excel, all three f F f will be the same letter.
## We can use R to search for them. Let use structure function to a factored Sex variable.
## You can also try to summary().

str(as.factor(e2m2$Sex))


levels(as.factor(e2m2$Sex))


## R printed out 4 levels and the levels values. What do you see? : "f", "F", "f "
## What different with the three f? The one is Uppercase and the other one has space aftet it.
## Look closely at the quote "". We found them!

## Let's move to the Age variables.
## We need to find which values are not numerical. Before, we see that R return some NA value when we
## tried to set Age as numeric. Then we need to find which produce NA values.

e2m2$Age[is.na(as.numeric(e2m2$Age))]

## Here we accessed the Age variable. Then we filter the content by [...] (do you remember?).
## inside the filter, we want content that produce NA values (with is.na(...)) when we set Age variable
## as numeric.
## R will print out the values we were searching for. What are those values?
## R will surely produce warning message at the bottom but we only need the values.
## As you can see, some age data were recorded as character: "seven", "nine", "one"
## 
## We can use the same logic to find the wrong value for Sex.
## We know that levels should only be male or female, and the right way to type them is "f" and "m".
## All lowercase and no space.
## Try the following code. (!= means not equal to)

e2m2$Sex[e2m2$Sex != "f" & e2m2$Sex !="m"]

# Get the number

length(e2m2$Sex[e2m2$Sex != "f" & e2m2$Sex !="m"])


## Now we found all wrong value, it is time to correct them
##
## ###############         Changing Values         #####################
## Keep in mind that we did not assign any data type to our variables yet.
## They still character or numeric

str(e2m2)

## Changing values of a variable is a straight forward manipulation. It just all about assigning
## the new value to the previous one.
## So, let's assign the right value "f" (lower case, no space) to the wrong value "F" and "f ".

e2m2$Sex[e2m2$Sex=="F"] <- "f"

e2m2$Sex[e2m2$Sex=="f "] <- "f"

## Now try to factor the Sex variable again and pirnt out its structure

as.factor(e2m2$Sex)

str(as.factor(e2m2$Sex))

levels(as.factor(e2m2$Sex))

## How many levels do we have now?         Which are they?         Are we good now?
## Now try to do the same thing to correct value of Age variable, and try as.numeric(...) again.

e2m2$Age[e2m2$Age=="one"] <- "1"

e2m2$Age[e2m2$Age=="two"] <- "2"

e2m2$Age[e2m2$Age=="four"] <- "4"

as.numeric(e2m2$Age)

## Is there any warning message again?
## Now everithing is ok, we need to assign the right format to our variables

e2m2$Sex <- as.factor(e2m2$Sex)

e2m2$Prev <- as.factor(e2m2$Prev)

e2m2$Age <- as.numeric(e2m2$Age)

e2m2$Date <- as.Date(e2m2$Date, format="%m/%d/%Y")

## Now try the str(...) function on our data again.

str(e2m2)

## We are good!

## What we have done so far is to correct already existing variables
## We can also use what we have seen to create new variables!
## It is just about assigning values to new variable name!...What?!
##
## In the example below, we are creating a new categorical variable named "Size"
## with three levels (small, medium, and large) from the Forearm variable.

e2m2$Size[e2m2$Forearm < 60] <- "Small"

e2m2$Size[e2m2$Forearm >= 60 & e2m2$Forearm < 75] <- "Medium"

e2m2$Size[e2m2$Forearm >= 75] <- "Large"

e2m2$Size <- as.factor(e2m2$Size)


View(e2m2)

##########################  Saving RData   #########################
### Now that we have cleaned the raw data, we can save it as RData and/or Write to a new csv file

save.image("Cleaned-e2m2.Rdata")

# This code write all data created or imported to RData Environment
# It will avoid you to do re-do all the steps
# For next use, just load the RData with load("Cleaned-e2m2.Rdata")

write.csv(e2m2,"Cleaned-e2m2.csv",row.names=F)

# This write specific data from the RData environment list to a csv file
# For next use, just import it like we did before

# Now check your R working folder to see these created file
# Now you can clean RData environment without loosing all steps done.
# So lets us try
# clean RData environment first. We will re-import it for the next steps

rm(list=ls())


## Now create a new variable named AgeClass from the Age variable, and factor it.
## for this, we need six categories as follow:
##     -"neonate" younger than or equal to 0.4 year
##     -"juvenile" between 0.4 to 1 year
##     -"subadt" between 1 to 2 year
##     -"adt1" between 2 to 4 year
##     -"adt2" older than 4 year

## Remember to use str(), summary(), range() and other useful function to check your variable.
## These function are very handy and already give useful hint about your data.

## We are done here!
## It is time to play with the data now


############################################################################################
#################################         Visualizing Data        ##########################
############################################################################################

## First, re-import RData, if have cleaned it

load("Cleaned-e2m2.Rdata")

# Look at your RData environment, and check it

View(e2m2)
str(e2m2)

# Youare good for the next steps

#####################     Arrange data   #####################

######## Base function, tapply().

### In the following, we will use these two base functions to get a monthly Parasite load.
# First Create a parasite load, acoording to data story, it was recorded as Positive fiels from 200 total fields.
# Calculate the parasite load as percentage of positive field ((Para/200)*100). Read the data story to understand.
# Create the new variable, named ParLoad


e2m2$ParLoad <- (e2m2$Para/200) 

# Check your database for the new variable.

str(e2m2$ParLoad)
range(e2m2$ParLoad)
summary(e2m2$ParLoad)

# Use tapply to calculate mean parasite load for each month, and Store it in a new data (PLoad)

?tapply

PLoad <- tapply(e2m2$ParLoad,factor(format(e2m2$Date,"%m")),mean)

# Look at the new data

View(PLoad)

# Try to understand what this function did.
# It simply calculate the mean of Parasite load (ParLoad) for each month. And we store it in PLoad (PLoad <- )
# tapply() is very usefull when you want to "summary value" in each catégorie of a catégorical variable.
# The value to be calculated needs to be numeric, and the grouping variable needs to be catégorical (factor)



########  Using packages
## For this we will need the package "dplyr"
## If you do not have this packages yet, install it with the function "install.packages(...)"

# install.packages("dplyr")

## Then load it with require() or library()

require(dplyr)

## The package dplyr contains usefull fonctions for data management.
## We will see three funtion from it: filter(), group_by(), summarise()

#########  group_by(...) 
## This function will produce new dataset  grouped by the value of a variable.Then, operations done on
## the data will be executed by group. It is very usefull when combined with the summarise() function.

#########  summarize(...)
## This function will produce a new dataset which return summary of the previous variables to new variables.
## These new variables are defined by the user.
## 
## Example: You want to summarize the mean of forearm by the sex of bats from our e2m2 dataset.

## First we group our data by Sex, and assign it to a new dataset name 
## (this avoid to overwrite our data source)
groupsex <- group_by(e2m2,Sex)

## Then we use the newly created data to get the mean of forearm and weight, we assign it again to a new dataset
stat_groupsex <- summarise(groupsex,
                           mForearm = mean(Forearm,na.rm=T),
                           mWeight = mean(Weight,na.rm=T),
                           count = n())

## to undestand What we did, try to view the stat_groupsex dataset

View(stat_groupsex)

## As you can see the operation was made by group, The sex variable.
## We created a variable "mean" by using Forearm mean " -> mean = mean(Forearm,na.rm=T)
## na.rm=T tell to remove NA values from the operation mean().
## We also created an new variable "count" that summarize the number of data. count=n()

## Try to do the same thing. summarize the mean of Weight by the Sex of bats.
## And later, try to summarize mean Forearm by Sex and Size of bats.
##
## Note that group_by() only work for factor variable.


## filter(...) will chose a part of the data by row according to a filter you have set.

## For example, you may need to work only with female bats from our data. The you can write.

femalebats <- filter(e2m2, Sex=="f")

## Check the data and scroll down. Only female data is shown

View(e2m2)

## R will create a new dataset which will contain only female bats data. You can add more filtering
## But now, now try to check the variable Sex structure

str(e2m2$Sex)

## R is still keeping the two leves male "m" and female "f"
## This may lead to some problem when sumarizing or plotting your data
## To avoid it, just add droplevels.data.frame() function to the filter() funcion
## Try the following code

femalebats <- droplevels.data.frame(filter(e2m2, Sex=="f"))

## And check the Sex variable structure

str(femalebats$Sex)

## Did you noticed? Only the female level is left

## Try create a data with all male bats.
## You can even add more filter logic.
## Try to create a data with all female with medium size
## how would you do that?


## These three function is very usefull when you want to arrange you data for a specific need, as report.
## Especially when you wannt a graphical view.
## Try to remember these function as they will be used in the following step
## Now it is time to have a graphical Visualization of our data

## But first we need to delete the the data we have created, except the e2m2.
## To do so, R use the rm() function.
## Look at R data gloval environment window before running the code.

rm(groupsex,stat_groupsex,PLoad,femalebats)

## Now look back to your R data environment. If you want to delete all data in
## R environment, then use rm(list=ls()).
## ls() function list all data in the R environment.
## Try

rm(list=ls())

######################################################################
#########################     Data Plotting    #######################
## There are some powerfull packages for graphic in R, like ggplo2. But you can already
## make nice graphics with built in function in R, which is the plot() function.
##

## Before we get int plot(), we will start with "ready to use" graphical function like:
## hist(), boxplot(), pie(), barplot().

# First reload the cleaned data again

load("Cleaned-e2m2.Rdata")

##########    Histogramm and Boxplot functions produce nice graphic.
# They give information about variable distribution.

### Histogramm hist().
# Histogram of Parasitemia

hist(e2m2$Para)

# check function details

?hist

# Switch to probability

hist(e2m2$Para,probability = T)

## Polished graph
hist(e2m2$Para,probability = T,
     main="Parasite load distribution",
     xlab="Parasite load")

######  Boxplot boxplot().

boxplot(e2m2$Forearm)

# Check out details about capability of boxplot functio.

?boxplot

# Now try to regroup the data to get get more interesting information
# For example by sex, use the "~" to do so.

boxplot(e2m2$Forearm~e2m2$Sex)

## A trick: instead of repeating the database name, attach it with attach() function

attach(e2m2)

# Now the database is loaded into the memory so you only have to type variable names
# Important! you have to detach it when you are done or swithing to different database.
# To detach dataset, use dettach() function.
# Try the same thing as above without database name (We have already attached it)

boxplot(Forearm~Sex)

## Add more features to the graph, so it will look nicer.
# Add the full names to the Sex group.
boxplot(Forearm~Sex,
        names=c("Females","Males"))

# Add color to the graph. The color is represented by number.
# Here we have two groups on x axis. Then we used two different number 2 and 8 (c(2,8)).
# Try to play with different color by changing the col values.

boxplot(Forearm~Sex,
        names=c("Females","Males"),
        col=c(2,8))

# Add main title and axis title, and set y limit to 100.

boxplot(Forearm~Sex,
        names=c("Females","Males"),
        col=c(2:3),
        main="Forearm by Sex",
        xlab="Sex",
        ylab="Forearm (mm)",
        ylim=c(0,100))

# With boxplot, you can add Interraction with "*"

boxplot(Forearm~Sex*Prev,
        col=c(2:3),
        main="Forearm by Sex",
        xlab="Sex",
        ylab="Forearm (mm)",
        ylim=c(0,100))

# Now you have to set the X axis value to look nicer.



## We are going to manipulate the e2m2 database for the next graphics
## We then need to detach e2m2 database first

detach(e2m2)

##

#### Barplot is usufull to present single value on each grouping variable.
# For example, the Monthly mean parasite load.
# First, calculate the monthly mean parasite load, like we did before.


PLoad <- tapply(e2m2$Para,factor(format(e2m2$Date,"%m")),mean)

# Then just plot it with barplot()

barplot(PLoad)

# Now you just have to polish the graph, use ?barplot to see it capability.


#####   You cal also use pie()

pie(PLoad)

####### Using dplyr and bar plot for monthly prévalence

# First, do a grouping by month

gp_pload <- group_by(e2m2,Months = factor(format(Date,"%m")))

# Calculate the proportion

st_pload <- summarise(gp_pload,
                      infected = length(Prev[Prev == 1]) / length(Prev[Prev == 0]),
                      Count=n())

# View the data

View(st_pload)

barplot(st_pload$infected,names.arg = st_pload$Months)

### This is it, now let us see plot() function
## First, remove temporary created data

rm(gp_pload,st_pload,PLoad)


###############    R Plot function plot()

attach(e2m2)


#### Check plot() capabilities
# Plotting Weight vs Forearm length

?plot()

# Then plot

plot(Weight~Forearm)

# Let polish the graph and separate male and female

plot(Weight~Forearm,
     main = "Weight/Forearm",
     ylab ="Weight (g)",
     xlab = "Forearm (mm)",
     col=Sex)


# Nicer, isn't it?

# plot() function will chose default parameter according to your variable type
# Here it uses point since your variables have numerous datapoint and it is not grouped.
# Let us try another set of variable

plot(Weight~Sex)

# Do you see, plot() uses boxplot since you are using a categorical variable against numerical.
# It automaticaly groups your graph
# Try different set of variables

## Let us get back to the Forearm vs age

plot(Forearm~Age,
          main = "Forearm/Age",
          ylab ="Forearm (mm)",
          xlab = "Age (year)",
          col=Sex,
     cex=0.5)
# You change the point type and size with "pch" and "cex" options

plot(Forearm~Age,
     main = "Weight/Forearm",
     ylab ="Weight (g)",
     xlab = "Forearm (mm)",
     col=Sex,
     pch=5,
     cex=.2)

# Try to play with option available in plot() function
# Now use plot() to show Weight vs Age

##### Let us play again
# plot() is the base of plotting, you can add another graph to it
# here we are going to add a fitting line to our plot

# First get fit a curve to our data, you will learn this later with next lectures
# For now we will just use a lowess fitting, a non liner fitting.

fitcurv <- lowess(Forearm~Age)

# then plot our base data

plot(Forearm~Age,
     main="Forearm/Age, lowess fitting",
     xlab="Age (year)",
     ylab="Forearm (mm)",
     type="p",
     pch=3,
     cex=0.7)

# Now add a line that fit the data, with red color

lines(fitcurv,col="red",cex=6)

# The line does not realy fit our data, why?
# Of course, we did not separate males and females
# Lets us do it
# First, detach the e2m2 data

detach(e2m2)

# Separate female and male data from the e2m2 database, we already know how to do it

datmal <- filter(e2m2,Sex=="m")
datfem <- filter(e2m2,Sex=="f")

# Then fit a lowess curve on each data

fitfem <- lowess(datfem$Forearm~datfem$Age)
fitmal <- lowess(datmal$Forearm~datmal$Age)
     
# then plot our base data, and separate female and male by color

plot(e2m2$Forearm~e2m2$Age,
     main="Forearm/Age, lowess fitting",
     xlab="Age (year)",
     ylab="Forearm (mm)",
     type="p",
     pch=3,
     cex=0.7,
     col=e2m2$Sex)

## Now add the two lines with different color
lines(fitfem,col="black",lwd=3)
lines(fitmal,col="red",lwd=3)

## Let us add more stuff show plot() capabilities

legend(x=0, y=95,
       legend=c("Males","Females"),
       col=c("red","black"),
       title="Lowess fit",
       lty=1,
       x.intersp = .5,
       y.intersp = .8)

## Nicer! isn't it?
### Dettach all data if not yet done


################ This ends data visualizing with integrated graphic function in R
######  Try to play with it and explore all possible option for each functions
######  R software is so versatile that possibilities are large. It will be hard to list them all
######  The only way you can master it is to play with it. Its is up to you to deepen you knowledge!

######  In the following is an example of graphical package for R, the ggplot2

##################################################################################
######  Quick Introduction to ggplot2, a powerfull graphical library for R.#############

## Clean your RData environment first

rm(list=ls())

# load the library if not done yet

require(ggplot2)

# Check its capabilities

?ggplot

# load the e2m2 Rdata if not done yet

load("Cleaned-e2m2.Rdata")

## Let's built our base plot with Age in y axis and Forearm, in x axis

ggplot(e2m2,aes(Age,Forearm))  ##Try this


## ggplot created a blank plot with our xy axis as Foearm and Weight.
## We now need to tel R to populate the plot with the function geom_
## Tehre is a list of geom in ggplot. It is up to you to choose which best fit for your data.
## We will start with point since we have two continuous variables.

ggplot(e2m2,aes(Age,Forearm)) +
  geom_point()

## try Line

ggplot(e2m2,aes(Age,Forearm)) +
  geom_line()

## put them in the same plot

ggplot(e2m2,aes(Age,Forearm)) +
  geom_point() +
  geom_line()

## Line does not realy make sens since it is just linking every data points
## Using geom_smooth() is nicer and more informative as it is giving a smoothing line.

## Add a smooth to the point plot
ggplot(e2m2,aes(Age,Forearm)) +
  geom_point() +
  geom_smooth()


############  Mapping aesthetic vs Fixed Value   ##################
ggplot(e2m2,aes(Age,Forearm)) +
  geom_point()

## The code above plot the lenght of bats Weight by the Forearm. You can clean the plot by grouping
## your variable and give clearer separation in the curve.
## ggplot offers you a list of aesthetics that can be used to polish your plot.
## 
## First, let's give a fixed value to our aesthetics

ggplot(e2m2,aes(Age,Forearm)) +
  geom_point(color="red", shape=3)

## Not better
## You can map your variable to the aesthetic to give a clear separation.

ggplot(e2m2,aes(Age,Forearm)) +
  geom_point(aes(color=Sex, shape=Sex))

## Note that, When you want to map the aesthetic to a variable, you need to use the aes().
## All aesthetic outside that function will be considered as fixed.

##################       Polishing the plot       ###############
## You can access plot element and set a custom setting to give nice look to your plot
## There are a big list of theming in ggplot. But We will just see "ggtitle() and scale()
## 

ggplot(e2m2,aes(Age,Forearm)) +
  geom_point(aes(color=Sex, shape=Sex)) +
  
  ggtitle("Forearm by Age with ggplot2")+
  
  scale_x_continuous(name="Age (Year)",
                     limits=c(0,5))+
  
  scale_y_continuous(name="Forearm (mm)",
                     limits=c(40,100))+
  
  scale_color_discrete(name="Sex",
                       breaks=c("f","m"),
                       label=c("female","male"))+
  
  scale_shape_discrete(name="Sex",
                       breaks=c("f","m"),
                       label=c("female","male")) +
  geom_smooth(aes(color=Sex),method = "loess") +
  theme_classic()



## Alway keep in mind that this tutorial is a way like other to handle data. You should try playing with data
## in R. It is the only way you can master it.

##############       Thanks A lot !!!!  ##################


