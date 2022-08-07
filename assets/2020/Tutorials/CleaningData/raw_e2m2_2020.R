#############################################################################
##   Data cleaning with in R, 2016-2020
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

######################################################################
############             Impoting data              ##################
######################################################################
##  Import the dataset named "raw_e2m2_FB.csv" to a new data name.

rawe2m2 <- read.csv("raw_e2m2.csv",stringsAsFactors = F,header = T,sep=";",dec=",")

############  Look at the dataset  ##############

View(rawe2m2)

# or

rawe2m2

##  Data overview 
dim(rawe2m2)

names(rawe2m2)

str(rawe2m2)

head(rawe2m2)



#############  FORMATTING DATA  #################

### Going variable by variable

########  Id variable  #######
## The Id variable reprents individual record, a duplicate in it's values should not be correct.
# Finding unwanted duplicate

length(unique(rawe2m2$Id))

rawe2m2$Id[duplicated(rawe2m2$Id)]

## Now look at your data with View and scroll to the Id with duplicate from the data view windows

View(rawe2m2)

## It is obvious that bat Id number 23 was repeated twice, isn't it? We can see that the two identical value has different Sex.
## The first is a male but the seond one is a female. Also, there is no bat Id number 24
## We can corretc it with that observation. What we want is to give the secod one the Id number 24

rawe2m2$Id[rawe2m2$Id == "bat_23" & rawe2m2$Sex == "F"] <- "bat_24"

## Now Check back
length(unique(rawe2m2$Id))

rawe2m2$Id[duplicated(rawe2m2$Id)]

##########  Site variable   ########
### If you check at site variable, how many site should it be in it?
### If we rely one the number, it goes to 6 Site_1 to Site_9
##k Let check how many are there

length(unique(rawe2m2$Site))

## Now let us see the levels

levels(factor(rawe2m2$Site))

## As we see Site_2 levels is not written the same way.
## we can correct this

rawe2m2$Site[rawe2m2$Site == "SITE_2"] <- "Site_2"

length(unique(rawe2m2$Site))

levels(factor(rawe2m2$Site))



########  Date variable   #########

rawe2m2$Date

## There is definitely something wrong with Date variable
## It is showing numbers!
## We usually record data on excel as the usual way, ex: Year-month-day.
## But sometime, it is changed to another format when we try to use diffrent file type.
## Here the date was changed to numaerical, which represent the time spent from the origin date
## The origin date can be different between software used.
## For excel windows "1899-12-30"
## For excel mac "1904-01-01"
## For R "1970-01-01"

## We can foramt this type of date in R, but we need to specify the wright origin

rawe2m2$Date <- as.Date(rawe2m2$Date,origin="1899-12-30")

## Now check it back

rawe2m2$Date


##########  Sex variable  ########
# Same thing here
# Check levels and unique values

length(unique(rawe2m2$Sex))

levels(factor(rawe2m2$Sex))

## as we can see, there are five levels. But there is something wrong with one. We do not know which sex is it.
## For this kind of value (when we do not know what it is) it is better to give it an NA value (non available)
## Let's then fix Sex variable

rawe2m2$Sex[rawe2m2$Sex == "F "] <- "F" ## Correcting value of F with space

rawe2m2$Sex[rawe2m2$Sex == "f"] <- "F" ## Correcting value of F lower case

rawe2m2$Sex[rawe2m2$Sex == ""] <- NA ## Correcting sex with missing value

## Lest check back again

length(unique(rawe2m2$Sex))

## There is still three unique value
## We should not be worried since now R is considering the wrong one as NA. We cannot fix NAs, they will be handle within function.
## Let check it

levels(factor(rawe2m2$Sex))



#########  Weight variable  #####
## as we can see, weight variable is recorded as numeric type

## Let try to change it to numeric then

as.numeric(rawe2m2$Weight)

## I think we are good here


#########  Rainfall variable  #####

rawe2m2$Rainfall
  
## What do you think it is wrong here?
## Of course, the units are still on the values!
## For this, we can fix it by removing al units
## Also, there is "NAmm" values in it.
## Let us first fix those with unit
## And we will be using a package "stringr". This package is usefull for string manipulation
# Load libraries

require(stringr)


rawe2m2$Rainfall <- str_replace(rawe2m2$Rainfall,"mm","")

## Check back again
## All mm are gone.

## Now we will set "NA" as string to correct NA

rawe2m2$Rainfall[rawe2m2$Rainfall=="NA"] <- NA

## check back

rawe2m2$Rainfall

#########  gmam and stes variables

rawe2m2$gmam

rawe2m2$stes

## They are fine for now

## Check back the data

str(rawe2m2)

## Lest then giveright format to our data
rawe2m2$Site <- as.factor(rawe2m2$Site)
rawe2m2$Sex <- as.factor(rawe2m2$Sex)
rawe2m2$Rainfall <- as.numeric(rawe2m2$Rainfall)

## Check it back

str(rawe2m2)

## Save each data step

## Formatted data

e2m2fromated <- rawe2m2



############ We are done with formating the Data!! #########
## But format is not the only type of messy data
## You also need to give consistency to it.
## We will only see one type of consistency since many cases could happen


############## Example ########
## let us check back to the last two variables
## These data are linked to Sex variable. Indeed, gmam (mammary gland should be for female only!)
## The same for stes (Testicule size)
## Let us check then

e2m2fromated$gmam[e2m2fromated$Sex=="M" & !is.na(e2m2fromated$gmam)]

e2m2fromated$stes[e2m2fromated$Sex=="F" & !is.na(e2m2fromated$stes)]

## Waht do you think about it?
## As we can see there are some female who have testicule size values!!
## It is obvious that the when recording the data, some male testicule value were placed into femelle mammary gland values!
## Now it is easy to fix them with simple line

## First we need to move these values

e2m2fromated$gmam[e2m2fromated$Sex=="F" & !is.na(e2m2fromated$stes)] <- e2m2fromated$stes[e2m2fromated$Sex=="F" & !is.na(e2m2fromated$stes)]

## Then give NA to the wrong value

e2m2fromated$stes[e2m2fromated$Sex=="F" & !is.na(e2m2fromated$stes)] <- NA


## Let us check back!!
e2m2fromated$gmam[e2m2fromated$Sex=="M" & !is.na(e2m2fromated$gmam)]

e2m2fromated$stes[rawe2m2$Sex=="F" & !is.na(e2m2fromated$stes)]

## We are good now!!


##############       Thanks A lot !!!!  ##################


