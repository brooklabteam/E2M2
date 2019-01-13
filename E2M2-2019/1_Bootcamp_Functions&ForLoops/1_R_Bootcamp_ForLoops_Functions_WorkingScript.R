

## Let's work together to re-analyze Christian's dataset:

setwd("/Users/caraebrook/iCloud-Drive-Dec-2018/Documents/R/R_repositories/E2M2/E2M2-2019/1_Bootcamp_Functions&ForLoops")
##############################################
# Load data and visualize
##############################################
## First, read in the .csv file.

dat <- read.csv(file="e2m2.csv", header = T, stringsAsFactors = F)
head(dat)
str(dat)

## What do we do if we want to print out the Forearm length of every bat in the dataset?

dat$Forearm
dat['Forearm']


##############################################
# For-loops
##############################################

## Can we use a for-loop to do this? Demonstrate.

## Remember the syntax:

## for (variable in vector) {
# do something
# }

length(dat$Id)

#for(i in 1:length(dat$Id)){


for(i in 1:1200){
  print(dat$Forearm[i])
}

for(p in 1:3){
  print(dat$Forearm[p])
}


i = 3
print(dat$Forearm[i])



##############################################
# If-Else Statements
##############################################

## What if we only want to print the forearm for the bats that are positive for infection?

## Remember the syntax:

## if (condition is TRUE) {
# do something
# } else {
# do different thing
# }

head(dat)
dat$Pos
unique(dat$Pos)

if('bats are positive for babesia'){
  'print forearm length'
}else{
  'print "bats not infected."'
}

interpret.bat <- function(data){
  if (data$Pos==1){
    print(data$Forearm)
  }else{
    print("bats not infected")
  }
}


for(i in 1:1200){
  interpret.bat(dat[i,])
}

library(dplyr)
library(plyr)
dat.list <- dlply(dat, .(Id))
interpret.bat(dat)


lapply(X=dat.list, FUN=interpret.bat)


for(i in 1:1200){
  if(dat$Pos[i]==1){
    print(dat$Forearm[i])
  }else{
    print("bats not infected.")
  } 
  }



##############################################
# Functions
##############################################


## In the above section, we learned how to for-loop over our data and apply useful built-in functions in R (i.e. 'mean').
## Fidy showed you a list of some of these fucntions in his lecture this morning. Try to build a vector and use some of these
## functions to accomplish tasks in R. Calculate the min, max, and mean of the bat forearms in Christian's data, and
## tell us how many there are.

min(dat$Forearm)
max(dat$Forearm)
mean(dat$Forearm)



## Next, can you write an alternative function to the default in R that calculates the mean of your vector? Try it.
## Remember that you will first need to highlight and run your function to save it into your working directory, and
## you'll then be able to use it later.


mean2 <- function(data){
  sum.of.values <- sum(data)
  number.of.values <- length(data)
  new.mean <- sum.of.values/number.of.values
  return(new.mean)
}

mean.forearm.positive <- function(data){
  sum.of.values <- sum(data$Forearm[data$Pos==1])
  number.of.values <- length(data$Forearm[data$Pos==1])
  new.mean <- sum.of.values/number.of.values
  return(new.mean)
}
mean.forearm.negative <- function(data){
  sum.of.values <- sum(data$Forearm[data$Pos==0])
  number.of.values <- length(data$Forearm[data$Pos==0])
  new.mean <- sum.of.values/number.of.values
  return(new.mean)
}

mean2(data = c(1,2))

mean(dat$Forearm)
mean2(dat$Forearm)

mean.forearm.positive(dat)

mean.forearm.negative(dat)



## Remember the syntax:

## All functions in R take on the following general syntax:
#function_name <- function(arguments) {
# body
#}

## And, finally, can you write a third function to get the mean of the vector, but without using any built-in 
## functions in R, specifically sum()? (Hint: Try using a for-loop)





#############################################################################
## Apply Functions

## Many of the built-in functions in R are secretly running C-code (ususally)
## for-loops under the hood, but translating them to English so that you 
## don't actually have to think about it. One of the most powerful
## tools you can use in R, is the "apply" family of functions, which can be used
## to perform a general function over a list of sub-datasets in a linearized fashion.

## For example, what if we want to calculate the mean of female bat forearms and 
## the mean of male bat forearms separately?

## Try it with indexing.

## Try it with a for-loop.

## Try it by applying a function.




