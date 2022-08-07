## Tutorial: Introduction to Model Fitting (Fire Model)
## E2M2: Ecological and Epidemiological Modeling in Madagascar
## January 13-22, 2018

## Cara Brook, 2018

## This tutorial introduces you to the concept of model fitting, here using
## least square techniques to fit both statistical and mechanistic models 
## to simulated data of forest coverage in Madagascar

## By the end of this tutorial, you should be able to:
##
##  * Construct a simple linear regression model and fit it, both manually and via R's built-in functions
##  * Construct a simple mechanistic model, both in discrete and continuous time
##  * Fit the mechanistic model via least squares

setwd("/Users/SarahGuth/Documents/Berkeley/E2M2/Model Fitting 2019")

rm(list=ls())

library(deSolve)

## load and plot your data
load("treedata.Rdata")

head(treedata)

with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 

################ STATISTICAL MODEL ###################

### First, let's try a statistical model on these data

## fit a simple linear regression using R's built-in functions.
m.stat = glm(formula = forest ~ yr, family = gaussian,  data=treedata)

## check out the paramters which this model estimates. 
summary(m.stat)

## we are interested in the slope (m) (-2293) and the y-intercept (b) (5152515.5)
## add the fit line to your graph
with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 
lines(treedata$yr, predict(m.stat), type="l", col="red", lwd=2)

## now write a function to create a new line, using the estimate m and b parameters.
mxb = function(m,x,b){
  y = m*x + b
  return(y)
}
tree.predict = mxb(m=(-2293), x = treedata$yr, b = 5152515.5)

## add to graph
with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 
lines(treedata$yr, predict(m.stat), type="l", col="red", lwd=2)
lines(treedata$yr, tree.predict, type="l", col="magenta", lwd=2)

## write a function to fit your new your new manual model by least squares
lst.sq = function(par, xval, data.forest){
  m.guess = par[1]
  b.guess = par[2]
  
  ##run your model with your guess parameters
  tree.predict = mxb(m=m.guess, x = xval, b = b.guess)
  
  ##compare with data via sum of squares.
  
  sum.sq = sum((tree.predict - data.forest)^2)
  
  return(sum.sq)
  
  
}

## and a wrapper to minimize this function
wrap.fit = function(m.guess, b.guess, xval, data.forest){
  
  par = c(m.guess, b.guess)
  
  
  out = optim(par=par, fn = lst.sq, xval=xval, data.forest=data.forest)
  
  m.fit = out$par[1]
  
  b.fit = out$par[2]

  return(list(m.fit, b.fit))  
}

m.stat.2 = wrap.fit(m.guess = -1000, b.guess = 8*10^5, xval=c(treedata$yr), data.forest = c(treedata$forest))

## and run the model with the estimates.
tree.predict2 = mxb(m=m.stat.2[[1]], x = treedata$yr, b = m.stat.2[[2]])

with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 
lines(treedata$yr, predict(m.stat), type="l", col="red", lwd=2)
lines(treedata$yr, tree.predict, type="l", col="magenta", lwd=2)
lines(treedata$yr, tree.predict2, type="l", col="cyan", lwd=2)

## We managed to recreate the same line many ways! 
## It tells us that forest is declining (negative slope) and gives us an estimate
## of the amount of forest cover at the beginning of our time series (b=~5e10^6 sq. km)
## But we still don't know anything about WHY or HOW that decline took place.

## Let's see if a mechanistic model can help us find out...

