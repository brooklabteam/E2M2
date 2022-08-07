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


rm(list=ls())

library(deSolve)

## load and plot your data
load("treedata.Rdata")

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
lst.sq = function(par,  xval, data.forest){
  m.guess = par[1]
  b.guess = par[2]
  
  ##run your model with your guess parameters
  tree.predict = mxb(m=m.guess, x = xval, b = b.guess)
  
  ##compare with data via sum of squares.
  
  sum.sq = sum((tree.predict - data.forest)^2)
  
  return(sum.sq)
  
  
}

## and a wrapper to minimize this function
wrap.fit = function(guess.slope, guess.int, xguess, data.forest){
  
  par = c(guess.slope, guess.int)
  
  
  out = optim(par=par, fn = lst.sq, xval=xguess, data.forest=data.forest)
  
  m.fit = out$par[1]
  
  b.fit = out$par[2]

  return(list(m.fit, b.fit))  
}

m.stat.2 = wrap.fit(guess.slope = -1000, guess.int = 8*10^5, xguess=c(treedata$yr), data.forest = c(treedata$forest))

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

## Let's see if a mechanistic model can help us find out.

################ MECHANISTIC MODEL ###################
## We think that the REASON behind the declining forest might have to do with slash
## and burn agriculture. We decide to build a model describing conversion of forest
## to savanna via slash and burn. Instead of m and b, the slope and intercept of a line,
## we now want to estimate the the rate of slash and burn. We already know the rate
## of forest growth (r=1.01)
## Let's build and run both a continuous time version mechanistic model.

## Continuous Time Mechanistic Model
ForestSavannaCont <- function(t,y,parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    N = For + Sav
    dFordt <- r*((K-N)/K)*For - slash*For*Sav
    dSavdt <- slash*((K-N)/K)*For*Sav 
    
    #Return forest to compare with data
    list(c(dFordt, dSavdt))
  })
}


## Intialize population
Mada.start <- c(For = 450000, #sq. km of forest in Mada
                Sav = 100000)  #sq. km of savanna in Mada  

## Set up time-steps and units - here go for slightly longer to allow model to equilibrate
times <- seq(1900,2010,by=1) 

## Define parameters 
values <- c(r=1.01,  # the number of sq. km of forest produced by one existing sq. km of forest per year. r=1.01 means that, one year from now, we expect that 1 sq. km of forest will have grown to 1.01 sq. km of forest
            slash = .00000045, # number of sq km of savanna produced per sq km of slashed forest each year. slash = 1 would mean successful 1:1 conversion. slash = .8 would mean some forest gets slashed and still goes back to forest
            K = 900000) # carrying capacity for forest. sq. km of land available for forest/savanna

## And run model
contModTree <- data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, parms = values))
names(contModTree) = c("yr", "forest")

## And plot with data.
with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 
lines(contModTree$yr, contModTree$fores, type="l", col="magenta", lwd=2)


## It looks like it converts forest too quickly!
## We'll need to adjust the rates. We can do this by fitting the models to the data.

################ FITTING MECHANISTIC MODELS ###################

## To do this, we will make three functions:
## (1) A function that runs our model (we already made this above).
## (2) A function that compares our model with the data via a statistical test (we'll
## use least squares)
## (3) A wrapper function that minimizes the difference between our model and the
## data (here, we'll use R's function 'optim' to minimize the function we write 
## under #2)

## Write our comparison function:
sum.sq.mech = function(par, xstart, times, data, r, K){
  ## run either your discrete or continuous model
    parms = c(r = r,
              slash = par,
              K = K)
    
    model.out = data.frame(lsoda(y = xstart, times = times, func = ForestSavannaCont, parms = parms))
    names(model.out) = c("yr", "forest")
    
  ## remember, that we ran our model longer than the time period for which we have data,
  ## so we can't compare the earlier timesteps. we can subset the model output to compare
  ## only those chunks for which we have data.
    
  ## find the year your data starts.
    min.dat = min(data$yr)
    
    model.out = model.out[model.out$yr >= min.dat,]
    
  ## and compare the output of the model with the data
    sum.sq = sum((model.out$forest - data$forest)^2)
    
    return(sum.sq)
    
}

##test it on one run of your model
sum.sq.mech(par=.00000045,
            xstart=Mada.start,
            times = times,
            r=1.01,
            K= 900000,
            data=treedata)
            
## Write our wrapper function
wrap.fit.mechanistic = function(xstart, times, r, K, data){
  
  ## make a sequence of guesses for the rate
  slash.list = seq(0, .000001, .0000001)
  
  ## then, search for the minimum
  sum.sq.lst = list()
  for (i in 1:length(slash.list)){
    sum.sq.lst[[i]] = sum.sq.mech(par=slash.list[i], xstart=xstart, r=r, K=K, times=times, data=data)
  }
  sum.sq.lst = c(unlist(sum.sq.lst))
  
  plot(slash.list, sum.sq.lst, type="b")
  
  ##find the minimum
  fit.slash = slash.list[sum.sq.lst==min(sum.sq.lst)]
  
  return(fit.slash)
}

## Fit the continuous model. 

fit.continuous = wrap.fit.mechanistic(xstart=Mada.start,
                                    times = times,
                                    r=1.01,
                                    K= 900000,
                                    data=treedata)


## Run your continuous model with your fitted parameters and plot with the data
new.values.cont= c(r=1.01, 
          slash = fit.continuous,
          K = 900000)

run.fit.continuous.new  <- data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, parms = new.values.cont))
names(run.fit.continuous.new) = c("yr", "forest")

with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 
lines(contModTree$yr, contModTree$forest, type="l", col="yellow", lwd=2)
lines(run.fit.continuous.new$yr, run.fit.continuous.new$forest, type="l", col="orange", lwd=2)
legend("bottomleft", legend=c("data", "model (guess)", "model (fitted)"), col=c("green", "yellow", "orange", lwd=2))

