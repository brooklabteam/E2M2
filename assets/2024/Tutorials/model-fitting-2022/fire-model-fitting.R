# E2M2 2022: Model Fitting Basics
# Project: 
# Authors: 

# Michelle V Evans
# Github: mvevans89
# Institut de Recherche pour le Developpement
# Email: evans.michelle@ird.fr

# Script originated Nov 2022

# Description of script and Instructions #####################

#' This is a tutorial that was developed as part of the E2M2 workshop held in December 2022. 
#' It is based on code developed by Cara Brook, further developed by Michelle Evans. 
#' If you have any questions, please send an email to Michelle at mv.evans.phd@gmail.com.

#' This tutorial introduces the concept of model fitting. We use the term "model fitting" 
#' broadly, referring to both statistical and mechanistic models. In this tutorial you will 
#' learn the least-squares technique to fit both types of models to simulated data of forest
#' coverage in Madagascar.

#' By the end of the tutorial, you should be able to:
#'  - Construct a simple linear regression model and fit it, both manually and via R's built-in functions
#'  - Construct a simple mechanistic model in continuous time
#'  - Fit a mechanistic model via least squares


# Packages & Options #########################################

options(stringsAsFactors = F, scipen = 999, digits = 4)

#plotting
library(ggplot2); theme_set(theme_bw())

#ordinary differential equations
library(deSolve)

#data manipulation
library(tidyr)
library(dplyr)

# Initial Exploration of the Data ############################

#simulated forest cover in Madagascar over time
treedata <- read.csv("treedata.csv")

head(treedata)
str(treedata)

#plot the data
base.plot <- ggplot(treedata, aes(x = yr, y = forest)) +
  geom_point(color = "forestgreen") +
  ylim(c(0,8e5)) +
  xlab("Year") +
  ylab(bquote('Forest Cover'(km^2)))

base.plot

# Statistical Model #####################################

## Fit Using glm ########################################

#fit a linear regression, y = mx + b
mod.stat <- glm(formula = forest ~ yr,
                family = gaussian,
                data  = treedata)

#investigate parameters estimated by model
summary(mod.stat)

#these parameters are the Intercept (b) and the coefficient of the year (slope, m)
#view coefficients
mod.stat$coefficients

#add fit trend-line to plot
base.plot +
  geom_abline(slope = mod.stat$coefficients[2], intercept = mod.stat$coefficients[1], color = "darkred")

#manually create predictions for a trend-line using existing `predict` function
base.predict <- data.frame(yr = treedata$yr,
                           prediction = predict(mod.stat)) 

head(base.predict)

#write a function to create a new line from the estimated m (slope) and b (intercept)
mxb_line <- function(m, x, b){
  y = m*x + b
  return(y)
}

#example of predicting with this
my.predict <- data.frame(yr = treedata$yr,
                         prediction = mxb_line(m = mod.stat$coefficients[2],
                                               x = treedata$yr, 
                                               b = mod.stat$coefficients[1])) 

head(my.predict)

#add both predicted lines to the plot to compare
base.plot + 
  geom_line(data = base.predict, aes(y = prediction), color = "darkred", size = 1) +
  geom_line(data = my.predict, aes(y = prediction), color = "dodgerblue", linetype = "dashed",
            size = 1)

## Manually fit by least-squares ############################

#' There are two parts to manually fitting a model, corresponding to two functions.
#' The first function takes the parameters for `m` and `b`, predicts y-values based
#' on these parameters, and calculates the sum of squares wit the true data. 
#' The second function wraps around the first function to minimize the sum of 
#' squares, resulting in the optimum estimate for the parameters.

### First Function #########################################

#' The first function to calculate a sum of squares for one set of parameters. 

lst_sq <- function(par,  xval, ytrue){
  m.guess <- par["m"]
  b.guess <- par["b"]
  
  ##run your model with your guess parameters
  y.pred <- mxb_line(m = m.guess, x = xval, b = b.guess)
  
  ##compare with data via sum of squares.
  sum.sq <- sum((y.pred - ytrue)^2)
  
  return(sum.sq)
  
}

#' as an example, let's compare the least squares for the optimal parameters
#' to some that we have manually chosen

#optimal parameters from linear regression
lst_sq(par = c("m" = -2293, "b" = 5152515), xval = treedata$yr, ytrue = treedata$forest)

#estimate SS with parameters lightly moved from optimum. a worse fit
lst_sq(par = c("m" = -2600, "b" = 5500000), xval = treedata$yr, ytrue = treedata$forest)

# we can visualize these fit lines too
base.plot +
  geom_abline(intercept = 5152515, slope = -2293) +
  geom_abline(intercept = 5500000, slope = -2600, color = "dodgerblue")

### Second Function ########################################

#' The second function is a wrapper that works to minimize the output 
#' (i.e. the sum of squares) from the previous function. By minimizing the 
#' output, we will find the optimum fit. This is done using the `optim` function 
#' in R, which searches the parameter space to find the values that result in 
#' the minimum sum of squares.

wrap_fit = function(guess.slope, guess.int, xguess, ydata){
  
  par <- c("m" = guess.slope, "b" = guess.int)
  
  
  out <- optim(par=par, fn = lst_sq, xval = xguess, ytrue = ydata)
  
  m.fit = out$par['m']
  
  b.fit = out$par['b']
  
  return(list(m.fit, b.fit))  
}

#' The two "guesses" that we give to this function are the initial 
#' starting values that the `optim` function will start the search at.
#' What happens if we give poor starting values?

poor.fit <- wrap_fit(guess.slope = 10, guess.int = 100, 
                     xguess = treedata$yr, ydata = treedata$forest)
poor.fit

base.plot +
  geom_abline(slope = poor.fit[[1]], intercept  = poor.fit[[2]])

#' we get a bad fit. the starting guess is important and can be estimated 
#' by plotting data. what if we use a good fit?

good.fit <- wrap_fit(guess.slope = -5000, guess.int = 800000, 
                     xguess = treedata$yr, ydata = treedata$forest)

good.fit

#predict from this model
good.predictions <- mxb_line(m = good.fit[[1]], x = treedata$yr, b = good.fit[[2]])

#plot this data
base.plot +
  geom_line(aes(x = treedata$yr, y = good.predictions))

### Compare all methods ####################################

#create a plot to compare all the linear regressions
base.plot +
  #original statistical model with glm
  geom_line(aes(x = treedata$yr, y = predict(mod.stat)), color = "darkred") +
  #model fit using optim
  geom_line(aes(x = treedata$yr, y = good.predictions), color = "blue", 
            linetype = "dotted", size = 1)

#' From this we can see that forest is declining, because the slope is 
#' negative, and also estimate he amount of forest cover at the beginning 
#' of our time series using any of our predict functions. Remember that the 
#' intercept returned by the model isn't the intercept at the beginning of 
#' our data in 1920, but would correspond to the amount of forest cover when 
#' x = 0 (which isn't very useful to us), so you need to predict it using the 
#' value of the year.

#' This is all great, but we don't know anything about why or how this decline 
#' took place. For this we need to use a mechanistic model.
#' 

# Mechanistic Model #################################################

#' Mechanistic models allow us to directly model the process that may 
#' be creating the data that we have. For example, we may hypothesize 
#' that the declining forest is due to slash and burn agriculture. We decide 
#' to build a model describing conversion of forest to savanna via slash and 
#' burn. Instead of m and b, the slope and intercept of a line, we now want 
#' to estimate the rate of slash and burn.

# define the model as a function that can be read by the ODE solver
ForestSavannaCont <- function(t,y,parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    N = For + Sav
    dFordt <- r*((K-N)/K)*For - gamma*For*Sav
    dSavdt <- gamma*((K-N)/K)*For*Sav 
    
    #Return forest to compare with data
    list(c(dFordt, dSavdt))
  })
}

#' to run the model, we need to provide the initial conditions (amount 
#' of forest and savanna), the time series we want to run it over, and values 
#' for the parameters

## Intialize population
Mada.start <- c(For = 450000, #sq. km of forest in Mada
                Sav = 100000)  #sq. km of savanna in Mada  

## Set up time-steps and units - here go for slightly longer to allow model to equilibrate
times <- seq(1900,2010,by=1) 

## Define parameters 
values <- c(r=1.01, 
            gamma = .00000045, 
            K = 900000) 

#run the model and plot
slash.mod1<- data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, 
                              parms = values))
names(slash.mod1) = c("yr", "forest", "savanna")

base.plot +
  geom_line(data = slash.mod1, aes(x = yr, y = forest))

## Fitting a mechanistic model ######################################

#' Like the statistical model, we use a workflow consisting of several 
#' functions to fit a mechanistic model:

#' 1. A function that runs our model (we already made this above `ForestSavannaCont`)
#' 2. A function that compares our model with the data via a statistical test 
#' (we'll use least squares)
#' 3. A wrapper function that minimizes the difference between our model and 
#' the data (here, we'll use R's function 'optim' to minimize the function we write under #2)

#' The comparison function takes a set of parameters, estimates the the forest 
#' cover for each year, and compares these predictions to the true data using 
#' the sum of least squares:

sum_sq_mech <- function(gamma, r, K, xstart, times, data){
  parms <- c(r = r,
             gamma = gamma,
             K = K)
  
  model.out <- data.frame(lsoda(y = xstart, times = times, 
                                func = ForestSavannaCont, parms = parms))
  names(model.out) <- c("yr", "forest", "savanna")
  
  # the model runs for longer than our data, so we need to make sure 
  # we are comparing the same years across the two datasets
  forest.preds <- model.out[match(data$yr, model.out$yr),]
  
  ## and compare the output of the model with the data
  sum.sq <- sum((forest.preds$forest - data$forest)^2)
  
  return(sum.sq)
}

#test on one value of parameters
sum_sq_mech(gamma = 0.00000045,
            r = 1.01,
            K = 900000,
            xstart = Mada.start,
            times = times,
            data = treedata)

#' create the wrapper function. it explores a range of values for the slashing
#' rate and choose the one which results in the lowest sum of squares
wrap_mech <- function(xstart, times, r, K, data){
  
  ## make a sequence of guesses for the rate
  slash.list = seq(from = 0, to = .000001, by = .0000001)
  
  ## then, search for the minimum
  sum.sq.lst = list()
  for (i in 1:length(slash.list)){
    sum.sq.lst[[i]] = sum_sq_mech(gamma=slash.list[i], xstart=xstart, r=r, K=K, times=times, data=data)
  }
  sum.sq.lst = c(unlist(sum.sq.lst))
  
  plot(slash.list, sum.sq.lst, type="b")
  
  ##find the minimum
  fit.slash = slash.list[sum.sq.lst==min(sum.sq.lst)]
  
  return(fit.slash)
}

#fit the model
fit.mech <- wrap_mech(xstart = Mada.start, times = times, r = 1.01, 
                      K = 900000, data = treedata)
#it returns a plot and the value of the best fit
fit.mech

#use the estimated parameter to fit the mechanistic model
param.fit <- c(r = 1.01, gamma = fit.mech, K = 900000)

mech.fit <- data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, 
                             parms = param.fit))
names(mech.fit) <- c("yr", "forest", "savanna")

#add this line to the plot, compare to the initial one that wasn't fit
base.plot +
  geom_line(data = slash.mod1, aes(x = yr, y = forest)) +
  geom_line(data = mech.fit, aes(x = yr, y = forest), color = "dodgerblue") +
  xlim(c(1920, 2010))
