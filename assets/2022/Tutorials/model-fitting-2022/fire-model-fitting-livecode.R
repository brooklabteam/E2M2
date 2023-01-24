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

#' This is the live-coding version, so there are some things that will need to be filled in for
#' the code to run (usually in 'CAPITALS'). There is a filled out script in `fire-model-fitting-livecode.R`.

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
#read in the data set in treedata.csv and save as an object called treedata
treedata <- 'CSV NAME'

#how can we inspect our data to see what is in each of the columns?
'FUNCTIONS TO INSPECT DATA (head and str)'

#plot the data
#what are the x and y variables?
base.plot <- ggplot(treedata, aes(x = 'X VARIABLE', y = 'Y VARIABLE')) +
  geom_point(color = "forestgreen") +
  ylim(c(0,8e5)) +
  xlab("Year") +
  ylab(bquote('Forest Cover'(km^2)))

base.plot

# Statistical Model #####################################

## Fit Using glm ########################################

#fit a linear regression, y = mx + b
#make sure you include the formula, family and data
mod.stat <- 'STAT MODEL FX'

#investigate parameters estimated by model
summary(mod.stat)

#' these parameters are the Intercept (b) and the coefficient of the year (slope, m)
#' how can we view the parameters (also called coefficients) using the operator symbol `$`
'INDEX COEFFICIENTS'

#add fit trend-line to plot
#fill in the proper coefficient for the slope and the intercept using abline
base.plot +
  geom_abline('FILL IN ABOVE COEFFICENTS')

#manually create predictions for a trend-line using existing `predict` function
base.predict <- data.frame(yr = treedata$yr,
                           prediction = 'ADD PREDCITIONS HERE') 

#' let's compare the predictions and original data using `head`
'COMPARISON CODE'

#write a function to create a new line from the estimated m (slope) and b (intercept)
#what arguments does the function need?
#what goes inside the function?
mxb_line <- function(m, x, b){
  y = 'FORMULA HERE'
  return(y)
}

#example of predicting with this
my.predict <- data.frame(yr = treedata$yr,
                         prediction = 'PREDICT USING MXB_LINE') 

head(my.predict)

#add both predicted lines to the base plot to compare
# plot one line using linetype = "dashed" so we can see it
base.plot + 
  geom_line('ADD MY PREDICT') +
  geom_line('ADD BASE PREDICT')

## Manually fit by least-squares #############################

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
  #what will be m, x, and b?
  y.pred <- mxb_line(m = , x = , b = )
  
  ##compare with data via sum of squares
  #how do we calculate the sum of squares?
  sum.sq <- 'SUMSQ EQUATION'
  
  return(sum.sq)
  
}

#' as an example, let's compare the least squares for the optimal parameters
#' to some that we have manually chosen

#optimal parameters from linear regression
#fill in with our optimal parameters from above (hint: use summary(mod.stat))
lst_sq(par = c("m" = 'BEST_M', "b" = 'BEST_B'), xval = treedata$yr, ytrue = treedata$forest)

#estimate SS with parameters lightly moved from optimum. a worse fit
#change the parameters and see what the sum of squares is on the return
lst_sq(par = c("m" = , "b" = ), xval = treedata$yr, ytrue = treedata$forest)

# we can visualize these fit lines too
#use those parameters above in geom_abline
base.plot +
  geom_abline('LINE WITH BEST FIT') +
  geom_abline('LINE WITH MOVED FROM OPTIMUM')

### Second Function ########################################

#' The second function is a wrapper that works to minimize the output 
#' (i.e. the sum of squares) from the previous function. By minimizing the 
#' output, we will find the optimum fit. This is done using the `optim` function 
#' in R, which searches the parameter space to find the values that result in 
#' the minimum sum of squares.

wrap_fit = function(guess.slope, guess.int, xguess, ydata){
  
  par <- c("m" = guess.slope, "b" = guess.int)
  
  # what function are we using to fit this?
  out <- optim(par=par, fn = 'SUM OF SQ FX' , xval = xguess, ytrue = ydata)
  
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
#' fill in with some "better" starting parameters

good.fit <- wrap_fit(guess.slope = 'GUESS', guess.int = 'GUESS', 
                     xguess = treedata$yr, ydata = treedata$forest)

good.fit

# How would we plot this line (geom_abline) on our base_plot?
'CODE TO ADD THE LINE TO BASE.PLOT'

#predict from this model
#what is m and b?
good.predictions <- mxb_line(m = 'GOOD M', x = treedata$yr, b = 'GOOD B')

#plot this data
base.plot +
  geom_line(aes(x = treedata$yr, y = good.predictions))

### Compare all methods ####################################

# create a plot to compare all the linear regressions
# we want to compare the output of predictions from mod.stat to our predictions
# from the fit using sum of squares


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
    #write out the ODE for each population (using the equations from the lecture)
    'ODEs GO HERE'
    
    #Return forest to compare with data
    list(c(dFordt, dSavdt))
  })
}

#' to run the model, we need to provide the initial conditions (amount 
#' of forest and savanna), the time series we want to run it over, and values 
#' for the parameters

## Intialize populations
#create a vector of starting values for each population
Mada.start <- 'VECTOR OF TWO STARTING VALUES'

## Set up time-steps and units - here go for slightly longer to allow model to equilibrate
# use `seq` to create a vector years from 1900-2010
times <- 'SEQUENCE OF YEARS'

## Define parameters 
# we have three parameters in our model that we need to provide values for
values <- 'STARTING PARAMETERS'

#run the model and plot
#fill in the arguments for lsoda (y, times, func, parms)
slash.mod1<- data.frame(lsoda('LSODA ARGUMENTS'))
names(slash.mod1) = c("yr", "forest", "savanna")

#add the fit from this to our base.plot
base.plot +
  geom_line()

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
  sum.sq <- 'SUM OF SQ EQUATION'
  
  #what should the function return? what are we using to measure model performance?
  return('RETURN')
}

#test on one value of parameters
sum_sq_mech(gamma = 0.0000000045,
            r = 1.01,
            K = 900000,
            xstart = Mada.start,
            times = times,
            data = treedata)

#' create the wrapper function. it explores a range of values for the slashing
#' rate and choose the one which results in the lowest sum of squares
wrap_mech <- function(xstart, times, r, K, data){
  
  ## make a sequence of guesses for the rate (use seq)
  slash.list = 'RATES TO TEST'
  
  ## then, search for the minimum
  sum.sq.lst = list()
  for (i in 1:length(slash.list)){
    sum.sq.lst[[i]] = sum_sq_mech(gamma=slash.list[i], xstart=xstart, r=r, K=K, times=times, data=data)
  }
  sum.sq.lst = c(unlist(sum.sq.lst))
  
  plot(slash.list, sum.sq.lst, type="b")
  
  ##how do we identify the minimum sum of squares?
  fit.slash = 'GET MINIMUM'
  
  return(fit.slash)
}


#fit the model
fit.mech <- wrap_mech(xstart = Mada.start, times = times, r = 1.01, 
                      K = 900000, data = treedata)
#it returns a plot and the value of the best fit
fit.mech

#use the estimated parameter to fit the mechanistic model
param.fit <- 'BEST FIT PARAMETERS'

mech.fit <- data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, 
                             parms = param.fit))
names(mech.fit) <- c("yr", "forest", "savanna")

#add this line to the plot, compare to the initial one that wasn't fit (slash.mod1)
base.plot +
  'ADD FIT LINES'
