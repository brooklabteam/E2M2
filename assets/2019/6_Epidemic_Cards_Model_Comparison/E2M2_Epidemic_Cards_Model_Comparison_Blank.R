## Epidemic Cards Tutorial: computer exercise
## E2M2: Ecological and Epidemiological Modeling in Madagascar
## January 12-21, 2019

## Cara Brook, 2019

## By now, you've been exposed to some of the variety of model
## types used to understand infectious disease data. We'll now work
## with a few of those forms to explain our own data and compare the
## fits to the data we generated yesterday.


## By the end of this tutorial, you should be able to:
##  * Differentiate between deterministic and stochastic AND discrete and continuous time models.
##  * Load, run, and understand all models of these sorts.
##  * Run a function that fits various models to data by minimizing the sum of squares.
##  * Plot these various fits and compare them

##  * Fit a simple linear regression model and one incorporating random effects
##  * Compare these models performance by AIC

library(plyr)
library(lme4)
library(deSolve)

rm(list=ls())
######################################################################
## Part 1: Visualizing your data

## Earlier, you played an epidemic card game and generated several time
## series of infecteds and susceptibles. 

## (1) Import these data and select the subset for which R0 = 2.
##     Save that subset as an object called: dat.R0
dat <-  read.csv("Epidemic_Cards.csv", header=T)



## (2) Plot each trial as a thin, dashed line (red = infecteds; green = susceptibles)
## like you did yesterday.




######################################################################
## Part 2: Modeling your data in multiple discrete ways

## (3) Brainstorm the different types of models that you might use to
## reconstruct these data. Which type of model did you already build?
## Why did the data sometimes vary from the model predictions? What are
## other types of models you might try to fit?


## (4) Make a vector of values called 'timesteps' that spans from 1 to 10 
## and increases by increments of 1.



## (5) Rewrite your discrete time function from this morning. As before, 
## call it "SIR.discrete.deterministic" and have it take in the input arguments
## of 'R0' and 'time' and return a dataframe of 'time', 'S', and 'I'.
## Save that function to your global environment.



## (6) Check that it is working by running it with the input arguments
## R0=2 and time = timesteps.



## (7) Write a stochastic version of the same model called "SIR.discrete.stochastic"
## with the same input arguments and output dataframe. (Hint, use the built-in R function 
## rnorm() with sd=.5 to add some noise to R0 within each timestep.)


## (8)  Check that it is working by running it with the input arguments
## R0=2 and time = timesteps.



## (9) Run the stochastic model from #8 again above. Does the output look the
## same? Why or why not?

## (10) Save the output from #8 and #9 as the objects 'discrete.determ' and 
## 'discrete.stoch' and plot their output as thick lines atop the plot from #2.
## Color the S class friom discrete.stoch as 'forestgreen' and the I class as 
## orange so you can see that they are different. Which fits the data better?




## (11) Here is the sum_sq function you used for fitting your deterministic 
## model this morning:
sum_sq <- function(par, real.dat, time){
  ## 1. You will first need to run your model at your chosen parameters.
  out.model <- SIR.discrete.deterministic(R0=par, time=time)
  
  ## 2. You will then need to add two columns labeled "I_predictions" and "S_predictions" to your dat.R0
  ##    dataset and fill with 'NAs'
  real.dat$I_predictions <- real.dat$S_predictions <- NA
  
  ## 3. You will then next to fill those columns with model predictions at the appropriate timesteps
  ##    like those from #21. This will require a double for-loop over the unique trials and timesteps.
  
  ## Infecteds first:
  for(i in 1:length(unique(real.dat$trial))){
    for(j in 1:length(unique(real.dat$timestep))){
      real.dat$I_predictions[real.dat$trial==i & real.dat$timestep==j] <- out.model$I[out.model$time==j]
    }
  }
  
  ## And Susceptibles:
  for(i in 1:length(unique(real.dat$trial))){
    for(j in 1:length(unique(real.dat$timestep))){
      real.dat$S_predictions[real.dat$trial==i & real.dat$timestep==j] <- out.model$S[out.model$time==j]
    }
  }
  
  ## 4. You will need to calculate the sum of squared differences between all model predictions (both S
  ##    and I) from the data and return it from the function.)
  
  ## First, get both S and I differences
  real.dat$S_differences <- real.dat$S_predictions - real.dat$susceptibles
  real.dat$I_differences <- real.dat$I_predictions - real.dat$infecteds
  
  ## Then, square them and add them:
  sum.of.sqs = sum(c(real.dat$S_differences, real.dat$I_differences)^2)
  
  return(sum.of.sqs)
  
}
## Adapt it by adding an input argument called 'model' that allows you to choose 
## whether to run either SIR.discrete.deterministic() or SIR.discrete.stochastic().




## (12) Now, use the 'BFGS' method in the function optim() to minimize the difference 
## between your model and the data, and save the returned values from each model
## as objects called 'out.optim.deter' and  'out.optim.stoch'.
## Guess the value of R0 as 2, meaning that you should use the following input variables: 
## optim(par = 2, fn=sum_sq, real.dat=dat.R0, model = "discrete.deterministic", time=seq(1,10,1), method="BFGS")
## optim(par = 2, fn=sum_sq, real.dat=dat.R0, model = "discrete.stochasitic", time=seq(1,10,1), method="BFGS")



## (13) Save out.optim.deter$par as the object R0.discrete.deter.fit and out.optim.stoch$par as R0.discrete.stoch.fit.
## How close are these values to R0=2? Why might they be different?



## (14) Compare out.optim.deter$value and out.optim.stoch$value. Which model
## has the lower sum of squares and better fits the data?


## (15) Re-fit the stochastic model. Do you get the same estimate for R0? Why or 
## why not?


## (16) Run your both your models with their new fitted parameters and write over the 
## output from step #10 with your new predictions. Add these to the plot.




######################################################################
## Part 3: Modeling your data continuously

## (17) Take a look at the following continuous time, deterministic function
## that models the same system as the discrete time model shown above. Adapt the
## function to make a new function called "SIR.continuous.stochastic()" that 
## incorporates stochasticity in the form sd=.5 on the R0 value outside of the sub-ODE.
## The deterministic function:
SIR.continuous.deterministic <- function(R0, time){
  params = list(R0 = R0, gamma=1)
  x.start = c(S=25, I=1)
  ## load the sub-ODE function that the line runs below:
  SIR.ODE <-function(t,x,parms){
    
    #now for the differential equations
    with(c(as.list(c(x,parms))), {
      N = S + I
      
      dS = -R0*((S/N)*I)
      dI = R0*((S/N)*I) - gamma*I
      
      
      
      return(list(c(dS,dI)))
    })
  }
  out <- as.data.frame(lsoda(x.start, time, SIR.ODE, params))
  return(out)
  
  
}



## (18) Test both models at R0=2 to see that they run:



## (19) Do you expect to get the same output every time from both? Try them again to see.



## (20) Edit the sum_sq function from #11 above to allow you to choose to run these
## two additional models too.



## (21) Fit both continuous time models with optim, using the edited sm_sq 
## function from 20. Save output, respectively as out.cont.deter and out.cont.stoch.



## (22) How do these fits compare with the discrete time fits above? Run the fitted
## models and add these new lines to the plot. 




######################################################################
## Part 4: BONUS: TSIR. Modeling your data statistically, but with mechanism

## We've stressed in this course that statistical models tell us about 
## 'correlation' while mechanstic models tell us about 'causation'...But there 
## are a few exceptions. We CAN understand something mechanistic from a 
## statistical model if our x and y regression terms have a mechanistic relationship.

## Recall that our discrete time model for this system is written as:

## S[t+1] = S[t] - (R0 * S[t]/N)*I[t]
## I[t+1]= (R0 * S[t]/N)*I[t]

## Since there is only 1 term in the I equation, we can re-imagine it as:
## y = m * x where:
## y = I[t+1], m = R0, and  x = (S[t]/N)*I[t])
## (b=0 in this case)

## (23) Here, we provide a function that takes a dataset from this model
## and converts it to xy coordinates. Use this function to save a new dataset called
## xy.dat, derived from the R0=2 subset of our original data. (Hint: you will need
## to 'apply' the function over all timesteps in each distinct trial).
## Function:
get.xy <- function(data){
  
  y <- data$infecteds[2:length(data$infecteds)]
  
  x1 <- data$infecteds[1:(length(data$infecteds)-1)]
  x2 <- data$susceptibles[1:(length(data$susceptibles)-1)]
  
  x = x1*x2/26
  
  
  new.dat <- cbind.data.frame(x,y, unique(data$trial))
  names(new.dat) <- c("x", "y", "trial")
  
  return(new.dat)
}

## A. First, make a list of the data separated by trial, using 'dlply'


## B. then apply a function to get the x and y predictors by each trial


## C. Then, and rejoin the list into a dataset


## D. make sure "trial" is a factor to use it as a random effect


## E. and sort by ascending values of Ix


## (24) Plot the values of x vs. y to visualize:


## (25) Fit a simple linear regression, lm(), to these data, and 
## save the output as the object, 'simple.lm'
## Plot a fit line from simple.lm in red.



## (26) Look at the second coefficient of the regression model. 
## What parameter does this signify?


## (27) Examine the AIC. Can we learn anything from this?


## (28) Now fit a model incorporating random effects of 'trial' to these
## data. Save the output as the object 're.lmer'. Plot a fit line
## from this model in 'green'



## Do you notice anything unusual about this line? Why or why not?

## (29) Compare the AIC of re.lmer with simple.lm (above). Which fits the data better?
## Which model would you choose?

