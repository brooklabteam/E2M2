#######################################
## Introduction 
#######################################
## Modeling vector-borne diseases R tutorial
## E2M2 - 2019
## written by: Amy Wesolowski, Ben Rice

## In this tutorial, we will construct a simple vector-borne disease model (humans: SIRS, vector: SI). First, we will make a function that is just an SIRS model:

##############################################################################
## Contents: 
# Line 24: SIRS model
# Line 67: SIRS human model, SI vector model
# Line 133: Exercise: Different recovery rates
# Line 163: Different waning immunity rates
# Line 174: Temperature dependent mortality
##############################################################################
  
## load your libraries first

library(deSolve)


##############################################################################
## SIRS Model
##############################################################################

# Write a function that defines the parameters in our SIRS model

sirs.model<-function(t, y, parms){
  with(c(as.list(y), parms),{
    dSdt <- -betaH*S*I + w*R
    dIdt <- betaH*S*I - gamma*I
    dRdt <- gamma*I - w*R
    list(c(dSdt, dIdt, dRdt))
  })}

## NOTES: This function has an additional parameter, w, which measures waning immunity: 
## i.e. the rate at which recovered individuals loose their immunity and become 
## susceptible again. Note: here we have just constructed just a basic SIRS model 
## -- without a vector component. Next, we will investigate how changes in the waning 
## immunity rate change the size of the outbreak. 

## We will run 3 simulations using the same initial conditions and other parameters 
## and only change w. 

## We will run each simulation for 1 year and compare the results.

## set all of your parameters and run your 3 models
## We decided to start with parameters:
    ## population of size 1000
    ## 1 infected individual
    ## 0 recovered individuals
## We are also going to start with:
    ## betaH (the transmission rate) = 0.0001 infected indvidual per unit time (ie day)
    ## gamma (the rate of recovery) = 1/14 (ie it takes two weeks to recover)
## But we are going to vary w (the rate of loss of immunity) and try 3 different rates:
    ## w = 1/10  (it takes  10 days to lose immunity and become susceptible)
    ## w = 1/50  (it takes  50 days to lose immunity and become susceptible)
    ## w = 1/100 (it takes 100 days to lose immunity and become susceptible)

initial.pops<-c(S = 1000, I = 1, R = 0)
parameters.1<-c(betaH = 0.0001, gamma = 1/14, w = 1/10)
parameters.2<-c(betaH = 0.0001, gamma = 1/14, w = 1/50)
parameters.3<-c(betaH = 0.0001, gamma = 1/14, w = 1/100)

## Now let's run the model for 365 days  (1 year), 1 day at a time
    ## Let's store the results we get fromm the 3 parameter scenarios we defined above 
    ## in a variable "result"

times<-seq(0,365,by=1)
result.1<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.1))
result.2<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.2))
result.3<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.3))

## TO DO: 
#### 1) Now, we can plot these results to see the effect of changing the waning immunity 
####    between the three simulations. What differences do you expect? What would you 
####    expect if w = 1/2?
#### 2) What would you expect if w = 1/200? 

## First, plot your results (the number of susceptible individuals over time):

## Use par(mfrow()) to allow showing 2 plots next to each other
par(mfrow=c(1,2))

## Create a blank plot, that will later contain our results 
plot(NA, NA, 
     xlim = c(range(times)), ylim = range(c(result.1$S, result.2$S, result.3$S)), 
     xlab = 'Time', ylab = 'Susc')

## Add the lines from our results and susceptible individuals and add different colors
#### Note that we are plotting the column result.1$S (where S is our susceptible inds)
lines(result.1$S, col = 'red')
lines(result.2$S, col = 'orange')
lines(result.3$S, col = 'blue')

## Add the legend
legend('topright', 
       legend = c('w=1/10', 'w=1/50', 'w=1/100'), 
       col = c('red', 'orange', 'blue'), 
       pch = 15, bty = 'n')

## Now let's plot the infected individuals over time:
plot(NA, NA, 
     xlim = c(range(times)), ylim = range(c(result.1$I, result.2$I, result.3$I)), 
     xlab = 'Time', ylab = 'Infect')
lines(result.1$I, col = 'red')
lines(result.2$I, col = 'orange')
lines(result.3$I, col = 'blue')
legend('topleft', 
       legend = c('w=1/10', 'w=1/50', 'w=1/100'), 
       col = c('red', 'orange', 'blue'), 
       pch = 15, bty = 'n')


#***************************************************************************************#
## FROM ABOVE:
## **TO DO: 
#### 1) Now, we can plot these results to see the effect of changing the waning immunity 
####    between the three simulations. What differences do you expect? What would you 
####    expect if w = 1/2?
#### 2) What would you expect if w = 1/200? 

## Add results from w = 2 and w = 1/200 to the plots from above:
## Note that this code will need to be adjusted:

initial.pops<-c(S = 1000, I = 1, R = 0)
parameters.1<-c(betaH = 0.0001, gamma = 1/14, w = 1/10)
parameters.2<-c(betaH = 0.0001, gamma = 1/14, w = 1/50)
parameters.3<-c(betaH = 0.0001, gamma = 1/14, w = 1/100)
#* Add parameters.4 and parameters.5 where w = 1/2 and w = 1/200

times<-seq(0,365,by=1)
result.1<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.1))
result.2<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.2))
result.3<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.3))
#* Add result.4, result.5 that use parameters.4, parameters.5

par(mfrow=c(1,2))
plot(NA, NA, 
     xlim = c(range(times)), ylim = range(c(result.1$S, result.2$S, result.3$S)), 
     xlab = 'Time', ylab = 'Susc')
#* add result.4 and result.5 to the plot of Susceptibles

lines(result.1$S, col = 'red')
lines(result.2$S, col = 'orange')
lines(result.3$S, col = 'blue')
#* add result.4 and result.5 
    ## (Hint: remember to choose different colors for result.4 and result.5
    #e.g. of colors in R can be found at http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf
    #e.g. chocolate4, coral3, blueviolet, deepskyblue, darkgreen, darkorchid, gold

legend('topright', 
       legend = c('w=1/10', 'w=1/50', 'w=1/100'), #*add w = 1/2, w = 1/200
       col = c('red', 'orange', 'blue'), #*add your colors from above
       pch = 15, bty = 'n')

#Now for Infecteds
plot(NA, NA, 
     xlim = c(range(times)), ylim = range(c(result.1$I, result.2$I, result.3$I)), 
     xlab = 'Time', ylab = 'Infect')
#* add result.4 and result.5 to the plot of Infected

lines(result.1$I, col = 'red')
lines(result.2$I, col = 'orange')
lines(result.3$I, col = 'blue')
#* add result.4 and result.5 and the colors you used in the Susceptible plot

legend('topleft', 
       legend = c('w=1/10', 'w=1/50', 'w=1/100'), #*add w = 1/2, w = 1/200
       col = c('red', 'orange', 'blue'), #*add your colors from above
       pch = 15, bty = 'n')

## What do we observe?

# Which w value gives the most infections at the end of the year?
# Which w value gives the least susceptible individuals at the end of the year?







##############################################################################
##############################################################################
##############################################################################
## SIRS human model, SI vector model
##############################################################################
##############################################################################
##############################################################################

## NOTE: Now we will use the same structure as an SIRS model for humans, 
## but also add in an SI model for vectors. 

## NOTE: now transmission for humans will depend on mosquitoes and visa versa!

vector.trans.model<-function(t, y, parms){
  with(c(as.list(y), parms),{
    dSHdt <- -betaH*SH*IV + w*RH
    dIHdt <- betaH*SH*IV - gamma*IH
    dRHdt <- gamma*IH - w*RH
    
    dSVdt <- mu - betaV*SV*IH - sigma*SV
    dIVdt <- betaV*SV*IH - sigma*IV
    list(c(dSHdt, dIHdt, dRHdt, dSVdt, dIVdt))
  })
}

## NOTES: This function has a number of additional parameters including: 
    ## transmission for humans (betaH)
    ## waning immunity (w)
    ## recovery rate humans (gamma)
    ## vector birth rate (mu)
    ## transmission for vectors (betaV)
    ## a death rate for vectors (sigma)

## Here, we are assuming that vectors do not recover! 
## They can only be in either the susceptible or infected class. 

## Now, we can run a simulation first assuming that the transmission between 
## humans and vectors and between vectors and humans is equal. 

## set your parameters and run your model:
## (refer to the vector-borne lecture slides if needed):
initial.human.vector.pops <- c(
  SH = 1000,    #Susceptible humans
  IH = 1,       #Infected humans
  RH = 0,       #Recovered humans
  SV = 1e6,     #Susceptible vectors
  IV = 0)       #Infected vectors
same.trans.vector.parameters <- c(
  betaH = 0.00005, #Transmission rate for humans
  betaV = 0.00005, #Transmission rate for vectors
  gamma = 1/14,    #Rate of recovery
  w = 1/60,        #Rate of loss of immunity (waning immunity)
  mu = 100,        #Birth rate of mosquitoes
  sigma = 1/14)    #Death rate of mosquitoes

result.same.trans<-data.frame(lsoda(
  y=initial.human.vector.pops, 
  times = times, 
  func = vector.trans.model, 
  parms = same.trans.vector.parameters))

## And we can plot the results from both the vector and human populations:

## plot your results
par(mfrow=c(1,2))
plot(NA, NA, xlim = range(times), ylim = range(c(result.same.trans$SH, result.same.trans$IH)), xlab = 'Days', ylab = '', main = 'Human')
legend('right', legend = c('Susceptible', 'Infected'), col = c('blue', 'red'), pch = 15, bty = 'n')
lines(result.same.trans$SH, col = 'blue', lwd = 2)
lines(result.same.trans$IH, col = 'red', lwd = 2)

plot(NA, NA, xlim = range(times), ylim = range(c(result.same.trans$SV, result.same.trans$IV)), xlab = 'Days', ylab = '', main = 'Vector')
legend('right', legend = c('Susceptible', 'Infected'), col = c('blue', 'red'), pch = 15, bty = 'n')
lines(result.same.trans$SV, col = 'blue', lwd = 2)
lines(result.same.trans$IV, col = 'red', lwd = 2)

## Now, we will change the transmission variables and assume that the transmission rate 
## between humans and vectors (betaV) is less than between vectors and humans (betaH)

## set the different transmission parameters (where betaV does not equal betaH) 
## and run your model:
diff.trans.vector.parameters <- c(
  betaH = 0.00005, #Note that betaH (transmission for humans)  = 0.00005, was 0.0005 before
  betaV = 0.00001, #Note that betaV (transmission for vectors) = 0.00001, was 0.0005 before
  gamma = 1/14, 
  w = 1/60, 
  mu = 100, 
  sigma = 1/14)

result.diff.trans<-data.frame(lsoda(
  y=initial.human.vector.pops, 
  times = times, 
  func = vector.trans.model, 
  parms = diff.trans.vector.parameters))

## plot your results
par(mfrow=c(1,2))
plot(NA, NA, xlim = range(times), ylim = range(c(result.diff.trans$SH, result.diff.trans$IH)), xlab = 'Days', ylab = '', main = 'Human')
legend('right', legend = c('Susceptible', 'Infected'), col = c('blue', 'red'), pch = 15, bty = 'n')
lines(result.diff.trans$SH, col = 'blue', lwd = 2)
lines(result.diff.trans$IH, col = 'red', lwd = 2)

plot(NA, NA, xlim = range(times), ylim = range(c(result.diff.trans$SV, result.diff.trans$IV)), xlab = 'Days', ylab = '', main = 'Vector')
legend('right', legend = c('Susceptible', 'Infected'), col = c('blue', 'red'), pch = 15, bty = 'n')
lines(result.diff.trans$SV, col = 'blue', lwd = 2)
lines(result.diff.trans$IV, col = 'red', lwd = 2)

## TO DO: 
#### 1) How are these results different? What was the effect of making betaV smaller? 

#### 2) [write your own code:] Change betaH = 0.001, betaV = 0.0001 Plot these results

#### 3) Why does the number of susceptible individuals increase over time 
####    if there are no births in the model? 




##############################################################################
## Exercise: Different recovery rates
##############################################################################

## TO DO (see code below):
#### 1) [modify code]: Now re-run the model where transmission values are the same 
####                   (betaH = betaV) assuming that humans recover every 30 days.
####                   (Hint: w = 1/30)

#### 2) [modify code]: Now re-run the model where transmission values are different 
####                   (betaH > betaV) assuming that humans recover every 30 days. 

#### 3) How are these results different? the same? 

## set the different transmission parameters and run your model 
rec.30.same.trans.vector.parameters <- c() ## need to write your parameters out too

#run the model:
rec.30.result.same.trans<-data.frame(lsoda(y=initial.human.vector.pops, times = times, func = vector.trans.model, parms = rec.30.same.trans.vector.parameters))

## set the different transmission parameters and run your model 
rec.30.diff.trans.vector.parameters <- c() ## need to write your parameters out

#run the model:
rec.30.result.diff.trans<-data.frame(lsoda(y=initial.human.vector.pops, times = times, func = vector.trans.model, parms = rec.30.diff.trans.vector.parameters))

## plot your results
par(mfrow=c(1,2))
plot(NA, NA, xlim = range(times), ylim = range(c(rec.30.result.diff.trans$SH, rec.30.result.diff.trans$IH)), xlab = 'Days', ylab = '', main = 'Human')
legend('right', legend = c('Susceptible', 'Infected'), col = c('blue', 'red'), pch = 15, bty = 'n')
lines(rec.30.result.diff.trans$SH, col = 'blue', lwd = 2)
lines(rec.30.result.diff.trans$IH, col = 'red', lwd = 2)

plot(NA, NA, xlim = range(times), ylim = range(c(rec.30.result.diff.trans$SV, rec.30.result.diff.trans$IV)), xlab = 'Days', ylab = '', main = 'Vector')
legend('right', legend = c('Susceptible', 'Infected'), col = c('blue', 'red'), pch = 15, bty = 'n')
lines(rec.30.result.diff.trans$SV, col = 'blue', lwd = 2)
lines(rec.30.result.diff.trans$IV, col = 'red', lwd = 2)

##############################################################################
## Optional extra exercise: Different waning immunity rates
## (Don't do during the tutorial depending on time)
##############################################################################

## TO DO 
#### 1) [write your own code]:  Now re-run both the model where transmission values are the same (betaH = betaV) and assuming that immunity wanes after 10 days
#### 2) [write your own code]:  Now re-run both the model where transmission values are the same (betaH = betaV) and assuming that immunity wanes after 100 days
#### 3) [write your own code]:  Now re-run both the model where transmission values are the different (betaH > betaV) and assuming that immunity wanes after 10 days
#### 4) [write your own code]:  Now re-run both the model where transmission values are the different (betaH > betaV) and assuming that immunity wanes after 100 days

##############################################################################
## Temperature dependent mortality  
##############################################################################

### NOTE: The relationship between mortality and temperature is NOT based on published 
###       literature. Here, we just made a function to illustrate how you could 
###       incorporate a temperature dependent mortality rate -- but you should base 
###       your equation on experimental evidence (ideally)

## We will first generate some synthetic temperature data
test.times<-seq(0,135, length = 367)
temp.variables<-20+5*cos(test.times*0.05)
max.temp = max(temp.variables, na.rm=T)
min.temp = min(temp.variables, na.rm=T)
sigma.values<-seq(1/2,1/45,length = (max.temp - min.temp)+2)

## This is our fake temperature dependent mortality rate
sigma_func<-function(t){
  temp = temp.variables[t+1]
  return(sigma.values[ceiling(max.temp - temp)+1])}

### Let's see what synthetic temperature and mosquito mortality data we have generated
par(mfrow=c(1,2))
plot(temp.variables, xlab = 'Day', ylab = 'Temperature')
plot(sigma_func(seq(0,365*1, by = 1)), xlab = 'Day', ylab = 'Mortality Rate', main = 'Vector')

### Now we can change the original SIR, SI human-vector model to include a temperature 
### dependent mortality rate (sigma_func):

### Define the new model:
vector.trans.w.sigma<-function(t, y, parms){
  with(c(as.list(y), parms),{
    sigma = sigma_func(ceiling(t)) #changing sigma from a fixed value to a function of time
    
    dSHdt <- -betaH*SH*IV + w*RH
    dIHdt <- betaH*SH*IV - gamma*IH
    dRHdt <- gamma*IH - w*RH
    
    dSVdt <- mu - betaV*SV*IH - sigma*SV
    dIVdt <- betaV*SV*IH - sigma*IV
    list(c(dSHdt, dIHdt, dRHdt, dSVdt, dIVdt))
  })
}

# Set the initial populations and vector parameters:
initial.pops <- c(
  SH = 1000, 
  IH = 1, 
  RH = 0, 
  SV = 1e6, 
  IV = 100)
vector.parms.w.temp <- c(
  betaH = 0.00005, 
  betaV = 0.00005, 
  gamma = 1/14, 
  w = 1/60, 
  mu = 100000)

# Run the new model that includes temperature dependent vector mortality:
times<-seq(0,365*1, by = 1)
result.w.temp<-data.frame(lsoda(y = initial.pops, times = times, func = vector.trans.w.sigma, parms = vector.parms.w.temp))

# For comparison: Re-run the old model without temperature dependent vector mortality
vector.parms <- c(betaH = 0.00005, betaV = 0.00005, gamma = 1/14, w = 1/60, mu = 100000, sigma = mean(sigma_func(temp.variables)))
result<-data.frame(lsoda(y = initial.pops, times = times, func = vector.trans.model, parms = vector.parms))

# Let's look:
par(mfrow=c(2,2))
plot(NA, NA, xlim = range(times), ylim = range(c(result.w.temp$SH, result.w.temp$IH)), xlab = 'Days', ylab = '', main = 'Human.misy.temp.dep')
legend('right', legend = c('Susceptible', 'Infected'), col = c('blue', 'red'), pch = 15, bty = 'n')
lines(result.w.temp$SH, col = 'blue', lwd = 2)
lines(result.w.temp$IH, col = 'red', lwd = 2)

plot(NA, NA, xlim = range(times), ylim = range(c(result.w.temp$SV, result.w.temp$IV)), xlab = 'Days', ylab = '', main = 'Vector.misy.temp.dep')
lines(result.w.temp$SV, col = 'blue', lwd = 2)
lines(result.w.temp$IV, col = 'red', lwd = 2)

plot(NA, NA, xlim = range(times), ylim = range(c(result$SH, result$IH)), xlab = 'Days', ylab = '', main = 'Human.tsisy.temp.dep')
lines(result$SH, col = 'blue', lwd = 2)
lines(result$IH, col = 'red', lwd = 2)
plot(NA, NA, xlim = range(times), ylim = range(c(result$SV, result$IV)), xlab = 'Days', ylab = '', main = 'Vector.tsisy.temp.dep')
lines(result$SV, col = 'blue', lwd = 2)
lines(result$IV, col = 'red', lwd = 2)

## What do we see?
## Why is the effect on the human population less than on the vector population?
#### What are the starting populations of humans and mosquitoes? Compare mu and S!
#### When does the mosquito population increase

