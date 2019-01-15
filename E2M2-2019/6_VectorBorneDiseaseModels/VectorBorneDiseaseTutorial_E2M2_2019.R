#######################################
## Introduction 
#######################################
## Modeling vector-borne diseases R tutorial
## E2M2 - 2019
## written by: Amy Wesolowski


## In this tutorial, we will construct a simple vector-borne disease model (humans: SIRS, vector: SI). First, we will make a function that is just an SIRS model:

##############################################################################
## Contents: 
# Line 24: SIRS model
# Line 67: SIRS human model, SI vector model
# Line 133: Exercise: Different recovery rates
# Line 164: Different waning immunity rates
##############################################################################
  
## load your libraries first

library(deSolve)

##############################################################################
## SIRS Model
##############################################################################

sirs.model<-function(t, y, parms){
  with(c(as.list(y), parms),{
    dSdt <- -betaH*S*I + w*R
    dIdt <- betaH*S*I - gamma*I
    dRdt <- gamma*I - w*R
    list(c(dSdt, dIdt, dRdt))
  })}

## NOTES: This function as an additional parameter, w, which measures waning immunity: i.e. the rate at which recovered individuals loose their immunity and become susceptible again. Note: here we have just constructed just a basic SIRS model -- without a vector component. Next, we will investigate how changes in the waning immunity rate change the size of the outbreak. We will run 3 simulations using the same initial conditions and other parameters and only change w. We will run each simulation for 1 year and compare the results. 

## set all of your parameters and run your 3 models
initial.pops<-c(S = 1000, I = 1, R = 0)
parameters.1<-c(betaH = 0.0001, gamma = 1/14, w = 1/10)
parameters.2<-c(betaH = 0.0001, gamma = 1/14, w = 1/50)
parameters.3<-c(betaH = 0.0001, gamma = 1/14, w = 1/100)

times<-seq(0,365,by=1)
result.1<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.1))
result.2<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.2))
result.3<-data.frame(lsoda(y=initial.pops, times = times, func = sirs.model, parms = parameters.3))

## TO DO: 
#### 1) Now, we can plot these results to see the effect of changing the waning immunity between the three simulations. What differences do you expect? What would you expect if w = 1/2?
#### 2) What would you expect if w = 1/200? 

## plot your results
par(mfrow=c(1,2))
plot(NA, NA, xlim = c(range(times)), ylim = range(c(result.1$S, result.2$S, result.3$S)), xlab = 'Time', ylab = 'Susc')
lines(result.1$S, col = 'red')
lines(result.2$S, col = 'orange')
lines(result.3$S, col = 'blue')
legend('topright', legend = c('w=1/10', 'w=1/50', 'w=1/100'), col = c('red', 'orange', 'blue'), pch = 15, bty = 'n')

plot(NA, NA, xlim = c(range(times)), ylim = range(c(result.1$I, result.2$I, result.3$I)), xlab = 'Time', ylab = 'Infect')
lines(result.1$I, col = 'red')
lines(result.2$I, col = 'orange')
lines(result.3$I, col = 'blue')
legend('topleft', legend = c('w=1/10', 'w=1/50', 'w=1/100'), col = c('red', 'orange', 'blue'), pch = 15, bty = 'n')

##############################################################################
## SIRS human model, SI vector model
##############################################################################

## NOTES: Now we will use the same structure as an SIRS model for humans, but also add in an SI model for vectors. Note: now transmission for humans will depend on mosquitoes and visa versa!
  
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

## NOTES: This function has a number of additional parameters including: transmission for humans (betaH), waning immunity (w), recovery rate humans (gamma), vector birth rate (mu), transmission for vectors (betaV), and a death rate for vectors (sigma). Here, we are assuming that vectors do not recover! They can only be in either the susceptible or infected class. 

## Now, we can run a simulation first assuming that the transmission between humans and vectors and between vectors and humans is equal. 

## set your parameters and run your model
initial.human.vector.pops <- c(SH = 1000, IH = 1, RH = 0, SV = 1e6, IV = 0)
same.trans.vector.parameters <- c(betaH = 0.00005, betaV = 0.00005, gamma = 1/14, w = 1/60, mu = 100, sigma = 1/14)

result.same.trans<-data.frame(lsoda(y=initial.human.vector.pops, times = times, func = vector.trans.model, parms = same.trans.vector.parameters))

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

## Now, we will change the transmission variables and assume that the transmission rate between humans and vectors (betaV) is less than between vectors and humans (betaH)

## set the different transmission parameters (where betaV does not equal betaH) and run your model 
diff.trans.vector.parameters <- c(betaH = 0.00005, betaV = 0.00001, gamma = 1/14, w = 1/60, mu = 100, sigma = 1/14)

result.diff.trans<-data.frame(lsoda(y=initial.human.vector.pops, times = times, func = vector.trans.model, parms = diff.trans.vector.parameters))

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
#### 1) How are these results different? 
#### 2) [write your own code:] Change betaH = 0.001, betaV = 0.0001 Plot these results
#### 3) Why does the number of susceptible individuals increase over time if there are no births in the model? 

##############################################################################
## Exercise: Different recovery rates
##############################################################################

## TO DO 
#### 1) [modify code]: Now re-run the model where transmission values are the same (betaH = betaV) assuming that humans recover every 30 days. 
#### 2) [modify code]: Now re-run the model where transmission values are the different (betaH > betaV) assuming that humans recover every 30 days. 
#### 3) How are these results different? the same? 

## set the different transmission parameters and run your model 
rec.30.same.trans.vector.parameters <- c() ## need to write your parameters out too
rec.30.result.same.trans<-data.frame(lsoda(y=initial.human.vector.pops, times = times, func = vector.trans.model, parms = rec.30.same.trans.vector.parameters))

## set the different transmission parameters and run your model 
rec.30.diff.trans.vector.parameters <- c() ## need to write your parameters out

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
## Exercise: Different waning immunity rates
##############################################################################

## TO DO 
#### 1) [write your own code]: Now re-run both the model where transmission values are the same (betaH = betaV) and assuming that immunity wanes after 10 days
#### 2) [write your own code]: Now re-run both the model where transmission values are the same (betaH = betaV) and assuming that immunity wanes after 100 days
#### 3) [write your own code]: Now re-run both the model where transmission values are the different (betaH > betaV) and assuming that immunity wanes after 10 days
#### 4) [write your own code]: Now re-run both the model where transmission values are the different (betaH > betaV) and assuming that immunity wanes after 100 days


