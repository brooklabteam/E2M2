## Epidemic Cards Tutorial: computer exercise
## E2M2: Ecological and Epidemiological Modeling in Madagascar
## November 27-December 1, 2016

## Cara Brook, 2016

## By now, you've been exposed to some of the variety of model
## types used to understand infectious disease data. We'll now work
## with a few of those forms to explain our own data.

## By the end of this tutorial, you should be able to:
##
##  * Make a simple SIR model in discrete time.
##  * Visualize data from a time series of Infected and Susceptible case counts
##  * Fit a simple SIR model (i.e. estimate the transmission parameter) with 
##    logistic regression techniques.

rm(list=ls())
######################################################################
## Part 0: Visualizing your data

## Earlier, you played an epidemic card game and generated several time
## series of infecteds and susceptibles. Let's take a look at these data.

## Import .csv file 
dat <-  read.csv("Epidemic_Cards.csv", header=T)

## Take a look at the form of the dataset
head(dat)

## Take a subset of these data, choosing just those with R0 value = 2
dat.R0 <- subset(dat, R0==2)

## First we define the length of our time series. Recall that we played the
## game for 10 rounds, so we'll have a time series of 10 timesteps.

time <- seq(1,10, by=1)

## And we define N
N <- 26

## Plot Infecteds over time, with a different line for each trial
## First create a blank plotting frame
plot(1, type= "n",  xlim = c(1,10), ylim = c(0,26), xlab="timestep", ylab = "individuals")
## Now let's plot the time series of infected individuals for each trial
numb_trials = length(unique(dat.R0$trial))

for(i in 1:numb_trials){
  points(time, dat.R0$infecteds[dat.R0$trial == i], type = 'b', col = 'red')
}


## and do the same for susceptibles

for(s in 1:numb_trials){
  points(time, dat.R0$susceptibles[dat.R0$trial == s], type = 'b', col = 'green')
}

## Now let's add a legend to the plot
legend("topright", legend=c("susceptible", "infected"), col= c("green", "red"), lty=1, pch=16, bty='n')

## Note the use of inline function definitions and the with() function in
## the code above. Do you remember how these work? If not, you may want to
## revisit Tutorial 3.

## Now let's plot R_e across the time series of the epidemic.
## We could plot this individually for each of our time series, but instead, 
## let's take the average proportion susceptible for each point in 
## the time series. 

## First, let's add a column to the data frame that gives the value of R_e:

dat.R0$R_e <- dat.R0$R0 * dat.R0$susceptibles / N

## Next, we'll calculated the across trials for each timestep 

numb_timesteps<-length(unique(dat.R0$timestep))
mean.R_e<-rep(NA, numb_timesteps)
for(i in 1:numb_timesteps){
  mean.R_e[i] = mean(dat.R0$R_e[which(dat.R0$timestep == i)])
}

## Now we can plot the mean R_e through time. On average, when should the epidemic
## start to decline? Remember that we have already defined a vector of times (line 73).
plot(time, mean.R_e, type = "b",ylim=c(0,2), ylab=expression(R[e]))
abline (h=1, lty=2, col="red")

## Let's use a similar approach to calculate and plot the mean FOI through time
## from the simulated data:


## First calculate the FOI for each round in each trial
dat.R0$FOI <- dat.R0$R0 * dat.R0$infecteds / N

## Then average across trials
mean.FOI <- with(dat.R0,
                 tapply(FOI,timestep,mean)
)

## Then plot the result
plot(time, mean.FOI, type = "s", col="red", ylim=c(0,1), ylab="Mean FOI")

## How do things change when R0 = 3? Go back above, and take a 
## new subset of the data with R0 = 3, and repeat!

######################################################################
##Part 1: Modeling the epidemic

## Now let's build a model to recapture these data!

## We learned this week that models can take a variety of
## forms. The primary axes of model differentiation are: (1) discrete vs. 
## continuous treatment of time and (2) deterministic vs. stochastic model
## forms. Here, we model our epidemic from Epidemic Cards the game in discrete time.

## The version of the model we'll work with here has 3 state variables (susceptible,
## infected, and recovered), but since our population size is constant, we only 
## need to keep track of two of these, since we could always subtract to get the third.
## This is reason why we only kept track of susceptibles and infecteds
## in our game while counting cards.

## Our state variables are given as:

## S - the number of susceptibles in the population
## I - the number of infected/infectious individuals in the population

## and we have parameter, the basic reproduction number (R0). As we have already learned,
## R0 for a pathogen gives the expected number of new infectious individuals engendered by
## one infectious individual in a completely susceptible population.

## Also, remember that the efective reproduction number (R_e) is the expected number of new
## infectious individuals engendered by one infectious individual and can be calculated
## as R_e = R0 * S/N, an expression that scales R0 by the proportion of contacted individuals
## that will be susceptible to infection.

## Recall that cards are placed in the 'infection pile' and stay there for exactly 1 round,
## which means we can represent the number of infectious indviduals at a given time as a product
## of R_e and the number of infectious individuals in the previous timestep.

## We can now write down a model that describes the epidemic cards simulation
## using the following set of discrete time equations:

## S[t+1] = S[t] - R_e*I[t]
## I[t+1]= R_e*I[t]

## where S[t] indicates the population susceptible at time t and S[t+1]
## indiciates the population susceptible at one timestep in the future (with
## similar notation for I); or, equivalently:

## S[t+1] = S[t] - (R0 * S[t]/N)*I[t]
## I[t+1]= (R0 * S[t]/N)*I[t]

## Now we have all the information we need to make our discrete time model.

## First we define the length of our time series. Recall that we played the
## game for 10 rounds, so we'll have a time series of 10 timesteps.

time <- seq(1,10, by=1)

## Then, we make an empty vector of I and S for the full duration of the
## time series.

model.I <- rep(NA, length(time))
model.S <- rep(NA, length(time))

## And then we "seed" this vector with our initial conditions of the numbers 
## infected and susceptible. We start with one infected individual in a 
## population that is otherwise completely susceptible.

model.I[1] <- 1
model.S[1] <- 25

## We also need to specify a value of R0; for the first two rounsd of the simulation, we
## allowed each "infected" card to "infect" up to 2 other people; thus,

R0 <- 2

## Now we use the discrete time model specified above, to iterate forward in time:

for (t in 1:(length(time)-1)){
  model.S[t+1] = model.S[t] - (R0 * model.S[t]/N) * model.I[t]
  model.I[t+1] = (R0 * model.S[t]/N) * model.I[t]
}

par(bty="L",lwd=3,pch=16,mar=c(4,5,1,1),cex.axis=1.3,cex.lab=1.3)
plot(model.I, type="b", col="red", ylab="indivduals", xlab="time", ylim=c(0,26), xlim=c(0,10))
lines(model.S, col="green",type='b')
legend("topright", legend=c("susceptible", "infected"), col= c("green", "red"), lty=1, pch = 16, bty='n')

## You can see that the Susceptible population declines across the course of 
## the epidemic. What would happen if births were added to this population?


## Now let's calculate the effective reproduction nnumber throughout the epidemic.
## Remember that R_e[t] = R0 * S[t]/N
## An epidemic will increase when R_e > 1 and decrease when R_e < 1

model.R_e <- R0 * model.S/N

## We can plot R_e through time as well. Where should the epidemic start to decline?
plot(time, model.R_e, type = "b",ylim=c(0,2))
abline(h=1, lty=2, col="red")

## We might also like to look at the force of infection (FOI). The FOI
## helps us understand the trajectory of an epidemic from the perspective of
## infection risk to susceptible individuals (as opposed to R_e, which is from
## the perspective of infected individuals). Here, FOI[t] = R0 * I[t]/N. 

FOI <- R0 * model.I/N

## We can plot this over time, as well. Because we only "count" the 
## number infected at the beginning of each timestep, FOI should be constant 
## across the duration of each time step, then change based on the new 
## proportion in the next timestep (hence the "stairs"):
plot(time, FOI, type = "s", col="red", ylim=c(0,1))


## How do things change when R0 = 3? 
## Go back above, change R0 = 3, and repeat!

######################################################################
##Section 2: How likely are we to recover the observed data, assuming
## our model is true?

## We've now visualized our data and modeled our data, but how well does the 
## model recapture the data? Let's look. First re-run your model and select
## the data subset above, for which R0=2.

## Then, let's plot our model and data together.
## First plot your model lines, as above.

plot(time, model.I, type="b", col="red", ylab="indivduals", xlab="time", ylim=c(0,26), xlim=c(0,10))
lines(time, model.S, col="green",type='b')
## Then overlay your different trial time-series of infecteds, as above
numb_trials = length(unique(dat.R0$trial))

for(i in 1:numb_trials){
  lines(time, dat.R0$infecteds[dat.R0$trial == i], lty = 2, col = 'red', lwd = 0.75)
}
## and do the same for susceptibles
for(s in 1:numb_trials){
  lines(time, dat.R0$susceptibles[dat.R0$trial == s], lty = 2, col = 'green', lwd = 0.75)
}

## Now let's add a legend to the plot
legend("topright", legend=c("susceptible - modeled", "infected - modeled", "susceptible - data", "infected - data"), col= c("green", "red", "green", "red"), lty=c(1,1,2,2), lwd=c(3,3,1.5,1.5), pch=c(16,16,NA, NA), bty='n')

## We see that there is variation from trial-to-trial in our data and that,
## sometimes, our model fits the data better than other times. Why might this
## be? Did we get the same time series of S and I every time we played a new trial?
## What does this tell us about the stochastic vs. deterministic nature of 
## the real world? 

## Usually we would not be so lucky as to know our model parameters (in this case, there is just
## one: R0) ahead of time. So let's fit a model with flexible parameters to our
## data and see how well it does.

######################################################################
##Section 3: Optimize parameters to 'fit' your model to your data

## Now, let's assume that we had only our data and wanted to produce a model that 
## 'fit' these data well without knowing R0 in advance. There are many methods used to
## fit models to data, and while this diversity can at first be daunting - like many
## things, this becomes easier with experience. We've already talked a lot about different
## model fitting techniques, and here, we'll just demonstrate one simple method that will just minimize the sum of squared errors.

## From our equations above, we know that:
## I[t+1]= (R0 * S[t]/N)*I[t]

## Because we made the assumption that we have no overlapping generations of 'Infecteds,'
## we see that in every 'new' timestep, Infecteds_new = R0*Susceptibles_old*Infecteds_old
## where 'old' timesteps are one before the 'new.'

## So our equation looks like this:
## I_new = R0*S_old/N*I_old

## First, let's separate our data into two subsets: 'new' and 'old'. Since they
## will be regressed against one another, we can represent the 'new' subset 
## as timesteps t=2 through t=10 and the 'old' subset as timesteps t=1 through
## t=9. Once subset, we can form a new dataframe with the new and old values
## for I and the old values for S to enable the glmfit as above. Be sure to keep
## the column for 'trial' so as to allow for the random effect shown above.
## We'll also keep the timestep of the 'old' dataframe, so we'll always be tracking
## the transmission that follows each timestep.

dat.R0.fit <- cbind.data.frame(subset(dat.R0,timestep!=1)$trial,subset(dat.R0,timestep!=10)$timestep,subset(dat.R0,timestep!=1)$infecteds,subset(dat.R0,timestep!=10)$infecteds, subset(dat.R0,timestep!=10)$susceptibles)
names(dat.R0.fit) <- c("trial","timestep", "new_infecteds", "old_infecteds", "old_susceptibles")


## We can't calculate transmission when the infecteds get down to 0 because
## no transmission is taking place, so let's select the subset of the data for which 
## we had at least one infected and one susceptibles as predictors (an 'old' infected
## and an 'old' susceptible). It's okay if we get zeros in our response variable ('new' 
## infecteds) because this is still usable data.
dat.R0.fit <- subset(dat.R0.fit, old_infecteds>0 & old_susceptibles>0)

## Remember that from our equations above, we know that:
## I[t+1]= (R0 * S[t]/N)*I[t]


## And now we want to estimate a new R0 value. We'll first need to construct
## an array of possible R0 values and then use the sum of squared errors to 
## identify which one is most likely

possible.R0<-seq(0,10,length=100)
sum.sq.error.values<-rep(NA, length=100)
for(i in 1:100){
  diff.btwn.data.and.est<-(dat.R0.fit$new_infecteds - ((possible.R0[i]*dat.R0.fit$old_susceptibles/N)*dat.R0.fit$old_infecteds))^2
  diff<-sum(diff.btwn.data.and.est)
  sum.sq.error.values[i] = diff
}

plot(possible.R0, sum.sq.error.values, xlab = 'possible R0', ylab = 'SSE')
## and identify the best fit R0
best.fit.R0<-possible.R0[which(sum.sq.error.values == min(sum.sq.error.values))]
print(best.fit.R0)


## we can then use our best.fit.R0 to see how well we did
model.I.new <- rep(NA, length(time))
model.S.new <- rep(NA, length(time))
model.I.new[1] <- 1
model.S.new[1] <- 25

## Our new value for R0 and our previous value for N are already specified, so 
## we can now iterate forward in time.

## Now we use the discrete time model specified above, to iterate forward in time:

for (t in 1:(length(time)-1)){
  model.S.new[t+1] = model.S.new[t] - (best.fit.R0 * model.S.new[t]/N) * model.I.new[t]
  model.I.new[t+1] = (best.fit.R0 * model.S.new[t]/N) * model.I.new[t]
}


## Now we plot our 'fitted' model with our data.
plot(time, model.I.new, type="b", col="red", ylab="indivduals", xlab="time", ylim=c(0,26), xlim=c(0,10))
lines(time, model.S.new, col="green",type='b')
## Time-series of infecteds (data), as above
numb_trials = length(unique(dat.R0$trial))

for(i in 1:numb_trials){
  lines(time, dat.R0$infecteds[dat.R0$trial == i], lty = 2, col = 'red', lwd = 0.75)}

## and do the same for susceptibles
for(s in 1:numb_trials){
  lines(time, dat.R0$susceptibles[dat.R0$trial == s], lty = 2, col = 'green', lwd = 0.75)}

## And our edited legend
legend("topright", legend=c("susceptible - fitted model", "infected - fitted model", "susceptible - data", "infected - data"), col= c("green", "red", "green", "red"), lty=c(1,1,2,2), lwd=c(3,3,1.5,1.5), pch=c(16,16,NA, NA), bty='n')

## How well did we do at recapturing our data, as compared with our 'true model', plotted above?


## What does this tell us about the model fit?

## Try again, this time fitting the subset of the data for which R0 = 3. How well do you do at fitting the
## correct value here? Is it closer or not as close to your true R0 value as when you fit the subset for 
## R0 = 2? Why might this be?






