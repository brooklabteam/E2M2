## Tutorial: Fitting and Comparing Age-Prevalence Models
## E2M2: Ecological and Epidemiological Modeling in Madagascar
## 6-14 January, 2020

## Cara Brook, 2020
## with code and data from Brook et al. 2017. Epidemics.

## This tutorial expands on the previous Introduction to Model Fitting tutorial.
## You will again use maximum likelihood techniques to fit models to data. Here,
## you will fit to dynamic age-prevalence data of diverse strains of rodent-borne
## Bartonella spp., collected from Rattus rattus in Madagascar.

## You will use model comparison via likelihood ratio test (LRT) and AIC to evaluate 
## the fit of two different model forms against each strain type and decide on the 
## optimal model structure to fit each type of data.


rm(list=ls())
######################################################################
## Part Zero: Visualizing your data and forming hypotheses

## First load your data. Data have already been transformed into age-prevalence form 
## from age-weight, via the von Bertallanffy equation. Age is in days, and associated
## metadata is attached. We provide two datasets of the same length. The same rats were
## assayed for two different strains of Bartonella.

load("age-prev-rat-bart.Rdata")

## Look at your data
head(dat.eliza)
head(dat.pho)

## Now plot the prevalence (i.e. proportion infected) per age class
## for each strain type. 

## Will need to load a few packages first.
require(plyr)
require(ggplot2)

## First, use plyr to summarize prevalence by age class for each dataset
dat.eliza.sum <- ddply(dat.eliza, .(age.class), summarize, prev=sum(sick)/length(sick), count = length(sick))
## Make sure that any NAs are listed as 0
dat.eliza.sum$prev[is.na(dat.eliza.sum$prev)] <- 0
## and add the strain identifier back in
dat.eliza.sum$strain <- "elizabethae"

dat.pho.sum <- ddply(dat.pho, .(age.class), summarize, prev=sum(sick)/length(sick), count = length(sick))
## Make sure that any NAs are listed as 0
dat.pho.sum$prev[is.na(dat.pho.sum$prev)] <- 0
## and add the strain identifier back in
dat.pho.sum$strain <- "phoceensis"

## Join them together
dat.sum <- rbind(dat.pho.sum, dat.eliza.sum)

## Attach the corresponding midpoint to each age class for plotting
dat.sum$midpt <- rep(c(7.5, 22.5, 37.5, 52.5, 67.5, 142.5), times=2)

## And feed to ggplot. Write over age.class on the x-axis, so that we see the actual raw
## values by age and plot the cumulative prevalence for that bin at the midpoint of the bin.

p <- ggplot(data=dat.sum) + 
  geom_point(aes(x=midpt, y=prev, size=count, fill=strain), pch=21) +
  xlim(0, 200) + theme_bw() + xlab("Age (Days)") + ylab("Bartonella spp. Prevalence")
print(p)

## Do the age-prevalence patterns for the two strains look the same or different?
## What do you notice about each strain's trajectory across ages? 

######################################################################
## Part One: Modeling Your Data

## Previously, we've fit models to time series of quantitites (i.e. savanna/fire,
## Epidemic Cards), but this time, our data are different. What are the states
## in this system? How should our new model be structured?

## Similarly! We have age-prevalence data now. Bartonella spp. are traditionally
## thought to be persistent infections, so we might expect that prevalence increases
## monotonically with age.

## In equation form, this would look like:
## (dS(a))/da=(1-I(a)) - λ(a)(1-I(a))
## (dI(a))/da=λ(a)(1-I(a))

## This is an S-I model, so we have no recovered class, allowing us to focus on
## one of the state variables only. Since our data give us the proportion infected
## by age, we'll model the I state variable, so the second equation.

## We can then do some calculus to solve the (dI(a))/da equation for I. Since we
## assume that all rats are born susceptible (so I(0) = 0), we can reduce the math
## to the following expression:

## I(a) = 1- exp[-inte(λ(a)da)]

## In order to model this, we make the age-specific Force of Infection (FOI) 
## piece-wise across age classes. This is sort of like adopting a discrete time
## model, and it simplifies our model fitting considerably. We can then model
## the age specific prevalence and write a corresponding likelihood function
## to ask how likely our data are, given our model.

## Let's start build our model from above and plot over the data. We'll provide guesses 
## for the lambda values for each age bin and then see how well our model matches our data.

## Each line of the model is annotated to help you understand what is going on
model.age.prev <- function(log.lambda, data, cate){ 
  #takes in a value for lambda from each age bin, the data itself, and a vector that gives the lower bound of every age bin
  dur=c(diff(cate), 0) #calculates the length of each age bin
  p=as.list(rep(NA,max(data$age))) #make a list to store the prevalence by age predicted from the model
  for(a in 1:max(data$age)){ #for-loops over all possible ages in our data range
    dummy1=a>cate  #creates dummy variables to identify which age class this individual falls in
    dummy2 = a>cate & !c(a>cate[-1], FALSE) #adds a false at the end and takes the inverse to make true only those values above you
    dummy1=c(a>cate, FALSE)[-1]# and dummy 1 gives you a true for all age bins below you
    inte=sum(dur*exp(as.numeric(log.lambda))*dummy1) + exp(as.numeric(log.lambda[dummy2]))*(a-cate[dummy2]) 
    #integrates over the amount of time in previous age classes + how long you've been in this one as a susceptible
    #to give you the hazard of infection
    p[[a]]=1-exp(-inte) #returns the probability of infection for a given age.
  }
  return(p) #returns prevalence by age for the length of all possible ages in our dataset
}

## Run with guesses and plot versus age on the same plot as before.
## Even though we grouped data above into 6 age classes, let's try just 3 for the 
## age-specific FOI in our model and see if we can still recover the same patterns we
## see in the data.

#then we get prevalence to plot with our other function, feeding it these values for log.lambda
prev <- model.age.prev(log.lambda=log(c(1e-02, 1e-02, 1e-02)), data=dat.pho, cate=c(0,15,30))
prev <- model.age.prev(log.lambda=log(c(.9e-02, .9e-02, .9e-02)), data=dat.pho, cate=c(0,15,30))
prev <- model.age.prev(log.lambda=log(c(.9e-02, .9e-02, .5e-02)), data=dat.pho, cate=c(0,15,30))
prev <- c(unlist(prev))

#then you make your model outputs and CIs into classes ggplot can handle that match above
age <- seq(1,max(dat.eliza$age))
dat.age.prev <- cbind.data.frame(age,prev)

#Now plot with data from above.
p + geom_line(data=dat.age.prev[1:200,], aes(age,prev))

## Do we need to run this model again with the B. phoceensis data? Why or why not?

## How well does this model fit the B. elizabethae data? The B. phoceensis
## data? Try varying the parameters to get a better fit. 
## How might we find the best fit of all?

######################################################################
## Part Two: How likely are your data, given your model?

## Here is our likelihood function. Each line of code is annotated to explain what
## is going on:

#first  write our likelihood function
loglikpc=function(par,data, cate){
  dur=c(diff(cate), 0) #calculates the length of each age bin
  ll=0 #sets inital log-likelihood to 0
  for(a in 1:length(data$age)){ #for-loops over every individual in the dataset
    dummy1=data$age[a]>cate #creates dummy variables to identify which age class this individual falls in
    dummy2 = data$age[a]>cate & !c(data$age[a]>cate[-1], FALSE) #adds a false at the end and takes the inverse to make true only those values above you
    #so dummy 2 gives you a true for the age bin you are in (meaning the lower bound of the age bin you are in)
    dummy1=c(data$age[a]>cate, FALSE)[-1]
    # and dummy 1 gives you a true for all age bins below you
    inte=sum(dur*exp(as.numeric(par))*dummy1) + exp(as.numeric(par[dummy2]))*(data$age[a]-cate[dummy2]) 
    #integrates over the amount of time in previous age classes + how long you've been in this one as a susceptible
    #to give you the hazard of infection
    p=1-exp(-inte) #and here we calculate the probability of becoming infected at your given age
    #which is essentially is the prevalance at a given age
    ll=ll+dbinom(data$sick[a],1,p,log=T) #now we assess for each individual: ask how likely is my infection status (1 = pos, 0=neg), given the 'true' probability of becoming infected at my age
  }
  return(-ll)
}

## Run this function for the above model and parameters on both datasets.
loglikpc(par = log(c(1e-02, 1e-02, 1e-02)), data = dat.eliza, cate = c(0,15,30))

## And do the same for the B. phoceensis data
loglikpc(par = log(c(1e-02, 1e-02, 1e-02)), data = dat.pho, cate = c(0,15,30))

## Which dataset did the model fit better? How can you tell?

## How might we find the best parameter fit of all?

######################################################################
## Part Three: Optimize your parameters to make your model fit your data

## Our model is structured to allow different values of lambda for each
## age class, so it's a bit silly to force them to be the same value 
## in each class. Let's use optim to find the best combination of parameters
## to fit this model to each set of data.

out.eliza <- optim(par=log(c(1e-02, 1e-02, 1e-02)),fn=loglikpc,
                   cate=c(0,15,30), method="BFGS", data=dat.eliza)


out.pho <- optim(par=log(c(1e-02, 1e-02, 1e-02)),fn=loglikpc,
                   cate=c(0,15,30), method="BFGS", data=dat.pho)

## Did you recover the same values of lambda for B. elizabethae and
## B. phoceensis? Why or why not?

## Now run your model again with these new optimized parameters and plot 
## the data from above.

##B. elizabethae
prev <- model.age.prev(log.lambda=out.eliza$par, data=dat.eliza, cate=c(0,15,30))
prev <- c(unlist(prev))

#then you make your model outputs and CIs into classes ggplot can handle that match above
age <- seq(1,max(dat.eliza$age))
fit.eliza.age.prev <- cbind.data.frame(age,prev)
fit.eliza.age.prev$strain <- "elizabethae"

## And B. phoceensis too
prev <- model.age.prev(log.lambda=out.pho$par, data=dat.pho, cate=c(0,15,30))
prev <- c(unlist(prev))

## Then you make your model outputs and CIs into classes ggplot can handle that match above
age <- seq(1,max(dat.pho$age))
fit.pho.age.prev <- cbind.data.frame(age,prev)
fit.pho.age.prev$strain <- "phoceensis"

fit.age.prev <- rbind(fit.eliza.age.prev,fit.pho.age.prev)

## Now plot with data from above.
fit.sub <- subset(fit.age.prev, age <= 200)
p + geom_line(data= fit.sub, aes(age,prev, colour=strain))

## How well does each model fit each data set?
loglikpc(par = out.pho$par, data = dat.pho, cate = c(0,15,30))

loglikpc(par = out.eliza$par, data = dat.eliza, cate = c(0,15,30))

## What might make these models fit even better?

######################################################################
## Part Four: Reconstruct Your Model to Better Match Your Data

## Notice that the optimized lambda values for the B. phoceensis
## model are all quite similar. Do we really need the added complexity
## of an age-varying lambda, or could we recapture the same dynamics
## with a constant value for lambda?

## Try it, optimizing with just one age class and compare to the fit
## in the model above.

out.pho2 <- optim(par=log(1e-02),fn=loglikpc,
                 cate=c(0), method="BFGS", data=dat.pho)

prev2 <- model.age.prev(log.lambda=out.pho2$par, data=dat.pho, cate=c(0))
prev2 <- c(unlist(prev2))

## Plot the two fits together
age <- seq(1,max(dat.pho$age))
fit.pho.age.prev$model <- "3 age class"
fit.pho.age.prev2  <- cbind.data.frame(age,prev2)
names(fit.pho.age.prev2) <- c("age", "prev")
fit.pho.age.prev2$strain <- "phoceensis"
fit.pho.age.prev2$model <- "1 age class"

pho.prev <- rbind(fit.pho.age.prev, fit.pho.age.prev2)

p + geom_line(data= pho.prev, aes(age,prev, linetype=model), colour="turquoise")

## It's hard to tell which line better represent the data.
## Let's compute the likelihood for each model and compare via AIC

## Calculate AIC
## AIC is computed as : -2LL+2k where LL is log-likelihood and k is the number of fitted parameters

## First comput the log-likelihood of each model fitted to the data
lik_m1 <- loglikpc(par = out.pho$par, data = dat.pho, cate = c(0,15,30))
lik_m2 <- loglikpc(par = out.pho2$par, data = dat.pho, cate = c(0)) 

AIC_m1 <- 2*lik_m1 + 2*length(out.pho$par)
AIC_m2 <- 2*lik_m2 + 2*length(out.pho2$par)


## Look at these two calculations of AIC:
AIC_m1 #155.4076
AIC_m2 #152.5627

## We want to get the smallest AIC to have the model that most closely fits the data.

## The model with 3 parameters (m1) has a bigger AIC than the simplest model (m1) with only one 
## age class, so we should just use the simpler, 1-age-class model. 

## If these were standard fitted models via regression or gam like those in 
## the epidemic card game, we could also compare via AIC.

############################

## What about B. elizabethae? Do we need all the age classes there?
## As before, run both models, plot, and compare.

out.eliza2 <- optim(par=log(1e-02),fn=loglikpc,
                  cate=c(0), method="BFGS", data=dat.eliza)

prev3 <- model.age.prev(log.lambda=out.eliza2$par, data=dat.eliza, cate=c(0))
prev3 <- c(unlist(prev3))

## Plot the two fits together
age <- seq(1,max(dat.eliza$age))
fit.eliza.age.prev$model <- "3 age class"
fit.eliza.age.prev2  <- cbind.data.frame(age,prev2)
names(fit.eliza.age.prev2) <- c("age", "prev")
fit.eliza.age.prev2$strain <- "elizabethae"
fit.eliza.age.prev2$model <- "1 age class"

eliza.prev <- rbind(fit.eliza.age.prev, fit.eliza.age.prev2)

p + geom_line(data= eliza.prev, aes(age,prev, linetype=model), colour="coral")

## This time, these lines are very different! What if we compare them formally?
## Let's compute AIC: -2LL+2k where LL is log-likelihood and k is the number of fitted parameters 

lik_m3 <- loglikpc(par = out.eliza$par, data = dat.eliza, cate = c(0,15,30))
lik_m4 <- loglikpc(par = out.eliza2$par, data = dat.eliza, cate = c(0)) 

k_m3 = length(out.eliza$par)
k_m4 = length(out.eliza2$par)

AIC_m3 <- 2*lik_m3 + 2*k_m3
AIC_m4 <- 2*lik_m4 + 2*k_m4

## Look at AIC:
AIC_m3 #114.6327
AIC_m4 #126.2653

## This time, the model with more parameters (m3) has a lower AIC than the simple model (m4)
## with only 1 age class. Typically, we think that changes in AIC greater than 2 points are 
## significant. This change of 20 points suggests that the added complexity
## of 3 ages classes creates a model with a significantly better fit to
## the data than that without. We definitely need all three age classes here.

## The model still fails to recapture the decline in prevalence at higher 
## ages. Might there be another way that we could adjust for that? How might
## we make prevalence actually decrease?

## We need a fundamentally different model! Instead of S-I, what if we tried
## to model S-I-S? The equations would be altered slightly to produce the 
## following instead:

## (dI(a))/da=λ(a)(1-I(a))-σ*I(a)

## where σ is recovery from "Infected" that returns you back to the 
## Susceptible class.

## When we solve for I, integrate, and make the assumption that S(0) = 0,
## we get the following:
## I(a)=(λ(a))/(λ(a)+ σ) (1-exp[-inte(λ(a)+  σ)da])

## Our piece-wise model then takes a slightly different form than above
## and includes a parameter for clearance:

model.age.prev.clear <- function(log.lambda, sigma, data, cate){ #takes in the optimized estimates for lambda from the dataset
  dur=c(diff(cate), 0)
  p=as.list(rep(NA,max(data$age))) #make a list to store the prevalenc by age predicted from the model
  for(a in 1:max(data$age)){ 
    dummy1=a>cate 
    dummy2 = a>cate & !c(a>cate[-1], FALSE) 
    dummy1=c(a>cate, FALSE)[-1]
    inte=sum(dur*exp(as.numeric(log.lambda))*dummy1) + exp(as.numeric(log.lambda[dummy2]))*(a-cate[dummy2]) + sum(dur*sigma*dummy1) + sigma*(a-cate[dummy2])
    discount_lambda = exp(as.numeric(log.lambda[dummy2]))/(exp(as.numeric(log.lambda[dummy2]))+sigma)
    p[[a]]=discount_lambda*(1-exp(-inte)) 
  }
  return(p) #returns prevalence by age for the length of the ages
} 

loglikpc.clear=function(par, data, cate){
  dur=c(diff(cate), 0) #calculates the length of each age bin
  ll=0 #sets inital log-likelihood to 0 
  for(a in 1:length(data$age)){ #for-loops over every individual in the dataset length(data$age)
    dummy1=data$age[a]>cate #creates dummy variables to identify which age class this individual falls in
    dummy2 = data$age[a]>cate & !c(data$age[a]>cate[-1], FALSE) #adds a false at the end and takes the inverse to make true only those values above you
    #so dummy 2 gives you a true for the age bin you are in (meaning the lower bound of the age bin you are in)
    dummy1=c(data$age[a]>cate, FALSE)[-1]
    # and dummy 1 gives you a true for all age bins below you
    inte=sum(dur*exp(as.numeric(par[1:(length(par)-1)]))*dummy1) + exp(as.numeric(par[dummy2][1]))*(data$age[a]-cate[dummy2]) + sum(dur*par[length(par)]*dummy1) + (par[length(par)])*(data$age[a]-cate[dummy2])
    #integrates over the amount of time in previous age classes + how long you've been in this one as a susceptible (need to specify the first value of par true for dummy 2 b/c otherwise it will take the recovery rate)
    #the integration is performed twice, once for lambda and once for sigma, based on integration rules. note that sigma is a scalar while lambda (in this case) is a vector but the second one has a constant recovery rate, so the dummy variables need to be separate from it...
    #to give you the hazard of infection
    discount_lambda = exp(as.numeric(par[dummy2][1]))/(exp(as.numeric(par[dummy2][1]))+par[length(par)])
    p=discount_lambda*(1-exp(-inte)) #and here we calculate the probability of becoming infected at your given age
    #which is essentially is the prevalance at a given age
    ll=ll+dbinom(data$sick[a],1,p,log=T) #now we assess for each individual: ask how likely is my infection status (1 = pos, 0=neg), given the 'true' probability of becoming infected at my age
    #why do we add the ll??? 
  }
  return(-ll)
}

## There is a corresponding likelihood function which we can use to 
## simultaneously optimize lambdas and sigma.

## Try it with your data and plot, as above.
out.eliza3 <- optim(par=c(log(c(4e-15,4e-15,4e-15)), .005),
                    fn=loglikpc.clear, cate = c(0,15,30),
                    method="Nelder-Mead", data=dat.eliza)
            
prev5 <- model.age.prev.clear(log.lambda = out.eliza3$par[1:3], sigma=out.eliza3$par[4],
                           data=dat.eliza, cate=c(0,15,30))

prev5 <- c(unlist(prev5))

## Plot the two fits together
age <- seq(1,max(dat.eliza$age))
fit.eliza.age.prev5  <- cbind.data.frame(age,prev5)
names(fit.eliza.age.prev5) <- c("age", "prev")
fit.eliza.age.prev5$strain <- "elizabethae"
fit.eliza.age.prev5$model <- "3 age classes with clearance"

eliza.prev <- rbind(eliza.prev, fit.eliza.age.prev5)

p + geom_line(data= eliza.prev, aes(age,prev, linetype=model), colour="coral")

## This actually allows our prevalence to drop back down!

## Compare the AICs of a=our previous model with 3 age classes and no clearance (m3) 
## with this model that has three age classes as well as a clearance term
lik_m5 <- loglikpc.clear(par = out.eliza3$par, data = dat.eliza, cate = c(0,15,30))
k_m5 <- length(out.eliza3$par) # this time, extra parameters are fitted

AIC_m5 = 2*lik_m5 + 2*k_m5

##Compare with m3 AIC:
AIC_m3 #114.6327
AIC_m5 #113.1223

## The model including clearance (m5) has a lower AIC by ~1 points, meaning that it 
## a marginally  better fit to your data than model 3.

## You've just fitted, evaluated, compared, and selected models from data.
## But, remember that you had to have the model in mind when collecting these data.
## What if you hadn't collected information on age?!? None of this would have been possible.
## Awareness of modeling outcomes in the process of study design is what we called
## 'Model Guided Field Work.'

## Challenge: How would you attempt to fit a unique FOI for each site and inside vs. outside
## locality, as is done in the paper?

## Is it possible for us to test an SIR model, too? Why or why not?





