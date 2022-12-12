#######################################
## Introduction 
#######################################
## Compartmental models and differential equation R tutorial
## Live version
## E2M2 - 2020
## written by: Jessica Metcalf, Amy Wesolowski, Cara Brook, Sophia Horigan

## In this tutorial, we will construct simple compartmental models and explore different structures of these models
##############################################################################
## Contents: 
# Modeling Population Growth
# Investigating continuous vs discrete time
# Structured population models
# Lotka-Volterra models of predator-prey dynamics
# Susceptible-Infected-Recovered models
##############################################################################
## load your libraries first

library(deSolve)

##############################################################################
## Modeling Population Growth - LIVE
##############################################################################

## The simplest formulation for population growth in discrete time is N_t+1=lambda N_t
## where N_t is the size of the population at time t. The World Bank provides time-series 
## of population size for all countries in the world, including for Madagascar; this is in 
## the file "WorldBankPop.csv". Keep in mind that this is not data from direct observations 
## - censuses have not occurred in all the years for which data is available! 
  
## First, bring in this data by identifying the right path to the folder you are working in
## (use the command 'getwd' to find out where you are, and then use 'setwd' to get R to look 
## in the right place) and have a look at its structure in R (just the first 10 columns): 

##Load your data from the world bank and peek at the first few rows
## 1. read in the population data


## 2. look at the first few rows of the data


## 3. Make a sub-vector called 'mada.pop' that is the row of your dataset which corresponds to Madagascar.

## 4. Change the cells with numbers to 'numeric' 


## 5. Make a second sub-vector called 'france.pop' that corresponds to France.

## 6. Change the cells with numbers to 'numeric'


## Now, we can calculate the growth rate, N_t+1/N_t, by taking a vector that goes from the
## 2nd until the last value in our population vectors (this is years 1961 until 2016, and 
## thus we call it N_t+1). We will divide this by the the vector that goes from the first 
## until the before-last column, which will correspond to 1960 to 2015. Relative to the 
## first vector, this is N_t - i.e., the first number in the N_t+1 vector is 1961; and the 
## first number in the N_t vector is 1960, and so on.  

## 7. First make a vector of Nt Mada

## 8. Next, make a vector of Nt+1 for Mada


## 9. Make a vector of Nt for France

## 10. Make a vector of Nt+1 for France


## 11. Now estimate population growth rate for Madagascar (what is the equation for lambda?)

## 12. do the same for France


## ??? Is the population growing or declining?


## growing! above 1.

## 13. Now plot (a) the population sizes for Madagascar and France from 1960-2016,
## side-by-side with the (b) population growth rates
## 13.1 Set up plot window to hold two graphs

## 13.2 plot Mada pop data as a solid line


## 13.3 add France data as dotted line

## 13.4 add legend


## 13.5 plot Mada pop growth rate as solid line


## 13.6 add France pop growth rate as a dotted line



## ?? The growth rates for Madagascar and France never drop below 1 (although France comes close). 
## What does this tell us? Which country's population is growing more quickly relative to its population size?

# The population is always increasing. sometimes it is increasing more slowly, sometimes more quickly. 
# Madagascar is increasing more quickly.


## TO DO: 
#### 1) Change the colour of the plots and the legend; 
#### 2) Experiment with plotting out different countries (Germany? Tanzania?)
#### 3) [maybe later] See if you can find which country in the entire data-base 
####    has the fastest growth, and when that occurred. You might want to use 
####    the 'which' command, with the argument 'arr.ind=TRUE' since this will 
####    index the matrix, i.e., tell you the row (which corresponds to the country) 
####    and the column (which corresponds to the year). Note that you will need 
####    to apply this to a matrix of growth rates!



##############################################################################
## Investigating continuous vs discrete time 
##############################################################################

## Is the model above structured like a continuous or a discrete time model?

## Let's make a discrete time population model. We will use a for-loop to project a 
## population 100 time-steps, starting with population size of 10, and with a population 
## growth rate equal to 1.1

## Set up your (a) vector "N" of future population sizes, (b) fill in your initial population size,
## and (c) define "lambda", the population growth rate.
N <- rep(NA,100)
N[1] <- 10
lambda=1.1


## Now make a discrete time model using a for-loop to iterate this population 100 years
## into the future.

for(t in 2:100){
  N[t] <- lambda*N[t-1]
}


## 1. Does anyone remember the equation to shorten this discrete time model? 
## If we do this, do we get the same result?



## Plot the results together.
par(mfrow=c(1,1))
plot(x=seq(1,100,1), y=N, pch=15, col="blue", ylab="N", xlab="time")
points(x=seq(1,100,1), y=Ndirect, pch=19, col="red", cex=0.5)

## What if we want to make a continuous time model? How do we write that?
## First, write the differential equation dN(t)/dt = rN in a function
## called 'ngrow'

ngrow <- function(t,y,parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    dNdt <- r*N
    list(c(dNdt))
  })
}


## To solve this, we're going to need to use the R library called "deSolve".
## We also need to define the starting variables and parameters: 

library(deSolve)              # Load libary to be used for numerical integration
Pop <- c(N=10)                # Define starting variables (P(0)=10, as above)
values <- c(r=log(lambda))    # Define parameters - use an approx of the lambda used in the discrete time ex above
time.out <- seq(0,100,by=0.1) # Set up time-steps for variables - here go for 100 years
ngrow(t=0,y=Pop,parms=values)  #try out the function for t=0

## The function outputs the time derivative of N at the time t. 
# This means that in order to get the (approximate) values of N at some time 
# delta.t in the future, we need to calculate the equivalent of N(t+delta.t)=N(t)+rN(t)*delta.t. 
# We can do this by hand: 

delta.t <- 0.1                    # Set a small value for delta.t (0.1 day)
pop.next <- Pop + unlist(ngrow(t=time,y=Pop,parms=values)) * delta.t
pop.next

## We could iterate this process, updating our values of Pop 
## and time with each iteration to get a discrete-time approximation 
## of the values of $P$ through time. However, differential equations describe 
## the rates of change. It turns out that using the discrete time approximation 
## of this process leads to rapid accumulation of error in our estimate of the state 
## variables through time, even if we set our value of delta.t to be very small. 
## Fortunately,  a number of algorithms have been developed to come up with better 
## (both more accurate and more computationally efficient) approximations to the values of the 
## state variables through time. One such is "lsoda" which we loaded above, and which we can run, 
## here storing the result in a variable called 'Ncontinuous': 


Ncontinuous <- lsoda(
y = Pop,               # Initial conditions for population
times = time.out,      # Timepoints for evaluation
func = ngrow,          # Function to evaluate
parms = values         # Vector of parameters
)

## We can check the first few lines of the output:
head(Ncontinuous)

## The data frame has 2 columns showing the values of time (labeled "time") 
## and the state variable ("N", which is population size) through time. Let's see what 
## has happened after 2 years, by printing it out:
subset(Ncontinuous,time==2)

## Now we can also plot the full time-series: 
plot(Ncontinuous[,"time"],               # Time on the x axis
     Ncontinuous[,"N"],                  # Number people on the y axis
     xlab = "Time in days",     # Label the x axis
     ylab = "Pop size",         # Label the y axis
     main = " ",                # Plot title
     xlim = c(0,100),           # 
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot
points(1:100,N, pch=19,col="red",cex=0.2)  #compare with our discrete time

# Note that we overlaid our discrete time model predictions (the variable N) for comparison.
## This is close, but not an exact match to our discrete time version, which makes sense, 
## since these are different processes.  


## TO DO: 
#### 1) See what happens to a population where r<0
#### 2) Change the time resolution (set by "time.out"") 
####   and see if/where the results change

                           
##############################################################################
## Structured population models - LIVE
##############################################################################

## So far, we have ignored structure in the population, and treated all individuals as if they were identical.
## We know this is not true. The simplest extension is to assume some form of stage structure, i.e., different states within the population. 
## Let us assume that the population has two stage classes, 'juveniles' and 'adults'. 
## We can describe the population using a vector of numbers of the two types of individuals n_t=(n_j,n_a) 
## where the j index indicates juveniles, and the a index indicates adults, n indicates 'numbers'. 

## 1. What four variables are in a structured population model? Assign a value to each.



## These parameters can be estimated by tracking individuals in the field and
## counting the number that survive between censuses (here, distinguishing 'juvenile' 
## and 'adult' stages), what number moved between stages, and the number of offspring 
## produced at each stage. From this, we obtain survival and aging probabilities, and 
## fertility. We place these parameters in a matrix which has as many rows and columns as
## we have stages. The columns define the stage that individuals start in, and the rows the 
## stage that individuals end up in: 

## 2. Create a 2x2 matrix of 0's


## 3. assign probabilities to the matrix


 

## 4. Take a look at the matrix 'A' in your R console:


## We can project the population forward by one time-step:  

## 5. create a list called n.start with a population of 5 juveniles and 10 adults


## 6. get pop one time-step later using matrix multiplication (%*%)


## 7. Look at the output in your R console



## We could write a loop, and project the population into the future (by applying matrix multipicaton to A and n.next, and so on), 
## and from that, figure out the growth rate of this population.
## There is  also a shortcut, using some core results from matrix population model theory (Caswell, 2001): 
## the population growth rate is defined by the dominant eigenvalue of the matrix. :

## 8. extract the population growth rate.This eigenvalue is the same as lambda

## 9. double the fertility


## 10. See how this changes the eigenvalue



## There are many more powerful results that yield estimates of the stable population structure, the reproductive value, 
## life expectancy, sensitivies, elasticities, etc, that have been used extensively in fields from conservation biology to 
## life history evolution. You can also implement density dependence, etc. There are some nice R packages that can help you do this 
## (e.g., 'popbio').



##############################################################################
## Lotka-Volterra models of predator-prey dynamics
##############################################################################

## Moving beyond considering just one species, a classic dynamic is that of predator-prey cycles.
## To caricature the process, we can imagine that fossa eat lemurs, until there aren't enough lemurs. 
## Then the fossa population crashes, which allows the lemur population to grow again, etc. 
## We can frame this in continuous time with the following equation:

LotVmod <- function (Time, State, Pars) {
with(as.list(c(State, Pars)), {
dx = x*(alpha - beta*y)
dy = -y*(gamma - delta*x)  #reverse to just echo the line above :)
return(list(c(dx, dy)))
})
}

## where x reflects the prey population (hares, rabbits, lemurs, etc), 
## and y reflects the predator population (wolves, foxes, fossa, etc). 
## The prey population grows at a linear rate (alpha) and gets eaten by the predator 
## population at the rate of (beta) per predator. The predator gains a certain amount 
## vitality by eating the prey at a rate (delta), while dying off at another rate (gamma). 
## Let's set some initial parameter values, and a starting population structure, and 
## time-vector, as previously: 

Pars <- c(alpha = 2, beta = .5, gamma = .2, delta = .6)
State <- c(x = 10, y = 10)
Time <- seq(0, 150, by = 1)

## We can then run the model and plot the results: 
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time)) #run model

## plot results
matplot(out[,1],out[,-1], type = "l",    #separate 1st col (time), plot against other 2
        xlab = "time", ylab = "population", lty=1)  
legend("topright", c("Lemurs", "Fossa"), lty = c(1,1), col = c(1,2), box.lwd = 0)

## We can start with different starting numbers of rabbits and fosa, and overlay the two plots (using the command 'add=TRUE' in matplot). 
## Ro distinguish the two runs, we use dashed lines for the second set of starting values using lty=2 
## (colours still represent predator and prey, as in the legend):  

State2 <- c(x = 5, y = 10)
out2 <- as.data.frame(ode(func = LotVmod, y = State2, parms = Pars, times = Time))
matplot(out[,1],out[,-1], type = "l", 
        xlab = "time", ylab = "population", lty=1)  #plot out our previous results
matplot(out2[,1],out2[,-1], type = "l", lty=2,add=TRUE)  #add the new runs
legend("topright", c("Lemurs", "Fossa"), lty = c(1,1), col = c(1,2), box.lwd = 0)

## ?? How does this change the dynamics?
## makes them lag in time

## After the initial transient (detectable by slightly higher or lower peaks for the 
## first few cycles), the two simulations start repeating the same pattern over and over. 

## We can try the same with different parameters, to understand what the parameters do, 
## again overlaying the new results for comparison. Here, we double the parameter alpha, 
## which determines prey reproduction: 
  
Pars2 <- c(alpha = 4, beta = .5, gamma = .2, delta = .6)
out3 <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars2, times = Time))
matplot(out[,1],out[,-1], type = "l", 
        xlab = "time", ylab = "population", lty=1)  
matplot(out3[,1],out3[,-1], type = "l", lty=2,add=TRUE)  
legend("topright", c("Lemurs", "Fossa"), lty = c(1,1), col = c(1,2), box.lwd = 0)



##############################################################################
## Susceptible-Infected-Recovered models - LIVE
##############################################################################

## Let's think about another two species model: dynamics for a directly 
## transmitted immunizing infection, with a short generation time (e.g., ~2 weeks), 
## like measles. We need to keep track of 3 states, or types of individuals, i.e., 
## susceptible individuals ('S'), infected individuals ('I') and recovered individuals ('R').
## We can write the function that defines this process, assuming that the total population 
## size is constant, which means we don't need to keep track of the ('R') compartment. 
 
## 1. write a function for the SIR model





## We define a starting population with 500,000 susceptible individuals and 20 infected 
## individuals; and a list of parameters that define the processes, beta, the transmission 
## rate, and gamma, the recovery rate (defined as 1/infectious period): 

## 2. Start with 500,000 susceptible individuals and 20 infected individuals  

## 3. generate a list of parameters that define beta (transmission rate), gamma (the recovery rate), and population size



## 4. calculate R0



## ?? What values of R0 are important?

## As in the population example, we can use the function lsoda, here taking as our time-unit
# 1 day, setting a small time-step, and running it out across 1 year (365 days): 
## 5. use the function lsoda to run the sir function, generating output every day for one year, and save it into a dataframe called ts.sir




## 6. Look at the top rows of the resulting data frame from this run of the model pop.SI


## 7. Plot susceptible as green, infected as red, and add a legend



## The epidemic burns itself out (i.e., at the end, the number of infected individuals 
## is zero), but, interestingly, not all susceptible individuals are infected (i.e., the 
## number of susceptibles at the end is > 0). This is the phenomenon of herd immunity. 

## Now we can re-run the analysis for a less infectious pathogen, i.e. beta is lower, but
## with the same population parameters:
## 8. create a list of parameters with lower values called lower.values



## 9. generate one year of results using lower.values, called lower.ts.sir




## 10. plot the results of lower.ts.sir




## ?? How are these results different? 

## ?? Why don't we really need to plot R?



## TO DO: 
#### 1) Change the parameter values to see what happens if people recover slower (lower gamma). 
#### 2) Change the number of initial infected individuals to 1 -- will the outbreak start earlier? Later? 
### 3) This model reflects a 'closed' population. There are \textbf{no new births} entering the population. Can you imagine how you might change the model to include this? What are the main dynamical consequences? 
### 4) You can download time-series for a number of infections including measles from https://www.tycho.pitt.edu; you could explore seasonality for different pathogens, or patterns through space and time, and across different city sizes  



## References 

# Cawell, H. Matrix Population Models. 2001. Sinauer

# Salguero-Gomez, R., et al. The COMPADRE Plant Matrix Database: an open online repository for plant demography. Journal of Ecology, 2015. 103(1): p. 202-218.

# Some of this material is extended from DAIDD and MMED courses. 

