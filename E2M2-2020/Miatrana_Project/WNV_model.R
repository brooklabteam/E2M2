## Make a simple mechanistic model of WNV dynamics in wild birds in Madagascar,
## in which Hypsipetes madagascariensis serves as a reservoir source population
## for fatal infections in the other birds.

#clear your working directory
rm(list=ls())

#load a few libraries
library(deSolve) #for running the ode
library(ggplot2)
library(reshape2)

#write your differential equation function
WNV.SIR <- function(t,x,parms){
  
  SH = x[1] #Susceptible H. madagascariensis
  IH = x[2] #Infectious H. madagascariensis
  RH = x[3] #Recovered H. madagascariensis
  SB = x[4] #Susceptible other birds
  IB = x[5] #Infectious other birds
  
  
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(as.list(parms),{
    
    N = SH + IH + RH + SB + IB
    
    dSHdt <- parms$bH*(1-N/parms$K)*(SH+IH+RH)-parms$betaH*SH*IH-parms$muH*SH
    dIHdt <- parms$betaH*SH*IH - parms$muH*IH - parms$sigmaH*IH
    dRHdt <- parms$sigmaH*IH - parms$muH*RH
    dSBdt <- parms$bB*(1-N/parms$K)*(SB+IB)-parms$betaB*SB*IH-parms$muB*SB
    dIBdt <- parms$betaB*SB*IH - parms$muB*IB - parms$alphaB*IB
    
    
    list(c(dSHdt,dIHdt,dRHdt,dSBdt, dIBdt))
  })
}

#list your parameters:
params = list(bH= .01,#births of Hypsipetes per capita per day. could later make this a seasonal process
                betaH = .002,#per capita infectious contacts per day for an infected Hypsipetes to a susceptible Hypsipetes
                sigmaH = 1/14, #recovery rate, 1/duration of infection (in days). here 1/14, like measles, for Hypsipetes
                muH = 1/1825, #natural death rate for Hypsipetes (1/natural lifespan in days = 1/1825, or 5 years)
                bB= .01,#births of other birds per capita per day. could later make this a seasonal process
                betaB = .002,#per capita infectious contacts per day for an infected Hypsipetes to a susceptible other bird
                muB = 1/1825, #natural death rate for other birds (1/natural lifespan in days = 1/1825, or 5 years)
                alphaB = 1/10, #infection-induced mortality rate for other birds. 1/duration of lethal infection (in days)
                K = 500 #carrying capacity for all birds in the system
)


#list your initial conditions. start with 100 Hypsipetes, 1 of which is infected and 100 other birds, none of which is infected
xstart = c(SH = 99, IH = 1, RH = 0, SB=100, IB=0)

times = seq(1, 365*5, 1) #run model for 5 years, taking timesteps in days

out = as.data.frame(lsoda(y = xstart, times = times, func =WNV.SIR, parms = params))
head(out)
#rename columns

#reshape data for plotting in ggplot
out.long <- reshape(out, idvar = "time", timevar = "state", times=c("susceptible_Hypsipetes", "infectious_Hypsipetes", "recovered_Hypsipetes", "susceptible_other_bird", "infectious_other_bird"), varying = list(2:6), v.names="count", direction = "long")
head(out.long)
#rename columns
names(out.long) <- c("time", "state", "count")
#remove rownames
rownames(out.long) <- c()

#and plot
p1 <- ggplot(data=out.long) + geom_line(aes(x=time, y=count, color=state))
print(p1)


#For you, Miatrana:
# experiment with the parameters. Also try plotting it as proportions instead
# of counts. Also, try making the plot pretty and changing the colors

