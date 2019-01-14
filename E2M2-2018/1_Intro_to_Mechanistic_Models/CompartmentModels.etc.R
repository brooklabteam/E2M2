## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----Popdata, include=TRUE-----------------------------------------------
setwd("/Users/jessicametcalf/Downloads/E2M2-Dropbox-2018/E2M2.Metcalf.Lectures/Rcode/") 
pop.data <- read.csv("WorldBankPop.csv")
head(pop.data[,1:10])

## ----GetMada-------------------------------------------------------------
chs <- which(pop.data[,1]=="Madagascar",arr.ind=TRUE)
mada.pop <- as.numeric(pop.data[chs,3:ncol(pop.data)])

## ----GetFrance-----------------------------------------------------------
chs <- which(pop.data[,1]=="France",arr.ind=TRUE)
france.pop <- as.numeric(pop.data[chs,3:ncol(pop.data)])

## ----GetGrowth-----------------------------------------------------------
growth.mada.pop <- mada.pop[2:length(mada.pop)]/mada.pop[1:(length(mada.pop)-1)]
growth.france.pop <- france.pop[2:length(france.pop)]/france.pop[1:(length(france.pop)-1)]

## ----PlotCountries-------------------------------------------------------
par(mfrow=c(1,2)) #create a matrix of plots with 1 row and 2 cols
plot(1960:2016,mada.pop, type="l",xlab="",  #plot mada pop data
     ylab="Population size (World Bank estimate)", 
     ylim=range(c(mada.pop,france.pop),na.rm=TRUE))  #define limits of the y axis
points(1960:2016,france.pop, type="l",lty=2) #plot France pop data
legend("topleft",legend=c("Madagascar","France"),lty=1:2,bty="n")

plot(1960:2015,growth.mada.pop, type="l",xlab="", ylab="N(t+1)/N(t)",
     ylim=range(c(growth.mada.pop,growth.france.pop),na.rm=TRUE))
points(1960:2015,growth.france.pop, type="l",lty=2)
abline(h=1)

## ----LoopNdiscrete1------------------------------------------------------
N <- rep(NA,100) #create a vector for the population
N[1] <- 10       #set N(0) to be 10, as described
lambda <- 1.1    #set lambda, as described

## ----LoopNdiscrete-------------------------------------------------------
#loop t from 2 to 100, noting that t=1 corresponds to N(0)
for (t in 2:100) N[t] <- N[t-1]*lambda  

## ----directNdiscrete1----------------------------------------------------
Ndirect <- N[1]*lambda^(0:99) 

## ----compare-------------------------------------------------------------
plot(0:99,N, pch=15,col="red", xlab="time", ylab="N")
points(0:99,Ndirect,pch=19,col="blue", cex=0.5)

## ----dPdt----------------------------------------------------------------
pgrow <- function(t,y,parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    dPdt <- r*P
    list(c(dPdt))
  })
}

## ----run.dPdt------------------------------------------------------------
library(deSolve)              # Load libary to be used for numerical integration
Pop <- c(P=10)                # Define starting variables (P(0)=10, as above)
values <- c(r=log(lambda))    # Define parameters - use an approx of the lambda used in the discrete time ex above
time.out <- seq(0,100,by=0.1) # Set up time-steps for variables - here go for 100 years
pgrow(t=0,y=Pop,parms=values)  #try out the function for t=0

## ----test.dPdt-----------------------------------------------------------
delta.t <- 0.1                    # Set a small value for delta.t (0.1 day)
pop.next <- Pop + unlist(pgrow(t=time,y=Pop,parms=values)) * delta.t
pop.next

## ----run.test.dPdt-------------------------------------------------------
Ncontinuous <- lsoda(
  y = Pop,               # Initial conditions for population
  times = time.out,      # Timepoints for evaluation
  func = pgrow,          # Function to evaluate
  parms = values         # Vector of parameters
  )

## ----check.test.dPdt-----------------------------------------------------
head(Ncontinuous)

## ----check2.test.dPdt----------------------------------------------------
subset(Ncontinuous,time==2)

## ----check3.test.dPdt----------------------------------------------------
plot(Ncontinuous[,"time"],               # Time on the x axis
     Ncontinuous[,"P"],                  # Number infected (I) on the y axis
     xlab = "Time in days",     # Label the x axis
     ylab = "Pop size",         # Label the y axis
     main = " ",                # Plot title
     xlim = c(0,100),           # 
     type = "l",                # Use a line plot
     bty = "n")                 # Remove the box around the plot
points(1:100,N, pch=19,col="red",cex=0.2)  #compare with our discrete time

## ----LoopNdiscreteSep1---------------------------------------------------
N <- rep(NA,100)                            #create a vector for the population
N[1] <- 10                                  #set N(0) to be 10, as described
b <- 0.8; s=0.21                            #set b and r as described
for (t in 2:100) N[t] <- N[t-1]*b+ N[t-1]*s #loop, noting that t=1 corresponds to N(0) 

## ----LoopNdiscreteSepStoch-----------------------------------------------
Nstoch <- matrix(NA,50,100) #create a matrix to store the results
#                           # of the simulation, every row is a different run
Nstoch[,1] <- 10            #set N(0) to be 10, for every run. 
for (j in 1:50){            #loop over sims
for (t in 2:100) {          #loop over time
    Nstoch[j,t] <- rpois(1,Nstoch[j,t-1]*b)+ rbinom(1,Nstoch[j,t-1],s)
}}

## ----LoopNStochPlot------------------------------------------------------
matplot(0:99,t(Nstoch), xlab="", ylab="",
        type="l",lty=1,col="grey")  #'matplot' plots matrices; and extension of plot
points(0:99,N, type="l",lwd=2)

## ----ParamsAgeStruct-----------------------------------------------------
sj=0.2  #juvenile survival
aj=0.8  #juvenile aging 
sa=0.8  #adult survival
f=1     #adult fertility

## ----Matrix--------------------------------------------------------------
A <- matrix(0,2,2)
A[1,1] <- sj*(1-aj)  #juvenile to juvenile = not aging and surving as a juvenile
A[2,1] <- sj*aj      #juvenile to adult =  aging and surving as a juvenile
A[2,2] <- sa         #adult to adult =  surviving as an adult
A[1,2] <- f          #adult to juvenile = survivng and reproducing  

## ----MatrixEst-----------------------------------------------------------
n.start <- c(5,10)    #start with pop of 5 juveniles and 10 adults
n.next <- A%*%n.start #get pop one time-step latter using matrix multiplication (%*%)
n.next                #more juveniles, but fewer adults...

## ----MatrixLambda--------------------------------------------------------
eigen(A)$value[1]      #extract the pop growth rate using the eigenvalue; <1 indicates pop shrinking

A[1,2] <- f*2          #double the fertility 
eigen(A)$value[1]      #extract the pop growth rate using the eigenvalue; now pop growing

## ----Compadre------------------------------------------------------------
# Code not run - if you have downloaded the Compadre data, then you can load it 
#load("/Users/cmetcalf/Downloads/COMPADRE_v.4.0.1.RData")

# and then...

# names(compadre)           #check what it contains
# head(compadre$metadata)   #look at the metadata to see how its structured 
                            #indexes here correspond to indexes in the list of mat
# compadre$mat[[1]]         #pull out the 1st set of matrices - 
                            #in this list A is the full matrix;
#                           # U is only survival and growth, 
                            # F is fertility, and C is clonal reproduction (if there is any)
# eigen(compadre$mat[[1]]$matA)$value[1]  #get the population growth rate for this species 
#                                         #publication (defined by first row metadata)
# tmp <- which(compadre$metadata$Country=="MDG") #Find all the matrices relating to Madagascar
# compadre$metadata[tmp,]                        #Look at their metadata
# compadre$mat[[tmp[1]]]                         #If you want to manipulate in more depth - index...


## ----MatrixLV------------------------------------------------------------
LotVmod <- function (Time, State, Pars) {
    with(as.list(c(State, Pars)), {
        dx = x*(alpha - beta*y)
        dy = -y*(gamma - delta*x)  #reverse to just echo the line above :)
        return(list(c(dx, dy)))
    })
}

## ----MatrixLVpars--------------------------------------------------------
Pars <- c(alpha = 2, beta = .5, gamma = .2, delta = .6)
State <- c(x = 10, y = 10)
Time <- seq(0, 150, by = 1)

## ----MatrixLVplot--------------------------------------------------------
out <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars, times = Time)) #run
 
matplot(out[,1],out[,-1], type = "l",    #separate 1st col (time), plot against other 2
        xlab = "time", ylab = "population", lty=1)  
legend("topright", c("Rabbits", "Foxes"), lty = c(1,1), col = c(1,2), box.lwd = 0)

## ----MatrixLVplot2-------------------------------------------------------
State2 <- c(x = 5, y = 10)
out2 <- as.data.frame(ode(func = LotVmod, y = State2, parms = Pars, times = Time))
matplot(out[,1],out[,-1], type = "l", 
        xlab = "time", ylab = "population", lty=1)  #plot out our previous results
matplot(out2[,1],out2[,-1], type = "l", lty=2,add=TRUE)  #add the new runs
legend("topright", c("Rabbits", "Foxes"), lty = c(1,1), col = c(1,2), box.lwd = 0)

## ----MatrixLVplot3-------------------------------------------------------
Pars2 <- c(alpha = 4, beta = .5, gamma = .2, delta = .6)
out3 <- as.data.frame(ode(func = LotVmod, y = State, parms = Pars2, times = Time))
matplot(out[,1],out[,-1], type = "l", 
         xlab = "time", ylab = "population", lty=1)  
matplot(out3[,1],out3[,-1], type = "l", lty=2,add=TRUE)  
legend("topright", c("Rabbits", "Foxes"), lty = c(1,1), col = c(1,2), box.lwd = 0)

## ----SIRode1-------------------------------------------------------------
sir <- function(t,y,parms){
  with(c(as.list(y),parms),{
    dSdt <- -beta*S*I/N
    dIdt <- beta*S*I/N - gamma*I
    list(c(dSdt,dIdt))
  })
}

## ----SIRode2-------------------------------------------------------------
pop.SI <- c(S = 500000,I = 20)  
values <- c(beta = 3.6,          # Transmission coefficient
            gamma = 1/5,         #recovery rate
            N=800000)            # population size (constant)

## ----SIRode3-------------------------------------------------------------
R0 <- values["beta"]/values["gamma"] # value of R0
R0

## ----SIRode4-------------------------------------------------------------
ts.sir <- data.frame(lsoda(
  y = pop.SI,                   # Initial conditions for population
  times = seq(0,365,by=0.1),   # Timepoints for evaluation
  func = sir,                   # Function to evaluate
  parms = values                # Vector of parameters
  ))

## ----SIRodePlot----------------------------------------------------------
plot(ts.sir[,"time"],ts.sir[,"I"], xlab="Time", ylab="Infected", type="l")

## ----FunctionSeasonSIR---------------------------------------------------
simTSIR <- function(I0=10,S0=1000,N=300000,      #starting no infected, susceptible, total pop size
                    beta=15,alpha=0.97,gamma=0.2,#transmission, alpha, and recover
                    vacc.cover=rep(0,26),        #one level of vaccination for every week of the year
                    births=rep((30*300/26),26),  #one level of birth for every week of the year
                    Tmax=2600){                  #maximum desired time-span. 
  #set up storage for the population; vectors of length Tmax
  Istore <- Sstore <- rep(NA,Tmax)
  Istore[1] <- I0
  Sstore[1] <- S0
  
  #set up the seasonality (allowing transmission to fluctuate over the course of the year)
  beta.seas <- beta*(1+gamma*cos(2*pi*(1:26)/26))    
  #create an index, so for each of the Tmax time-steps, we know what season it is
  seas.index <- rep(1:26,Tmax/26)  
  
  #loop over the time-series calculating new infecteds and susceptibles.   
  for (t in 2:Tmax){
    Istore[t] <-  beta.seas[seas.index[t]]*(Istore[t-1]^alpha)*Sstore[t-1]/N
    Sstore[t] <- Sstore[t-1]+(births[seas.index[t]]*(1-vacc.cover[seas.index[t]]))-Istore[t]
  }
  return(list(Istore=Istore,Sstore=Sstore))
}

## ----plotTSIRsim---------------------------------------------------------
tmp <- simTSIR() ## by just leaving the arguments blank, R will use all the defaults.  
plot(tmp$Istore[2200:2400], type="l", xlab="Time", ylab="Incidence")  
        #(leave out the first 2200 time-steps to get rid of transients)
tmp1 <- simTSIR(births=rep((15*300/26),26)) ##half the birth rate, and rerun
points(tmp1$Istore[2200:2400], type="l",col=2)  #add this to the previous plot

## ----measlesData---------------------------------------------------------
meas = read.table("meas.csv", sep = ",", header = TRUE) 
names(meas)

## ----measlesDataPlot-----------------------------------------------------
plot(meas$time, meas$London, type = "b", xlab="Time (biweeks)", ylab="N cases")

## ----suscRecon-----------------------------------------------------------
cum.reg <- smooth.spline(cumsum(meas$B), cumsum(meas$London), df=5) 
D <- - resid(cum.reg) #The residuals
plot(cumsum(meas$B), cumsum(meas$London), type="l",
        xlab="cumulative births", ylab="cumuative incidence")
lines(cum.reg)
abline(a=0,b=1)

## ----under.report--------------------------------------------------------
rr <- predict(cum.reg, deriv=1)$y
summary(rr)

## ----correctIncidence----------------------------------------------------
Ic <- meas$London/rr
Dc <- D/rr

## ----makeVarialbles------------------------------------------------------
seas = rep(1:26, 21)[1:545]
lInew = log(Ic[2:546])
lIold = log(Ic[1:545])
Dold = Dc[1:545]

## ----makeS0--------------------------------------------------------------
N <- 3.3E6
Smean <-  seq(0.02, 0.2, by=0.001)*N
offsetN <- rep(-log(N), 545)

## ----makeS0llik----------------------------------------------------------
llik <- rep(NA, length(Smean))

## ----profile Likelihood--------------------------------------------------
for(i in 1:length(Smean)){
  lSold = log(Smean[i] + Dold)
  glmfit = glm(lInew ~ -1 +as.factor(seas) + lIold + offset(lSold+offsetN))
    llik[i] = glmfit$deviance / 2
  }
par(mfrow=c(1,1))
plot(Smean/3.3E6, llik, ylim=c(min(llik),35), xlab="Sbar", ylab="neg log-lik", type="l")

## ----profileLikelihoodBest-----------------------------------------------
lSold = log(Smean[which(llik==min(llik))] + Dold)
glmfit = glm(lInew ~ -1 +as.factor(seas) + lIold + offset(lSold+offsetN))
summary(glmfit)

## ----profile Likelihood for plot-----------------------------------------
beta=exp(glmfit$coef[1:26])
ubeta=exp(glmfit$coef[1:26]+summary(glmfit)$coef[1:26, 2])
lbeta=exp(glmfit$coef[1:26]-summary(glmfit)$coef[1:26, 2])

plot(beta, type="l", xlab="Biweek in year", ylab=expression(beta))
points(ubeta, type="l",lty=3)
points(lbeta, type="l",lty=3)



## ----simTSIR2------------------------------------------------------------
SimTsir2=function(beta, alpha, B, N,  
                  inits = list(Snull = 0, Inull = 0), type = "det"){
    type <- charmatch(type, c("det", "stoc"), nomatch = NA)
    if(is.na(type))
        stop("method should be \"det\", \"stoc\"")
    IT <- length(B)
    s <- length(beta)
    lambda <- rep(NA, IT)
    I <- rep(NA, IT)
    S <- rep(NA, IT)
    I[1] <- inits$Inull
    lambda[1] <- inits$Inull
    S[1] <- inits$Snull
    for(i in 2:IT) {
        lambda[i] = beta[((i - 2) %% s) + 1] * S[i - 1] * (I[i - 1]^alpha)/N
        if(type == 2) {
                I[i] = rpois(1, lambda[i])
            }
        if(type == 1) {
            I[i] = lambda[i]
}
        S[i] =S[i - 1] + B[i] - I[i]
    }
    return(list(I = I, S = S))

}

## ----simTSIR2show--------------------------------------------------------
sim <- SimTsir2(beta=exp(glmfit$coef[1:26]), alpha=glmfit$coef[27],
                B=meas$B, N=N, 
                inits=list(Snull=Dc[1]+Smean[which(llik==min(llik))], Inull=Ic[1]))
plot(sim$I, type="b", xlab="Time (biweeks)", ylab="Number of cases simulated")
lines(exp(lInew), col="red")
legend("topleft",
       legend=c("sim", "Ic"),
       lty=c(1,1),
       pch=c(1,0),
        col=c("black", "red"))

