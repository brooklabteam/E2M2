################ MECHANISTIC MODEL ###################
## We think that the REASON behind the declining forest might have to do with slash
## and burn agriculture. We decide to build a model describing conversion of forest
## to savanna via slash and burn. Instead of m and b, the slope and intercept of a line,
## we now want to estimate the the rate of slash and burn. We already know the rate
## of forest growth (r=1.01)
## Let's build and run both a continuous time version mechanistic model.

## Continuous Time Mechanistic Model
ForestSavannaCont <- function(t,y,parms){
  # The with() function gives access to the named values of parms within the
  # local environment created by the function
  with(c(as.list(y),parms),{
    N = For + Sav
    dFordt = r*((K-N)/K)*For - slash*For*Sav
    dSavdt = slash*((K-N)/K)*For*Sav 
    
    #Return forest to compare with data
    list(c(dFordt, dSavdt))
  })
}


## Set starting population
Mada.start <- c(For = 450000, #sq. km of forest in Mada
                Sav = 100000)  #sq. km of savanna in Mada  

## Set up time-steps and units - here go for slightly longer to allow model to equilibrate
times <- seq(1900,2010,by=1) 

## Define parameters 
parms <- c(r=1.01,  # the number of sq. km of forest produced by one existing sq. km of forest per year. r=1.01 means that, one year from now, we expect that 1 sq. km of forest will have grown to 1.01 sq. km of forest
            slash = .00000045, # number of sq km of savanna produced per sq km of slashed forest each year. slash = 1 would mean successful 1:1 conversion. slash = .8 would mean some forest gets slashed and still goes back to forest
            K = 900000) # carrying capacity for forest. sq. km of land available for forest/savanna

## And run model
contModTree <- data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, parms = parms))
names(contModTree) = c("yr", "forest")

## And plot with data.
with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 
lines(ContMod$yr, ContMod$forest, col = "blue", lwd = 2)


## It looks like it converts forest too quickly!
## We'll need to adjust the rates. We can do this by fitting the models to the data.

################ FITTING MECHANISTIC MODELS ###################

## To do this, we will make three functions:
## (1) A function that runs our model (we already made this above).
## (2) A function that compares our model with the data via a statistical test (we'll
## use least squares)
## (3) A wrapper function that minimizes the difference between our model and the
## data (here, we'll use R's function 'optim' to minimize the function we write 
## under #2)

## Write our comparison function:
sum.sq.mech = function(slash, Mada.start, times, data, r, K){
  ## run your continuous model
  parms = c(r = r,
            slash = slash,
            K = K)
  
  model.out = data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, parms = parms))
  names(model.out) = c("yr", "forest")
  
  ## remember, that we ran our model longer than the time period for which we have data,
  ## so we can't compare the earlier timesteps. we can subset the model output to compare
  ## only those chunks for which we have data.
  
  ## find the year your data starts.
  min.dat = min(data$yr)
  
  model.out = model.out[model.out$yr >= min.dat,]
  
  ## and compare the output of the model with the data
  sum.sq = sum((model.out$forest - data$forest)^2)
  
  return(sum.sq)
  
}

##test it on one run of your model
sum.sq.mech(slash=.00000045,
            Mada.start=Mada.start,
            times = times,
            r=1.01,
            K= 900000,
            data=treedata)

## Write our wrapper function
wrap.fit.mechanistic = function(Mada.start, times, r, K, data){
  
  ## make a sequence of guesses for the rate
  slash.list = seq(0, .000001, .0000001)
  
  ## then, search for the minimum
  sum.sq.lst = list()
  for (i in 1:length(slash.list)){
    sum.sq.lst[[i]] = sum.sq.mech(slash=slash.list[i], Mada.start=Mada.start, r=r, K=K, times=times, data=data)
  }
  sum.sq.lst = c(unlist(sum.sq.lst))
  
  plot(slash.list, sum.sq.lst, type="b")
  
  ##find the minimum
  fit.slash = slash.list[sum.sq.lst==min(sum.sq.lst)]
  
  return(fit.slash)
}

## Fit the continuous model. 

fit.continuous = wrap.fit.mechanistic(Mada.start=Mada.start,
                                      times = times,
                                      r=1.01,
                                      K= 900000,
                                      data=treedata)


## Run your continuous model with your fitted parameters and plot with the data
new.parms= c(r=1.01, 
                   slash = fit.continuous,
                   K = 900000)

run.fit.continuous.new  <- data.frame(lsoda(y = Mada.start, times = times, func = ForestSavannaCont, parms = new.parms))
names(run.fit.continuous.new) = c("yr", "forest")

with(treedata,
     plot(yr, forest, type="b", col="green", ylim=c(0, 800000), ylab = "sq. km forest", xlab = "", lwd=2)) 
lines(contModTree$yr, contModTree$forest, type="l", col="yellow", lwd=2)
lines(run.fit.continuous.new$yr, run.fit.continuous.new$forest, type="l", col="orange", lwd=2)
legend("bottomleft", legend=c("data", "model (guess)", "model (fitted)"), col=c("green", "yellow", "orange", lwd=2))

