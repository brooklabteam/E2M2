# fit the model to these new 

rm(list=ls())
#set wd
#.libPaths("~/R/x86_64-pc-linux-gnu-library/3.4")

library(epitools)
library(dplyr)
library(plyr)
library(cowplot)
library(deSolve)
library(mgcv)
library(lubridate)
library(matrixcalc)
library(Matrix)
library(dplyr)
library(plyr)
library(ggplot2)
library(mvtnorm)

edits
#load data
#setwd("/Users/caraebrook/Documents/R/R_repositories/bat-age-seroprev/JAE_Sept2018/revisions_Dec2018/refits/")
#load("seasonal_prev.Rdata")



#and remove the irrelevant columns
#dat.tot <- dplyr::select(dat.new, -(month), -(year), -(titer), -(titer_resid), -(doy))
#head(dat.tot)

homewd = "/Users/carabrook/Developer/bat-sero-dynamics"

dat <- read.csv(file = paste0(homewd, "/data/EidGhVfit.csv"), header = T, stringsAsFactors = F)
dat <- read.csv(file = paste0(homewd, "/data/EidNiVfit.csv"), header = T, stringsAsFactors = F)
names(dat)[names(dat)=="bat_species"] <- "species"
names(dat)[names(dat)=="antigen"] <- "type"
names(dat)[names(dat)=="seropos"] <- "prev"
names(dat)[names(dat)=="sampleid"] <- "sampleID"



#helper functions
buildTMat <- function(c, Npop, age.classes, surv.biwk, surv.juv.biwk, mu.sick,	beta, recov, sigma, sero, wane, add.inf.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  #if (length(beta)==1) beta <- rep(beta,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
                             (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
                        waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #psi = 1-exp(-sero*sum(I_m)) #force of seroconversion, also depends on density of I class
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi)))
  psi = 1-exp(-sero*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,4,4) # MSIR
  
  Tmat <- matrix(0,4*nage,4*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi-psi
    mat1[3,2] <- foi
    mat1[3,3] <- 1-recov
    mat1[4,2] <- psi
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[2,4] <- sigma[j]
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],4);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
      # surv[4] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*4+1):(j*4+4),((j-1)*4+1):(j*4)] <- mat1*surv*age.rate[j]
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat2*surv*(1-age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, mu.sick,	beta, recov, sigma, sero, wane, add.inf.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #psi = 1-exp(-sero*sum(I_m)) #force of seroconversion, also depends on density of I class
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #vector
  psi = 1-exp(-sero*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,4,4) # MSIR
  
  Tmat <- matrix(0,4*nage,4*nage) #MSIR for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]-psi
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,2] <- psi
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[2,4] <- sigma[j]
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],4);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){ #we'll let seropositives die more too
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
      # surv[4] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*4+1):(j*4+4),((j-1)*4+1):(j*4)] <- mat1*surv*age.rate[j]
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat2*surv*(1-age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*4+1):(j*4),((j-1)*4+1):(j*4)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRN <- function(c, Npop, age.classes, surv.biwk, surv.juv.biwk, mu.sick, boost, beta, recov, sigma, wane, add.inf.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  #if (length(beta)==1) beta <- rep(beta,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi)))
  Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,5,5) # MSIRN
  
  Tmat <- matrix(0,5*nage,5*nage) #MSIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi
    mat1[3,2] <- foi
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[4,5] <- Boost
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1-Boost
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],5);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
      #surv[4] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*5+1):(j*5+5),((j-1)*5+1):(j*5)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildTMat_MSIRN_age <- function(c, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, mu.sick, boost, beta, recov, sigma, wane, add.inf.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,5,5) # MSIRN
  
  Tmat <- matrix(0,5*nage,5*nage) #MSIRN for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[4,5] <- Boost
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1-Boost
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],5);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
      #surv[4] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*5+1):(j*5+5),((j-1)*5+1):(j*5)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}

buildTMat_MSIRNR2_age <- function(c, biwk, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, mu.sick, boost, beta, recov, sigma, sigma2, wane, add.inf.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #assign rate of titer elevation based on weeks of pregnancy/lactation - 
  #needs to start before the birth pulse, and presumably decline after
  #we'll take biweeks #18-26
  #if(biwk>=18){ #peak births
  # Boost = boost
  #}else{
  # Boost = 0
  #}
  # 
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,c,c) # MSIRNP
  
  Tmat <- matrix(0,c*nage,c*nage) #MSIRNP for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    #if the biweek is within the rep pulse, we say transmission from N to R is certain
    
    
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1-Boost
    mat1[5,6] <- sigma2 #fits an a different sigma2 by which you wane from R2 back into N post-boosting
    mat1[6,5] <- Boost
    mat1[6,6] <- 1-sigma2
    
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],c);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
      #surv[4] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}

#now try where N-class mums elevate titers (become R) seasonally...
buildTMat_MSIRNP_age <- function(c, biwk, Npop, age.classes, surv.biwk, age.brk, surv.juv.biwk, mu.sick, boost, beta, recov, sigma, sigma2, wane, add.inf.mort, age.rate, slope.wane){
  
  s <- nage <- length(age.classes)
  
  #assign rate of titer elevation based on weeks of pregnancy/lactation - 
  #needs to start before the birth pulse, and presumably decline after
  #we'll take biweeks #18-26
  if(biwk>=18){ #peak births
    Boost = boost
  }else{
    Boost = 0
  }
  # 
  #put the mortality rates end-to-end
  mort.vect <- c((1-surv.juv.biwk), rep((1-surv.biwk), (nage-1))) 
  
  #first, set up beta:
  if (length(beta)==1) beta <- rep(beta,nage) # 
  if (length(beta)==s) beta <- beta
  if (length(beta) > 1 & length(beta)<s){
    beta.list <- list()
    for (i in 1:length(beta)){
      beta.list[[i]] = rep(beta[i], age.brk[i])
    }
    beta = c(unlist(beta.list))
  } 
  #if (length(sero)==1) sero <- rep(sero,nage) # we assume the pop is perfectly mixed so this is irrelevant
  #if (length(boost)==1) boost <- rep(boost,nage) # we assume the pop is perfectly mixed so this is irrelevant
  if (length(sigma)==1) sigma <- rep(sigma,nage)
  if (length(mu.sick)==1) mu.sick <- rep(mu.sick, nage)
  if (length(age.rate)==1) age.rate <- rep(age.rate,nage)
  if (length(wane)==1) wane <- rep(wane,nage)
  waning.maternal = wane
  
  #then adjust the aging rate for age class 1 to reflect duration of mat. immunity
  #age.rate[1] = waning.maternal[1]
  
  
  #use logisitc regression - now with justin's amendment (i.e. correct for going in and out of classes)
  #try without
  #waning.maternal.here<-1-(exp(pmin(-wane+slope.wane*c(0,age.classes),20))/
  #                          (1+exp(pmin(-wane+slope.wane*c(0,age.classes),20))))
  #waning.maternal <- (waning.maternal.here[1:length(age.classes)]-
  #                     waning.maternal.here[2:(1+length(age.classes))])/waning.maternal.here[1:length(age.classes)]
  #waning.maternal[waning.maternal.here[1:length(age.classes)]<1e-8] <- 1
  
  #assume density dependent transmission - means you need the total infectious pop, so you need to isolate that first thing
  #first, transform to get all the infecteds together
  Npop_epi <- transform.vect(vec=Npop, s=s, c=c)
  I_m = mat_split(matrix(Npop_epi, ncol=1), r=s, c=1)[,,3] #always. for all model forms, I comes #3
  #foi = 1-exp(-beta*sum(I_m))
  #Boost = 1-exp(-boost*sum(I_m))
  foi = 1-exp(-beta*(sum(I_m)/sum(Npop_epi))) #will produce a vector. if there was structure in the age-contacts we would feed it a contact matrix too, but here we assume age classes are evenly mixed
  # Boost = 1-exp(-boost*(sum(I_m)/sum(Npop_epi))) #force of seroconversion, also depends on density of I class
  
  
  mat1 <- matrix(0,c,c) # MSIRNP
  
  Tmat <- matrix(0,c*nage,c*nage) #MSIRNP for every age cohort
  for (j in 1:nage) {
    
    #fill in epi matrix for each age class. this gives probability of transmission per biweek for a MSIR/MSIRS/MSRIR model
    #if the biweek is within the rep pulse, we say transmission from N to R is certain
    
    
    mat1[] <- 0
    mat1[1,1] <- 1- waning.maternal[j]
    mat1[2,1] <- waning.maternal[j]
    mat1[2,2] <- 1-foi[j]
    mat1[3,2] <- foi[j]
    mat1[3,3] <- 1-recov
    mat1[4,3] <- recov
    mat1[4,4] <- 1-sigma[j] #already included here. simply give sigma a value if you want to allow for waning immunity
    mat1[5,4] <- sigma[j]
    mat1[5,5] <- 1-Boost
    mat1[5,6] <- sigma2
    mat1[6,5] <- Boost
    mat1[6,6] <- 1-sigma2
    
    
    #put in surv. we first say it is equal across all infectious classes
    surv <- rep(1-mort.vect[j],c);
    
    #if we decide to add infection-induced mortality, we can do so here by replacing survival for the infecteds
    #in this case, looks like sick are only slighly more likely to add. option here will make them definitely die
    if (add.inf.mort==TRUE){
      surv[3] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
      #surv[4] <- pmax(1-(mu.sick[j]*mort.vect[j]),0)
    } 
    #print(c(j,surv[3]))
    #fill in Tmatrix
    if (j!=nage) { 
      #if you are not in the last age class, you go into next age class
      #multiply infection transitions times survival and aging rates
      Tmat[(j*c+1):(j*c+c),((j-1)*c+1):(j*c)] <- mat1*surv*age.rate[j]
      
      #and those that do not age are also subject to their own survival and transition rates:
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv*(1-age.rate[j])
      
      
      #except for the maternally immune - we don't let them transition between infection states
      #but maybe we do!
      #here, we ammend the epidemic transition matrix accordingly
      #mat2 <- mat1; mat2[1,1] <- 1; mat2[2,1] <- 0
      
      #and we fill in all the maternally immune correspondingly
      #when age.rate=1, this is 0, meaning that there is no survival across maternally immune categories--
      #you change class when you age up - this essentially says that for those that don't age up, if maternally immune, they stay in their current class. i think i disagree...
      #but we let that maternally immune class have a different aging rate (<1 year)
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat2*surv*(1-age.rate[j])
      #Tmat[((j-1)*5+1):(j*5),((j-1)*5+1):(j*5)] <- mat1*surv*(age.rate[j])
      
    } else {
      #stay in this age class if you are at the peak age
      Tmat[((j-1)*c+1):(j*c),((j-1)*c+1):(j*c)] <- mat1*surv
      
      
    }
    
  }
  
  return(Tmat)
  
}
buildFMatrix <- function(age.classes, adult_fec, surv.biwk, biwk){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.2
  }else if (biwk==3|biwk==25){
    new.fec = adult_fec*surv.biwk*.1
  }else if (biwk==4|biwk ==24){
    new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,4*nage,4*nage)
  
  for (j in 1:nage) {
    Fmat[1,((j-1)*4+1):(j*4)] <- c(0,0,fert.biwk[j],fert.biwk[j]) #the mat immune
    Fmat[2,((j-1)*4+1):(j*4)]<- c(fert.biwk[j],fert.biwk[j],0,0) #the susceptible (mom didn't get sick)
  }
  
  return(Fmat)
}
buildFMatrix_MSIRN <- function(age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3|biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,5*nage,5*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    if(N_stat=="matAB"){
      #try it assuming N moms produce maternally immune pups
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j]) #the mat immune offspring (so mom = I, R, and N)
      Fmat[2,((j-1)*5+1):(j*5)]<-  c(fert.biwk[j],fert.biwk[j],0,0, 0) #the susceptible offspring (mom = M (there are none in this class b/c too young) and S)
      
    }else if (N_stat=="matSus"){
      Fmat[1,((j-1)*5+1):(j*5)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0) #the mat immune offpring (so mom = I and R but not N)
      Fmat[2,((j-1)*5+1):(j*5)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j]) #the susceptible offspring (so mom= M, S, and N)
    }
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRNP <- function(c, age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3|biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,c*nage,c*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    #if(N_stat=="matAB"){
    #try it assuming N moms produce maternally immune pups
    # Fmat[1,((j-1)*c+1):(j*c)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j]) #the mat immune (so mom = I and R but not N)
    #Fmat[2,((j-1)*c+1):(j*c)]<- c(fert.biwk[j],fert.biwk[j],0,0, 0) #the susceptible (mom didn't get sick, so mom= N and S)
    
    #}else if (N_stat=="matSus"){
    #M, S, I, R, NM, NF, P
    
    Fmat[1,((j-1)*c+1):(j*c)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0,fert.biwk[j]) #the mat immune (so mom = I and R but not N)
    Fmat[2,((j-1)*c+1):(j*c)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j], 0) #the susceptible (mom didn't get sick, so mom= N and S)
    #}
    
    
    
  }
  
  return(Fmat)
}
buildFMatrix_MSIRNR2 <- function(c, age.classes, adult_fec, surv.biwk, biwk, N_stat){ 	#one for every age class
  #base your fertility on the biweek of the year
  #pop will grow a bit because higher chance of surviving to birth when births come earlier
  #but these total to the annual fec rate
  
  if(biwk==1){ #peak births
    new.fec = adult_fec*surv.biwk*.3
  }else if (biwk==2|biwk==26){
    new.fec = adult_fec*surv.biwk*.25
  }else if (biwk==3|biwk==25){
    new.fec = adult_fec*surv.biwk*.1
    #}else if (biwk==4|biwk ==24){
    #new.fec = adult_fec*surv.biwk*.05
  }else{
    new.fec = 0
  }
  s <- nage <- length(age.classes)
  
  fert.biwk <- c(0,rep(new.fec,(s-1)))
  
  
  
  #make matrix the same size as the transition
  Fmat <- matrix(0,c*nage,c*nage)
  
  for (j in 1:nage) { #no fertility in first age class (mat immune)
    #if(N_stat=="matAB"){
    #try it assuming N moms produce maternally immune pups
    # Fmat[1,((j-1)*c+1):(j*c)] <- c(0,0,fert.biwk[j],fert.biwk[j], fert.biwk[j]) #the mat immune (so mom = I and R but not N)
    #Fmat[2,((j-1)*c+1):(j*c)]<- c(fert.biwk[j],fert.biwk[j],0,0, 0) #the susceptible (mom didn't get sick, so mom= N and S)
    
    #}else if (N_stat=="matSus"){
    #M, S, I, R, NM, NF, P
    
    Fmat[1,((j-1)*c+1):(j*c)] <- c(0,0,fert.biwk[j],fert.biwk[j], 0,fert.biwk[j]) #the mat immune (so mom = I and R but not N)
    Fmat[2,((j-1)*c+1):(j*c)]<- c(fert.biwk[j],fert.biwk[j],0,0, fert.biwk[j], 0) #the susceptible (mom didn't get sick, so mom= N and S)
    #}
    
    
    
  }
  
  return(Fmat)
}
transform.vect <- function(vec, s, c){
  vec2 <- t(commutation.matrix(r=s,c=c))%*%vec
  return(vec2)
}
find.biweek = function(t, times){
  
  biwks <- sort(unique(round(revtrunc(times),4)))
  this.wk = round(revtrunc(times[t]), 4)
  this.biwk <- which(biwks==this.wk)
  
  
  return(this.biwk)
  
  
}
get.age.struct = function(pop, s){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=s, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  age.mat <- Reduce('+', age.mat.list)
  return(age.mat)
}
get.age.struct.M = function(pop, c){ #function that collapses population vector in disease state form to age structure
  age.mat <- mat_split(M=matrix(pop, ncol=1), r=c, c=1)
  age.mat.list <- c()
  for (i in 1:dim(age.mat)[3]){
    age.mat.list[[i]] <- age.mat[,,i]
  }
  #zage.mat <- Reduce('+', age.mat.list)
  age.mat <- lapply(age.mat.list, sum)
  age.mat.dat <- c( c(1:length(age.mat)), unlist(age.mat))
  names(age.mat.dat) <- c("age", "pop")
  return(age.mat)
}
prev.by.age <- function(dat){
  dat$n_age = sum(dat$count)
  dat$prev = dat$count/dat$n_age
  #and seroprev
  dat$seropos = sum(dat$count[dat$class=="M"], dat$count[dat$class=="R"])
  dat$seroprev = dat$seropos/dat$n_age
  return(dat)
}
revtrunc = function(x){
  newx = x - floor(x)
  return(newx)
}
spline.fit = function(data){
  #replace any Inf vals only the okay values
  #data <- data[complete.cases(data),]
  data$seroprevalence[data$seroprevalence==Inf] <- 1
  data$seroprevalence[data$seroprevalence<0] <- 0
  # if(length(data$age)>1){
  spline1 = with(data, smooth.spline(x=age, y=seroprevalence)) 
  #  return(spline1)
  # }else{
  #  return(NA)
  #  }
}
pred.fit <- function(spline1, ages.new){
  prev.comp = predict(spline1, x=ages.new)
  return(prev.comp$y)
}
get.seas.seroprev = function(dat){
  #get doy in bat calendar
  #convert to biweek
  #calc seroprev by biweek
  
  dat$doy <- yday(dat$date)
  #now correct for birthday of bat in questiob
  if(unique(dat$species=="Pteropus rufus")){
    bday <- yday("2015-10-01") #320
  } else if(unique(dat$species=="Eidolon dupreanum")){
    bday = yday("2015-11-01")
  }
  #write over for doy
  dat$new_doy = NA
  for (i in 1:length(dat$doy)){
    if(dat$doy[i] > bday){
      dat$new_doy[i] <- dat$doy[i] - bday
    } else if (dat$doy[i] <= bday){
      dat$new_doy[i] <- dat$doy[i] + (365-bday)
    }
  }
  
  dat$doy <- dat$new_doy
  dat <- dplyr::select(dat, -(new_doy))
  
  #now convert to biweek
  dat$biwk <- NA
  brk <- seq(0,365, by=14)
  for (i in 1:length(dat$biwk)){
    tmp <- brk[dat$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    dat$biwk[i] <- which(brk==tmp2)
  }
  
  #now calc seroprev by biweek
  biwk.sero <- ddply(dat, .(biwk), summarize, seropos = sum(prev), sero_lci = sum(prev_lci), sero_uci = sum(prev_uci), n=length(prev))  
  biwk.sero$seroprev = biwk.sero$seropos/biwk.sero$n
  biwk.sero$seroprev_lci = biwk.sero$sero_lci/biwk.sero$n
  biwk.sero$seroprev_uci = biwk.sero$sero_uci/biwk.sero$n
  biwk.sero$biwk = as.numeric(biwk.sero$biwk)
  
  return(biwk.sero)
}
get.mod.seas = function(mod.out){
  #convert time to doy
  mod.out$doy = mod.out$time*365
  #then to biwk
  brk <- seq(0,365, by=14) #doy breaks by biweek
  brk <- brk[-length(brk)]
  #brk.seq = 1:length(brk)
  mod.out$biwk <- NA
  for (i in 1:length(mod.out$biwk)){
    tmp <- brk[mod.out$doy[i] >=  brk]
    tmp2 <- tmp[length(tmp)] 
    mod.out$biwk[i] <- which(brk==tmp2)
  }
  #then return
  return(mod.out)
}
age.sero.plot = function(dat){
  with(dat, plot(age, seroprevalence, type = "b", ylim=c(0,1)))
}
get.seroprev.dat = function(data, vis_split, cutoff){
  
  visbin = seq(3,floor(max(data$age)), vis_split)
  #but breakdown the early years into more
  visbin = c(c(0,.5, 1, 2),   visbin)
  data$age_year <- NA
  
  for (i in 1:length(data$age_year)){
    tmp = visbin[data$age[i] > visbin]
    data$age_year[i] <- tmp[length(tmp)] 
  }
  
  if(cutoff=="mean"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev)/length(prev), count=length(prev))  
  }else if (cutoff=="uci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_uci)/length(prev_uci), count=length(prev_uci)) 
  }else if (cutoff=="lci"){
    dat.sum2 <- ddply(data, .(age_year), summarize,  prevalence=sum(prev_lci)/length(prev_lci), count=length(prev_lci)) 
  }
  
  #you want the midpoint between either end of the age class now
  vect.age.yr = sort(unique(data$age_year))
  dat.sum2$age_plot = NA
  
  for (i in 1:(length(dat.sum2$age_year)-1)){
    dat.sum2$age_plot [i] = ((dat.sum2$age_year[i + 1] - dat.sum2$age_year[i])/2) + dat.sum2$age_year[i]
  }
  dat.sum2$age_plot[length(dat.sum2$age_plot)] =  (ceiling(max(data$age)) - dat.sum2$age_year[length(dat.sum2$age_year)])/2  + dat.sum2$age_year[length(dat.sum2$age_year)]
  
  names(dat.sum2)  <- c("real_age", "prevalence", "count", "age_year") #<- names(dat.sum2.tmp)
  
  # dat.sum2 <-  rbind(dat.sum2.tmp, dat.sum2)
  
  dat.sum3 = dat.sum2
  dat.sum3$class = "seropositive"
  rownames(dat.sum3) <- c()
  
  return(dat.sum3)
  
}
build.pop.mat = function(surv, surv_juv, s, adult_fec){
  pop.mat = matrix(0,  nrow=(s-1), ncol = (s-1))
  diag(pop.mat) = surv
  diag(pop.mat)[1] = surv_juv
  col_s = c(rep(0, s-2), surv)
  pop.mat = cbind(pop.mat, col_s)
  row1 = c(0, rep((adult_fec*surv), (s-1))) #bats reproduce for the first time at the end of the second year of life. good for E. dup and P. ruf
  pop.mat = rbind(row1,pop.mat)
  return(pop.mat)
  
}#for stable age distribution
mat_split <- function(M, r, c){
  nr <- ceiling(nrow(M)/r)
  nc <- ceiling(ncol(M)/c)
  newM <- matrix(NA, nr*r, nc*c)
  newM[1:nrow(M), 1:ncol(M)] <- M
  
  div_k <- kronecker(matrix(seq_len(nr*nc), nr, byrow = TRUE), matrix(1, r, c))
  matlist <- split(newM, div_k)
  N <- length(matlist)
  mats <- unlist(matlist)
  dim(mats)<-c(r, c, N)
  return(mats)
}#for matrix slicing. give it the number of rows and columns you want in the resulting matrices
stack.age = function(dat, s){
  dat.new= list()
  for (i in 1:s){
    dat.new[[i]] = dat[i,]
  }
  dat.new = c(unlist(dat.new))
  return(dat.new)
}
stack.class = function(mat.dat, c){
  
  new.dat = list()
  for (i in 1:c){
    new.dat[[i]] = mat.dat[i,]
  }
  new.dat = c(unlist(new.dat))
  return(new.dat)
}

#and age-structured beta
sim.met.MSIR.age <- function(burnin, sim_pop, yrs, ntyr, age.brk, s, beta, sero, recov, mort, mort_juv, adult_fec, wane, slope.wane, sigma, mu.sick, add.inf.mort, model){
  #first pull out your data...
  c=4
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  # plot(stab.struct, xlab="Age", ylab="Proportion")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")
  
  
  #print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #print(i)
    Tmat <- buildTMat_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk=age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick,	beta=beta, sigma=sigma, recov=recov, sero=sero, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    #print(i)
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #print(biwk1)
    
    #feed into fertility matrix:
    Fmat <- buildFMatrix(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=20, c=4)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R")
  dat.tot$N = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  #par(mfrow = c(1,1))
  #with(dat.tot, plot(time, N, type="l")) #ylim =c(0,1.2*max(N))))
  #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  
  # #par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # 
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$N
  prop.tot$S = prop.tot$S/prop.tot$N
  prop.tot$I = prop.tot$I/prop.tot$N
  prop.tot$R = prop.tot$R/prop.tot$N
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #      plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R"), col = c("violet", "mediumseagreen","tomato","cornflowerblue"), lwd=1, cex=.3)
  # 
  
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #      plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R"), col = c("violet", "mediumseagreen","tomato","cornflowerblue"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R"))
  
  # colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  dat.tot$proportion = dat.tot$count/sim_pop
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >=(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #remove age at 0
  eq.sero.tot = subset(eq.sero.tot, age >0)
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  list.by.T[[1]] = c()
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = do.call("rbind", list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  #ggplot(data=out.mat) + geom_point(aes(age, seroprevalence, color=biwk), pch=1)
  #names(eq.sero.tot)[1] ="biwk"
  #plot one of them to get a sense
  #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRN.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, mort, mort_juv, adult_fec, wane, boost, slope.wane, sigma, mu.sick, add.inf.mort, model, N_stat){
  #first pull out your data...
  c=5
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRN_age(c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRN(age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N = prop.tot$N/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5]), cbind(dat.tot[,1], dat.tot[,6])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R", "N"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot[,1],prop.tot[,2]), cbind(prop.tot[,1], prop.tot[,3]), cbind(prop.tot[,1], prop.tot[,4]), cbind(prop.tot[,1], prop.tot[,5]), cbind(prop.tot[,1], prop.tot[,6])))
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S", "I", "R", "N"))
  
  
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R", "N"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >=(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #remove age at 0
  eq.sero.tot = subset(eq.sero.tot, age >0)
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  list.by.T[[1]] = c()
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = do.call("rbind", list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  
  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  
  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  
  #now, for simplicity, collapse to just age-seroprevalence cross-sectionally
  out.mat$age <- trunc(out.mat$age)
  #and because you can't have partial bats
  out.mat$seropos <- trunc(out.mat$seropos)
  out.mat$N <- trunc(out.mat$N)
  age.ser.dat <- ddply(out.mat, .(age), summarise, seropos = sum(seropos), N= sum(N))
  age.ser.dat$seroprevalence <- age.ser.dat$seropos/age.ser.dat$N
  
  ggplot(age.ser.dat) + geom_point(aes(x=age, y=seroprevalence, size=N)) +
    geom_line(aes(x=age, y=seroprevalence))                             
  
  #now convert to an actual data table, where each bat is given an age, and a serostatus
  age.split <- dlply(age.ser.dat, .(age))
  
  build.data <- function(df){
    
    dat = cbind.data.frame(age = rep(df$age, df$N), sero = 0)
    dat$index = rownames(dat)
    #randomly change a subset of seronegs to seropos
    seropos.ind <- sample(dat$index, size = df$seropos, replace = F)
    dat[seropos.ind,]$sero <- 1

    return(dat)
  }
  
  age.split <-lapply(age.split, build.data)
  
  ind.dat <- data.table::rbindlist(age.split)
  
  # now, sample at random from this distribution and see how many bats you 
  # need to distinguish your statistics
  
  ind.dat$index <- rownames(ind.dat)
  
  sample.ind <- sample(ind.dat$index, size = field_sample, replace = F)
  
  sub.sample <- ind.dat[as.numeric(sample.ind),]
  
  #and summarise and plot 
  sub.df <- ddply(sub.sample, .(age), summarise, seropos = sum(sero), N= length(sero))
  sub.df$seroprevalence <- sub.df$seropos/sub.df$N
  
  ggplot(sub.df) + geom_point(aes(x=age, y=seroprevalence, size=N)) +
    geom_line(aes(x=age, y=seroprevalence))                             
  # sample to distinguish patterns
  
  
  #return both seasonal and non-seasonal data
  #return(out.mat)#, seas.prev))
  return(age.ser.dat)
}
sim.met.MSIRNP.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, mort, mort_juv, adult_fec, wane, boost, slope.wane, sigma, sigma2, mu.sick, add.inf.mort, model, N_stat){
  #first pull out your data...
  c=6
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  P_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, P_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNP_age(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, sigma2=sigma2, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNP(c=c, age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N", "P")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N = prop.tot$N/prop.tot$pop
  prop.tot$P = prop.tot$P/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5]), cbind(dat.tot[,1], dat.tot[,6]), cbind(dat.tot[,1], dat.tot[,7])))
  
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R", "N", "P"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R", "N", "P"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot[,1],prop.tot[,2]), cbind(prop.tot[,1], prop.tot[,3]), cbind(prop.tot[,1], prop.tot[,4]), cbind(prop.tot[,1], prop.tot[,5]), cbind(prop.tot[,1], prop.tot[,6]), cbind(prop.tot[,1], prop.tot[,7])))
  
  
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S", "I", "R", "N", "P"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S", "I", "R", "N", "P"))
  
  colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy", 'P' = "deeppink")
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R", "N", "P"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >=(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  #seropositives are R, M, and P
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M" | eq.dat.sero$class=="P"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M" | eq.dat.sero$class=="P"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #remove age at 0
  eq.sero.tot = subset(eq.sero.tot, age >0)
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  list.by.T[[1]] = c()
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = do.call("rbind", list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  
  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  
  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  
  #names(eq.sero.tot)[1] ="biwk"
  #plot one of them to get a sense
  #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}
sim.met.MSIRNR2.age <- function(burnin, sim_pop, yrs, ntyr,  s, beta, age.brk, recov, mort, mort_juv, adult_fec, wane, boost, slope.wane, sigma, sigma2, mu.sick, add.inf.mort, model, N_stat){
  #first pull out your data...
  c=6
  # dat.tmp <- subset(dat.age1, species==species1 & type==type1)
  #then run model for 100 years 
  #can take more timesteps than the data
  times <-   seq(0, yrs, by =1/ntyr) # then, subtract last few biweeks so you end before the last birthpulse starts. 
  #means ctting last 4 biweeks
  #times <- times[1:(length(times)-4)]
  
  
  
  #first, we take our juvenile and adult survival rates and use them to get the stable age structure
  mat1 = build.pop.mat(surv=(1-mort), surv_juv=(1-mort_juv), s=(s), adult_fec = adult_fec)
  
  stab.struct = Re(eigen(mat1)$vector[,1])
  stab.struct <- stab.struct/sum(stab.struct)
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  #out of curiosity...
  lambda = round(max(Re(eigen(mat1)$value)), 8)
  print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
  #then, we use this stable age structure to make a bat population
  #gives counts of bats per age year
  bat.mat = stab.struct*sim_pop
  
  #introduce a few infecteds and run it out to equilibrium before you grab the data
  R_init = rep(0, s)
  I_init = rep(0, s); I_init[3] = 5 #comment out if you just want to check demography
  M_init = rep(0, s)
  N_init = rep(0, s)
  R2_init = rep(0, s)
  S_init = bat.mat - I_init 
  
  
  N_tot = cbind(M_init, S_init, I_init, R_init, N_init, R2_init)
  
  N_pop = vec(N_tot) #by disease status
  M_pop = vec(t(N_tot)) #by age class
  
  #make matrix to store your population as you go
  N_pop_ts <- matrix(NA, ncol = length(times), nrow(N_pop))
  N_pop_ts[,1] <- M_pop # by age class
  
  
  stab.struct <- get.age.struct.M(M_pop, c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")
  
  
  # print("start ts")
  #iterate SIR such that this gets transitioned each timestep. And change the births as you go
  for (i in 1:(length(times)-1)){
    #build matrices anew each time because foi depends on # infected
    #calculate biweek from timestep here:
    biwk1 <- find.biweek(t=i, times=times)
    #  print(i)
    Tmat <- buildTMat_MSIRNR2_age(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, sigma2=sigma2, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    
    
    #Tmat <- buildTMat_MSIRN_age_seas(biwk=biwk1, c=c, Npop= N_pop_ts[,i], age.classes=1:s, age.brk = age.brk, surv.biwk = (1-mort)^(1/ntyr), surv.juv.biwk = (1-mort_juv)^(1/ntyr), mu.sick=mu.sick, boost=boost,	beta=beta, sigma=sigma, recov=recov, wane= wane, add.inf.mort=add.inf.mort, age.rate=1/ntyr, slope.wane=slope.wane)
    # print(i)
    
    # print(biwk1)
    
    #feed into fertility matrix since births are dependent on biwk
    Fmat <- buildFMatrix_MSIRNR2(c=c, age.classes=1:s, adult_fec =adult_fec, surv.biwk = (1-mort)^(1/ntyr), biwk = biwk1, N_stat = N_stat)
    
    #make trans mat
    transMat <- Tmat + Fmat 
    
    #move forward in time
    nt1<-(transMat) %*% N_pop_ts[,i]
    N_pop_ts[,i+1] <- nt1
  }
  
  stab.struct <- get.age.struct.M(pop=N_pop_ts[,ncol(N_pop_ts)], c=c)
  stab.struct <- c(unlist(stab.struct))
  stab.struct <- stab.struct/(sum(stab.struct))
  #plot(stab.struct, xlab="Age", ylab="Proportion", type="b")  
  
  #transform whole vector to be split by class instead of age
  N_pop_ts <- transform.vect(vec=N_pop_ts, s=s, c=c)
  
  
  N_split_1 = mat_split(N_pop_ts, r=s, c=ncol(N_pop_ts))
  
  #transform array into list
  N_split = list()
  for (i in 1:dim(N_split_1)[3]){
    N_split[[i]] = N_split_1[,,i]
  }
  
  #now you have a list of state variables.
  #then take column sums of each and plot the class totals over time - note that there are 26 steps per year
  N_total = lapply(X=N_split, FUN=colSums)
  
  #plot both by age and total
  #times=seq(0,yrs,by =1/ntyr)
  dat.tot = cbind.data.frame(times,N_total)
  names(dat.tot) = c("time", "M", "S", "I", "R", "N", "R2")
  dat.tot$pop = rowSums((dat.tot[,2:ncol(dat.tot)]))
  #dat.tot$N = round(dat.tot$N, 2)
  # par(mfrow = c(1,1))
  # with(dat.tot, plot(time, pop, type="l")) #ylim =c(0,1.2*max(N))))
  # #with(dat.tot, plot(time[1:100], N[1:100], type="l", ylim =c(0,2000)))
  # 
  # par(mfrow = c(3,2))
  # 
  # with(dat.tot,
  #      plot(times, M, col="violet", type="l", xlab = "years", ylim =c(0,1.2*max(M))))
  # 
  # with(dat.tot,
  #      plot(times, S, col="mediumseagreen", type="l", xlab = "years", ylim =c(0,1.2*max(S))))
  # 
  # with(dat.tot,
  #      plot(times, I, col="tomato", type="l", xlab = "years", ylim =c(0,1.2*max(I))))
  # 
  # with(dat.tot,
  #      plot(times, R, col="cornflowerblue", type="l", xlab = "years", ylim =c(0,1.2*max(R))))
  # 
  # with(dat.tot,
  #      plot(times, N, col="navy", type="l", xlab = "years", ylim =c(0,1.2*max(N))))
  # 
  # with(dat.tot,
  #      plot(times, pop, col="black", type="l", xlab = "years", ylim =c(0,1.2*max(pop))))
  # 
  
  #eventually, could plot as proportions - these will stay stable at equilibrium
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N = prop.tot$N/prop.tot$pop
  prop.tot$P = prop.tot$R2/prop.tot$pop
  
  # par(mfrow = c(1,1))
  # with(prop.tot,
  #     plot(times, S, col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times, I, col="tomato", type="l"))
  # with(prop.tot, lines(times, R, col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times, M, col="violet", type="l"))
  # with(prop.tot, lines(times, N, col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  # #check annual dynamics - nice gentle peaks over a few weeks
  # with(prop.tot,
  #     plot(times[(length(times)-27*2): (length(times))], S[(length(times)-27*2): (length(times))], col="mediumseagreen", type="l", ylim = c(0, 1), ylab = "proportion", xlab = "years"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], I[(length(times)-27*2): (length(times))], col="tomato", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], R[(length(times)-27*2): (length(times))], col="cornflowerblue", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], M[(length(times)-27*2): (length(times))], col="violet", type="l"))
  # with(prop.tot, lines(times[(length(times)-27*2): (length(times))], N[(length(times)-27*2): (length(times))], col="navy", type="l"))
  # legend("topright", legend= c("M", "S", "I", "R", "N"), col = c("violet", "mediumseagreen","tomato","cornflowerblue", "navy"), lwd=1, cex=.3)
  # 
  # 
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5]), cbind(dat.tot[,1], dat.tot[,6]), cbind(dat.tot[,1], dat.tot[,7])))
  
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R", "N", "R2"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R", "N", "R2"))
  
  #colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy")
  #ggplot(data=dat.tot) + geom_line(aes(x=time, y=count, color=class)) + scale_color_manual(values=colz)
  
  prop.tot = data.frame(rbind(cbind(prop.tot[,1],prop.tot[,2]), cbind(prop.tot[,1], prop.tot[,3]), cbind(prop.tot[,1], prop.tot[,4]), cbind(prop.tot[,1], prop.tot[,5]), cbind(prop.tot[,1], prop.tot[,6]), cbind(prop.tot[,1], prop.tot[,7])))
  
  
  names(prop.tot) = c("time", "proportion")
  prop.tot$class = rep(c("M", "S", "I", "R", "N", "R2"), each= length(times))
  prop.tot$class = factor(prop.tot$class, levels=c("M", "S", "I", "R", "N", "R2"))
  
  colz = c('M'="violet", 'S' = "mediumseagreen", 'I' = 'tomato', 'R' = "cornflowerblue", 'N' = "navy", 'R2' = "deeppink")
  #ggplot(data=prop.tot) + geom_line(aes(x=time, y=proportion, color=class)) + scale_color_manual(values=colz) + coord_cartesian(ylim = c(0,1))
  
  # ### SEASONAL SEROPREV
  # #seas.tot <- dat.tot #take only last year at equilibrium
  # #seas.tot <- subset(seas.tot, time > (yrs-1))
  # 
  # #calc seroprevalence and plot with time
  # seas.tot$time = revtrunc(seas.tot$time)
  # seas.tot$sero = 0
  # seas.tot$sero[seas.tot$class=="M" | seas.tot$class=="R"] <- seas.tot$count[seas.tot$class=="M" | seas.tot$class=="R"]
  # 
  # seas.prev = ddply(seas.tot, .(time), summarise, seropos= sum(sero), n= sum(count))
  # seas.prev$seroprevalence = seas.prev$seropos/seas.prev$n
  # seas.prev = dplyr::select(seas.prev, -(n), -(seropos))
  # with(seas.prev, plot(time, seroprevalence, type="b", ylim=c(0,1)))
  # 
  ### AGE-SEROPREV HERE
  #within your stacked matrix, place them end to end, so all a1 are on top of all a2
  age.dat = lapply(N_split, stack.age, s=s)
  age.dat = do.call("rbind", age.dat)
  
  age.dat.tot = stack.class(age.dat, c=c)
  age.dat.tot = cbind.data.frame(rep(times, s*c), age.dat.tot)
  names(age.dat.tot) = c("times", "count")
  age.dat.tot$class = rep(c("M", "S", "I", "R", "N", "R2"), each= length(times)*s)
  age.dat.tot$age = rep(rep(seq(0,(s-1),1), each=length(times)), c)
  
  #and plot age classes over time
  age.dat.tot$age = factor(age.dat.tot$age)
  #ggplot(data=age.dat.tot) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #just seropositives 
  # ggplot(data=subset(age.dat.tot, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  #now we take the last year - at equilibrium - and we look at seasonal seroprevalence in that year
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >=(yrs-1))
  
  #plot to see what things look like across a year
  #ggplot(data=eq.dat) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  #ggplot(data=subset(eq.dat, class=="M" | class=="R")) + geom_line(aes(x=times, y=count, color=age, linetype=class))
  
  
  ####STABLE AGE STRUCTURE
  #check that the age-freq distribution is stable and same as beginning
  #age.struc = eq.dat[eq.dat$times==5,]
  age.struc = eq.dat[eq.dat$times==yrs,]
  age.struc = ddply(age.struc, .(age), summarize, tot = sum(count))
  age.struc$n = sum(age.struc$tot)
  age.struc$proportion = age.struc$tot/age.struc$n 
  age.struc$age= as.numeric(as.character(age.struc$age))
  #with(age.struc, plot(age,proportion, type="b"))
  
  
  #then take age-seroprevalence at each time point in this year and plot as a panel
  eq.dat.sero = eq.dat
  eq.dat.sero$age = as.numeric(as.character(eq.dat.sero$age))
  eq.dat.sero$age = eq.dat.sero$age + revtrunc(eq.dat.sero$times)
  eq.dat.sero$sero = 0
  #seropositives are R, M, and P
  eq.dat.sero$sero[eq.dat.sero$class=="R" | eq.dat.sero$class=="M" | eq.dat.sero$class=="R2"] <- eq.dat.sero$count[eq.dat.sero$class=="R" | eq.dat.sero$class=="M" | eq.dat.sero$class=="R2"]
  eq.sero.tot <- ddply(eq.dat.sero, .(times, age), summarize, seropos = sum(sero), N= sum(count))
  eq.sero.tot$seroprevalence = eq.sero.tot$seropos/eq.sero.tot$N
  
  eq.sero.tot  %>% mutate_if(is.factor, as.character) ->  eq.sero.tot
  eq.sero.tot$age = as.numeric(eq.sero.tot$age)
  eq.sero.tot$seroprevalence = as.numeric(eq.sero.tot$seroprevalence)
  #remove age at 0
  eq.sero.tot = subset(eq.sero.tot, age >0)
  
  #make a list of age-seroprev by biweek across the year
  #drop the first which is really the last
  list.by.T = dlply(eq.sero.tot, .(times))
  list.by.T[[1]] = c()
  #convert them into a matrix by biweek to keep track for fitting
  out.mat = do.call("rbind", list.by.T)
  rownames(out.mat) = c()
  uni.times = unique(out.mat$times)
  biwk.seq = 1:26
  out.mat$biwk = NA
  for (i in 1:length(biwk.seq)){
    out.mat$biwk[out.mat$times==uni.times[i]] <- biwk.seq[i]  
  }
  out.mat$biwk <- as.factor(out.mat$biwk)
  
  #assign months to each biwk and plot like that (based on the species)
  if (mort==.207){
    month.list <- c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }else if (mort==.456){
    month.list <- c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")
    month.list <- rep(month.list, each=2)
    month.list <- c(month.list[1], month.list, month.list[length(month.list)])
  }
  #need to make month.list same length as biwk 
  uni.biwk = as.numeric(as.character(unique(out.mat$biwk)))
  out.mat$month <- NA
  for (i in 1:length(uni.biwk)){
    out.mat$month[out.mat$biwk==uni.biwk[i]] = month.list[i]
  }
  
  if (mort==.207){
    out.mat$month <- factor(out.mat$month, levels=c("Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct"))
  }else if (mort==.456){
    out.mat$month <- factor(out.mat$month, levels=c("Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov"))
  }
  
  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~month, ncol=3) 
  
  ggplot(data=out.mat) + geom_point(aes(x=age, y=seroprevalence, col=biwk), pch=1) + theme_bw() + theme(panel.grid = element_blank()) +
    geom_line(aes(x=age, y=seroprevalence, col=biwk)) +  facet_wrap(~biwk, ncol=4) 
  
  #names(eq.sero.tot)[1] ="biwk"
  #plot one of them to get a sense
  #plot.age.seroprev = list.by.T[[5]]
  
  #age.sero.plot(plot.age.seroprev)
  
  #if you want to, plot all of them and save
  
  #return both seasonal and non-seasonal data
  return(out.mat)#, seas.prev))
}

#likelihood function
#this for a biweekly observation process
log.lik.wane <- function(par, data, cutoff, burnin, sim_pop, yrs, ntyr,  s, beta, sero, boost, recov, mort, mort_juv, adult_fec, wane, slope.wane, sigma, mu.sick, add.inf.mort, model){
  
  #run your model of choice
  if(model=="MSIR"){
    out = sim.met.MSIR(burnin=burnin,
                       sim_pop=sim_pop,
                       yrs=yrs,
                       ntyr=ntyr,
                       s=s,
                       beta=exp(par[2]),
                       sero=sero,
                       recov=recov,
                       mort=mort, #annual
                       mort_juv=mort_juv, #annual
                       adult_fec=adult_fec, #annual
                       wane=exp(par[1]),
                       slope.wane=slope.wane,
                       sigma=sigma,#biwks
                       model=model,
                       mu.sick = mu.sick,
                       add.inf.mort = add.inf.mort)
    
  } else if (model =="MSIRS"){
    out = sim.met.MSIR(burnin=burnin,
                       sim_pop=sim_pop,
                       yrs=yrs,
                       ntyr=ntyr,
                       s=s,
                       beta=exp(par[2]),
                       sero=sero,
                       recov=recov,
                       mort=mort, #annual
                       mort_juv=mort_juv, #annual
                       adult_fec=adult_fec, #annual
                       wane=exp(par[1]),
                       slope.wane=slope.wane,
                       sigma=exp(par[3]),
                       model=model,
                       mu.sick = mu.sick,
                       add.inf.mort = add.inf.mort)
    
  } else if (model=="MSRIR"){
    out = sim.met.MSIR(burnin=burnin,
                       sim_pop=sim_pop,
                       yrs=yrs,
                       ntyr=ntyr,
                       s=s,
                       beta=exp(par[2]),
                       sero=exp(par[3]),
                       recov=recov,
                       mort=mort, #annual
                       mort_juv=mort_juv, #annual
                       adult_fec=adult_fec, #annual
                       wane=exp(par[1]),
                       slope.wane=slope.wane,
                       sigma=sigma,
                       model=model,
                       mu.sick = mu.sick,
                       add.inf.mort = add.inf.mort)
    
  }else if (model =="MSIRN"){
    out = sim.met.MSIRN(burnin=burnin,
                        sim_pop=sim_pop,
                        yrs=yrs,
                        ntyr=ntyr,
                        s=s,
                        beta=exp(par[2]),
                        boost=boost,
                        recov=recov,
                        mort=mort, #annual
                        mort_juv=mort_juv, #annual
                        adult_fec=adult_fec, #annual
                        wane=exp(par[1]),
                        slope.wane=slope.wane,
                        sigma=exp(par[3]),
                        model=model,
                        mu.sick = mu.sick,
                        add.inf.mort = add.inf.mort)
  }else if (model =="MSIRNR"){
    out = sim.met.MSIRN(burnin=burnin,
                        sim_pop=sim_pop,
                        yrs=yrs,
                        ntyr=ntyr,
                        s=s,
                        beta=exp(par[2]),
                        recov=recov,
                        mort=mort, #annual
                        mort_juv=mort_juv, #annual
                        adult_fec=adult_fec, #annual
                        wane=exp(par[1]),
                        slope.wane=slope.wane,
                        sigma=exp(par[3]),
                        boost=exp(par[4]),
                        model=model,
                        mu.sick = mu.sick,
                        add.inf.mort = add.inf.mort)
  }
  #first, run your model at the desired parameters for estimation (par)
  #sometimes, they don't produce anything at all - in which case, just return a giant ll here and move on
  
  #print the guess par
  print(exp(par))
  
  if(is.na(unique(out$seropos)[1])){
    ll=-1*10^9#write over with real number that is huge
    print(-ll)
    return(-ll)
  }else{
    
    
    #fit a line to these data to get the hazard by age:
    out$biwk <- as.numeric(as.character(out$biwk))
    out.lst <- dlply(out, .(biwk))
    spline.list <- lapply(out.lst, spline.fit)
    
    #separate age data
    
    #data = subset(data, !is.na(age)) #already done
    data = arrange(data, age)
    #sort data by age and make a prediction for each age by biweek
    ages.new = sort(unique(data$age))
    
    pred.list = lapply(spline.list, pred.fit, ages.new=ages.new)
    #plot.ageseroprev <- function(dat){
    #with(dat, plot(age, seroprevalence, ylim=c(0,1), type="b"))
    #}
    #lapply(out.lst, plot.ageseroprev)
    
    ages.pred = rep(ages.new, times = length(pred.list))
    biwks = rep(seq(1:max(out$biwk)), each = length(ages.new ))
    predictions = c(unlist(pred.list))
    
    out.comp = cbind.data.frame(biwks, ages.pred, predictions) # the same ages as found in the dataset are represented for every biwk
    names(out.comp) = c("biwk", "age", "seroprev") 
    #gives you predicted age-seroprevalence for each biweek of the year
    out.comp$seroprev[out.comp$seroprev<0] =0
    out.comp$seroprev[out.comp$seroprev>1] =1
    
    #data already has biwk
    
    
    #slim data down to essentials
    data = dplyr::select(data, sampleID, site, date, species, sex, class, age, type, prev, prev_uci, prev_lci, biwk, fit)
    #data = subset(data, fit==1) already done
    
    
    ll=0 #sets initial log-likelihood to 0
    
    #then, compare with data by age and biweek
    for(a in 1:length(data$age)){ 
      #now run through data and compare with model
      #for-loops over every individual in the dataset, starting above the lower bound - 
      #compares to age-seroprev for those with age data and seasonal seroprev for those without
      #if (!is.na(data$age[a])){
      pprev <- out.comp$seroprev[out.comp$age==data$age[a] & out.comp$biwk==data$biwk[a]]  
      # } else if(is.na(data$age[a])){
      #  pprev <- out.seas$seroprev[out.seas$biwk==data$biwk[a]]  
      #}
      #print(pprev)
      if(cutoff=="mean"){
        ll=ll+dbinom(data$prev[a],1,pprev,log=T)    
        print(ll)
      }else if (cutoff=="lci"){
        ll=ll+dbinom(data$prev_lci[a],1,pprev,log=T)    
      }else if(cutoff=="uci"){
        ll=ll+dbinom(data$prev_uci[a],1,pprev,log=T)    
      }
      
    }
    print(-ll)
    #calculate and return the negative log-likelihood 
    return(-ll)
  }}
log.lik.wane.age <- function(par, data, cutoff, burnin, sim_pop, yrs, ntyr, age.brk,  s, beta, sero, boost, recov, mort, mort_juv, adult_fec, wane, slope.wane, sigma, mu.sick, add.inf.mort, model, N_stat){
  
  #run your model of choice
  if(model=="MSIR"){
    out = sim.met.MSIR.age(burnin=burnin,
                           sim_pop=sim_pop,
                           yrs=yrs,
                           ntyr=ntyr,
                           s=s,
                           beta=exp(par[2:length(par)]),
                           sero=sero,
                           recov=recov,
                           age.brk=age.brk,
                           mort=mort, #annual
                           mort_juv=mort_juv, #annual
                           adult_fec=adult_fec, #annual
                           wane=exp(par[1]),
                           slope.wane=slope.wane,
                           sigma=sigma,#biwks
                           model=model,
                           mu.sick = mu.sick,
                           add.inf.mort = add.inf.mort)
    
  } else if (model =="MSIRS"){
    out = sim.met.MSIR.age(burnin=burnin,
                           sim_pop=sim_pop,
                           yrs=yrs,
                           ntyr=ntyr,
                           s=s,
                           beta=exp(par[2:(length(par)-1)]),
                           sero=sero,
                           recov=recov,
                           age.brk=age.brk,
                           mort=mort, #annual
                           mort_juv=mort_juv, #annual
                           adult_fec=adult_fec, #annual
                           wane=exp(par[1]),
                           slope.wane=slope.wane,
                           sigma=exp(par[length(par)]),
                           model=model,
                           mu.sick = mu.sick,
                           add.inf.mort = add.inf.mort)
    
  } else if (model=="MSRIR"){
    out = sim.met.MSIR.age(burnin=burnin,
                           sim_pop=sim_pop,
                           yrs=yrs,
                           ntyr=ntyr,
                           s=s,
                           beta=exp(par[2:(length(par)-1)]),
                           sero=exp(par[length(par)]),
                           recov=recov,
                           age.brk=age.brk,
                           mort=mort, #annual
                           mort_juv=mort_juv, #annual
                           adult_fec=adult_fec, #annual
                           wane=exp(par[1]),
                           slope.wane=slope.wane,
                           sigma=sigma,
                           model=model,
                           mu.sick = mu.sick,
                           add.inf.mort = add.inf.mort)
    
  }else if (model =="MSIRN"){
    out = sim.met.MSIRN.age(burnin=burnin,
                            sim_pop=sim_pop,
                            yrs=yrs,
                            ntyr=ntyr,
                            s=s,
                            beta=exp(par[2:(length(par)-1)]),
                            boost=boost,
                            recov=recov,
                            mort=mort, #annual
                            mort_juv=mort_juv, #annual
                            adult_fec=adult_fec, #annual
                            wane=exp(par[1]),
                            age.brk=age.brk,
                            slope.wane=slope.wane,
                            sigma=exp(par[length(par)]),
                            model=model,
                            mu.sick = mu.sick,
                            add.inf.mort = add.inf.mort,
                            N_stat=N_stat)
  }else if (model =="MSIRNR"){
    out = sim.met.MSIRN.age(burnin=burnin,
                            sim_pop=sim_pop,
                            yrs=yrs,
                            ntyr=ntyr,
                            age.brk=age.brk,
                            s=s,
                            beta=exp(par[2:(length(par)-2)]),
                            recov=recov,
                            mort=mort, #annual
                            mort_juv=mort_juv, #annual
                            adult_fec=adult_fec, #annual
                            wane=exp(par[1]),
                            slope.wane=slope.wane,
                            sigma=exp(par[length(par)-1]),
                            boost=exp(par[length(par)]),
                            model=model,
                            mu.sick = mu.sick,
                            add.inf.mort = add.inf.mort,
                            N_stat=N_stat)
  }else if (model =="MSIRNP"){
    out = sim.met.MSIRNP.age(burnin=burnin,
                             sim_pop=sim_pop,
                             yrs=yrs,
                             ntyr=ntyr,
                             s=s,
                             beta=exp(par[2:(length(par)-3)]),
                             boost=exp(par[length(par)-1]),
                             recov=recov,
                             mort=mort, #annual
                             mort_juv=mort_juv, #annual
                             adult_fec=adult_fec, #annual
                             wane=exp(par[1]),
                             age.brk=age.brk,
                             slope.wane=slope.wane,
                             sigma=exp(par[length(par)-2]),
                             sigma2=exp(par[length(par)]),
                             model=model,
                             mu.sick = mu.sick,
                             add.inf.mort = add.inf.mort,
                             N_stat=N_stat)
  }else if (model =="MSIRNR2"){
    out = sim.met.MSIRNR2.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=exp(par[2:(length(par)-3)]),
                              boost=exp(par[length(par)-1]),
                              recov=recov,
                              mort=mort, #annual
                              mort_juv=mort_juv, #annual
                              adult_fec=adult_fec, #annual
                              wane=exp(par[1]),
                              age.brk=age.brk,
                              slope.wane=slope.wane,
                              sigma=exp(par[length(par)-2]),
                              sigma2=exp(par[length(par)]),
                              model=model,
                              mu.sick = mu.sick,
                              add.inf.mort = add.inf.mort,
                              N_stat=N_stat)
  }
  #first, run your model at the desired parameters for estimation (par)
  #sometimes, they don't produce anything at all - in which case, just return a giant ll here and move on
  
  #print the guess par
  print(exp(par))
  
  if(is.na(unique(out$seropos)[1])){
    ll=-1*10^9#write over with real number that is huge
    print(-ll)
    return(-ll)
  }else{
    
    
    #fit a line to these data to get the hazard by age:
    out$biwk <- as.numeric(as.character(out$biwk))
    out.lst <- dlply(out, .(biwk))
    spline.list <- lapply(out.lst, spline.fit)
    
    #separate age data
    
    #data = subset(data, !is.na(age)) #already done
    data = arrange(data, age)
    #sort data by age and make a prediction for each age by biweek
    ages.new = sort(unique(data$age))
    
    pred.list = lapply(spline.list, pred.fit, ages.new=ages.new)
    #plot.ageseroprev <- function(dat){
    #with(dat, plot(age, seroprevalence, ylim=c(0,1), type="b"))
    #}
    #lapply(out.lst, plot.ageseroprev)
    
    ages.pred = rep(ages.new, times = length(pred.list))
    biwks = rep(seq(1:max(out$biwk)), each = length(ages.new ))
    predictions = c(unlist(pred.list))
    
    out.comp = cbind.data.frame(biwks, ages.pred, predictions) # the same ages as found in the dataset are represented for every biwk
    names(out.comp) = c("biwk", "age", "seroprev") 
    #gives you predicted age-seroprevalence for each biweek of the year
    out.comp$seroprev[out.comp$seroprev<0] =0
    out.comp$seroprev[out.comp$seroprev>1] =1
    
    #data already has biwk
    
    
    #slim data down to essentials
    #data = dplyr::select(data, sampleID, site, date, species, sex, class, age, type, prev, prev_uci, prev_lci, biwk, fit)
    #data = subset(data, fit==1) already done
    
    
    ll=0 #sets initial log-likelihood to 0
    
    #then, compare with data by age and biweek
    for(a in 1:length(data$age)){ 
      #now run through data and compare with model
      #for-loops over every individual in the dataset, starting above the lower bound - 
      
      
      pprev <- out.comp$seroprev[out.comp$age==data$age[a] & out.comp$biwk==data$biwk[a]]  
      
      #print(pprev)
      if(cutoff=="mean"){
        ll=ll+dbinom(data$prev[a],1,pprev,log=T)    
        print(ll)
      }else if (cutoff=="lci"){
        ll=ll+dbinom(data$prev_lci[a],1,pprev,log=T)    
      }else if(cutoff=="uci"){
        ll=ll+dbinom(data$prev_uci[a],1,pprev,log=T)    
      }
      
    }
    print(-ll)
    #calculate and return the negative log-likelihood 
    return(-ll)
  }}


#need to have it specify where the age-breaks are and where the beta changes
wrap.msir.fit.wane.age <- function(model, dat1, cutoff1, species1, type1, vis_split, s, sim_pop, yrs, ntyr, recov, sero, wane, sigma, sigma2, boost, mort, mort_juv, N_stat,
                                   adult_fec, mu.sick, add.inf.mort, beta, age.brk, slope.wane, burnin, do.plot, do.save, obs.proc, filename){
  
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & type==type1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  dat.tmp <- subset(dat.tmp, fit==1)
  
  #identify your parameter(s) for fitting based on your model of choice
  if(model=="MSIR" | model=="MSI"){
    par = c(wane, beta)  
  }else if(model=="MSIRS" | model=="MSIRN"){
    par = c(wane, beta, sigma)  
  }else if (model=="MSRIR"){
    par = c(wane, beta, sero)
  }else if (model=="MSIRNR"){
    par = c(wane, beta, sigma, boost)
  } else if (model=="MSIRNP"){
    par = c(wane, beta, sigma, boost, sigma2)
  }
  print(par)
  #estimate parameters
  
  out.lik <- optim(par=log(par),fn=log.lik.wane.age, method="Nelder-Mead", data=dat.tmp, model=model,  burnin=burnin, sigma=sigma, boost=boost,
                   cutoff=cutoff1, s=s, sim_pop=sim_pop, ntyr=ntyr, yrs=yrs, wane=wane, recov=recov, slope.wane=slope.wane, sero=sero, age.brk=age.brk,
                   mort=mort, mort_juv=mort_juv, adult_fec=adult_fec, mu.sick=mu.sick, add.inf.mort=add.inf.mort, N_stat=N_stat,
                   control=list(maxit=1000))
  
  
  
  
  #write over parameters with new optimized version
  if(model=="MSIR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:length(par)])
  }else if(model=="MSIRS" | model=="MSIRN"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-1)])
    sigma = exp(out.lik$par[length(par)])
  }else if (model=="MSRIR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-1)])
    sero = exp(out.lik$par[length(par)])
  }else if (model=="MSIRNR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-2)])
    sigma = exp(out.lik$par[length(par)-1])
    boost = exp(out.lik$par[length(par)])
  } else if (model=="MSIRNP"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-3)])
    sigma = exp(out.lik$par[length(par)-2])
    boost = exp(out.lik$par[length(par)]-1)
    sigma2 = exp(out.lik$par[length(par)])
  } 
  
  #save you likelihood and convergence history
  llik=out.lik$value
  convergence = out.lik$convergence
  
  #run your model again, now with estimated, optimized par
  if(model=="MSIR" | model=="MSIRS" | model=="MSRIR" ){
    #wrap into pars and run new fitted model 
    model.out  = sim.met.MSIR.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=beta,
                                  age.brk = age.brk,
                                  sero=sero,
                                  recov=recov,
                                  mort=mort, #annual
                                  mort_juv=mort_juv, #annual
                                  adult_fec=adult_fec, #annual
                                  wane=wane, #in biweeks
                                  slope.wane=slope.wane,
                                  sigma=sigma,#biwks
                                  model=model,
                                  mu.sick = mu.sick,
                                  add.inf.mort = add.inf.mort)  
  } else if (model =="MSIRNR" | model =="MSIRN"){
    model.out  = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=beta,
                                   age.brk = age.brk,
                                   boost = boost,
                                   recov=recov,
                                   mort=mort, #annual
                                   mort_juv=mort_juv, #annual
                                   adult_fec=adult_fec, #annual
                                   wane=wane, #in biweeks
                                   slope.wane=slope.wane,
                                   sigma=sigma,#biwks
                                   model=model,
                                   mu.sick = mu.sick,
                                   add.inf.mort = add.inf.mort,
                                   N_stat=N_stat) 
    
  } else if (model =="MSIRNP"){
    model.out  = sim.met.MSIRNP.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=beta,
                                    age.brk = age.brk,
                                    boost = boost,
                                    recov=recov,
                                    mort=mort, #annual
                                    mort_juv=mort_juv, #annual
                                    adult_fec=adult_fec, #annual
                                    wane=wane, #in biweeks
                                    slope.wane=slope.wane,
                                    sigma=sigma,#biwks
                                    sigma2=sigma2,
                                    model=model,
                                    mu.sick = mu.sick,
                                    add.inf.mort = add.inf.mort,
                                    N_stat=N_stat) 
  }
  
  #now plot your model output with the data
  
  #get data in age-seroprev form
  dat.sum = get.seroprev.dat(data=dat.tmp[!is.na(dat.tmp$age),], vis_split = vis_split, cutoff = cutoff1)
  
  
  #and plot
  if(do.plot==TRUE){
    #plot with data
    max.a = ceiling(max(dat.sum$age_year))
    p1 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
      geom_line(data=model.out, aes(x=age, y=seroprevalence), color = "royalblue") + #scale_colour_manual(values=c("dodgerblue", "tomato")) +# guide= guide_legend(override.aes = list(linetype = c("blank", "blank", "solid", "solid"), shape = c(NA, NA, 1,1)))) +
      theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
      coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
      theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
    print(p1)
    
    #and also plot a version where the model gets binned over .5-year increments
    model.out$age_year <- floor(model.out$age)
    new.mod <- ddply(model.out, .(age_year), summarize, seroprevalence= mean(seroprevalence))
    #new.mod <- bin.seroprev.mod(model.out)
    p2 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
      geom_line(data=new.mod, aes(x=age_year, y=seroprevalence), color = "royalblue") + #scale_colour_manual(values=c("dodgerblue", "tomato")) +# guide= guide_legend(override.aes = list(linetype = c("blank", "blank", "solid", "solid"), shape = c(NA, NA, 1,1)))) +
      theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
      coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
      theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
    print(p2)
    
    
    if(do.save==TRUE){
      ggsave(file = filename,
             units="mm",  
             width=80, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
    
  }
  
  #return your model stats, including estimated par
  
  
  #we penalize only for the estimated parameters in AIC
  k = as.numeric(length(par))
  AIC = 2*k + 2*(as.numeric(llik))
  
  #check AICc for small sample size
  n = length(dat.tmp$age) #should be 109 so should be fine but good to check
  AICc = 2*(k)*(n/(n-k-1)) + 2*(as.numeric(llik))
  
  
  fit.stat <- data.frame(cbind(species1, type1, mort, mort_juv, adult_fec,  recov, sigma, sigma2, wane,  boost, sero, mu.sick, llik, convergence))
  names(fit.stat) <- c("species", "type", "mort", "juv_mort", "adult_fec", "recov", "sigma", "sigma2", "wane",  "boost", "sero", "mu_sick","neg_llik", "convergence")
  fit.stat$AIC = AIC
  fit.stat$AICc = AICc
  fit.stat$k = k 
  fit.stat$n = n
  fit.stat$cutoff = cutoff1
  fit.stat$model=model
  rownames(fit.stat) = c()
  
  #and the chain of betas separately
  return(list(fit.stat, beta))
}


wrap.msir.fit.wane.age.CI <- function(model, dat1, cutoff1, species1, type1, vis_split, s, sim_pop, yrs, ntyr, recov, sero, wane, sigma, sigma2, boost, mort, mort_juv, N_stat,
                                      adult_fec, mu.sick, add.inf.mort, beta, age.brk, slope.wane, burnin, do.plot, do.save, obs.proc, filename){
  
  #first pull out your data...
  dat.tmp <- subset(dat1, species==species1 & type==type1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  dat.tmp <- subset(dat.tmp, fit==1)
  
  #identify your parameter(s) for fitting based on your model of choice
  if(model=="MSIR" | model=="MSI"){
    par = c(wane, beta)  
  }else if(model=="MSIRS" | model=="MSIRN"){
    par = c(wane, beta, sigma)  
  }else if (model=="MSRIR"){
    par = c(wane, beta, sero)
  }else if (model=="MSIRNR"){
    par = c(wane, beta, sigma, boost)
  } else if (model=="MSIRNP" | model=="MSIRNR2"){
    par = c(wane, beta, sigma, boost, sigma2)
  }
  print(par)
  #estimate parameters
  
  out.lik <- optim(par=log(par),fn=log.lik.wane.age, method="Nelder-Mead", data=dat.tmp, model=model,  burnin=burnin, sigma=sigma, boost=boost,
                   cutoff=cutoff1, s=s, sim_pop=sim_pop, ntyr=ntyr, yrs=yrs, wane=wane, recov=recov, slope.wane=slope.wane, sero=sero, age.brk=age.brk,
                   mort=mort, mort_juv=mort_juv, adult_fec=adult_fec, mu.sick=mu.sick, add.inf.mort=add.inf.mort, N_stat=N_stat,
                   control=list(maxit=1000), hessian = TRUE)
  
  
  #solve the hessian to get CIs
  hess <- solve(out.lik$hessian)
  prop_sigma <-sqrt(diag(hess))
  upper<-out.lik$par+1.96*prop_sigma
  lower<-out.lik$par-1.96*prop_sigma
  ConfidenceInterval<-exp(data.frame(value=out.lik$par, upper=upper, lower=lower))
  
  ConfidenceInterval$par <- NA
  #assign par names par below
  
  
  
  #write over parameters with new optimized version
  if(model=="MSIR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:length(par)] <- "beta"
    
  }else if(model=="MSIRS" | model=="MSIRN"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-1)])
    sigma = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-1)] <- "beta"
    ConfidenceInterval$par[length(par)] <- "sigma"
    
  }else if (model=="MSRIR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-1)])
    sero = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-1)] <- "beta"
    ConfidenceInterval$par[length(par)] <- "sero"
    
  }else if (model=="MSIRNR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-2)])
    sigma = exp(out.lik$par[length(par)-1])
    boost = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-2)] <- "beta"
    ConfidenceInterval$par[length(par)-1] <- "sigma"
    ConfidenceInterval$par[length(par)] <- "boost"
    
  } else if (model=="MSIRNP"| model=="MSIRNR2"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-3)])
    sigma = exp(out.lik$par[length(par)-2])
    boost = exp(out.lik$par[length(par)]-1)
    sigma2 = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-3)] <- "beta"
    ConfidenceInterval$par[length(par)-2] <- "sigma"
    ConfidenceInterval$par[length(par)-1] <- "boost"
    ConfidenceInterval$par[length(par)] <- "sigma2"
    
  } 
  
  #save you likelihood and convergence history
  llik=out.lik$value
  convergence = out.lik$convergence
  
  #run your model again, now with estimated, optimized par
  if(model=="MSIR" | model=="MSIRS" | model=="MSRIR" ){
    #wrap into pars and run new fitted model 
    model.out  = sim.met.MSIR.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=beta,
                                  age.brk = age.brk,
                                  sero=sero,
                                  recov=recov,
                                  mort=mort, #annual
                                  mort_juv=mort_juv, #annual
                                  adult_fec=adult_fec, #annual
                                  wane=wane, #in biweeks
                                  slope.wane=slope.wane,
                                  sigma=sigma,#biwks
                                  model=model,
                                  mu.sick = mu.sick,
                                  add.inf.mort = add.inf.mort)  
  } else if (model =="MSIRNR" | model =="MSIRN"){
    model.out  = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=beta,
                                   age.brk = age.brk,
                                   boost = boost,
                                   recov=recov,
                                   mort=mort, #annual
                                   mort_juv=mort_juv, #annual
                                   adult_fec=adult_fec, #annual
                                   wane=wane, #in biweeks
                                   slope.wane=slope.wane,
                                   sigma=sigma,#biwks
                                   model=model,
                                   mu.sick = mu.sick,
                                   add.inf.mort = add.inf.mort,
                                   N_stat=N_stat) 
    
  } else if (model =="MSIRNP"){
    model.out  = sim.met.MSIRNP.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=beta,
                                    age.brk = age.brk,
                                    boost = boost,
                                    recov=recov,
                                    mort=mort, #annual
                                    mort_juv=mort_juv, #annual
                                    adult_fec=adult_fec, #annual
                                    wane=wane, #in biweeks
                                    slope.wane=slope.wane,
                                    sigma=sigma,#biwks
                                    sigma2=sigma2,
                                    model=model,
                                    mu.sick = mu.sick,
                                    add.inf.mort = add.inf.mort,
                                    N_stat=N_stat) 
  } else if (model =="MSIRNR2"){
    model.out  = sim.met.MSIRNR2.age(burnin=burnin,
                                     sim_pop=sim_pop,
                                     yrs=yrs,
                                     ntyr=ntyr,
                                     s=s,
                                     beta=beta,
                                     age.brk = age.brk,
                                     boost = boost,
                                     recov=recov,
                                     mort=mort, #annual
                                     mort_juv=mort_juv, #annual
                                     adult_fec=adult_fec, #annual
                                     wane=wane, #in biweeks
                                     slope.wane=slope.wane,
                                     sigma=sigma,#biwks
                                     sigma2=sigma2,
                                     model=model,
                                     mu.sick = mu.sick,
                                     add.inf.mort = add.inf.mort,
                                     N_stat=N_stat) 
  }
  
  #now plot your model output with the data
  
  #get data in age-seroprev form
  dat.sum = get.seroprev.dat(data=dat.tmp[!is.na(dat.tmp$age),], vis_split = vis_split, cutoff = cutoff1)
  
  
  #and plot
  if(do.plot==TRUE){
    #plot with data
    max.a = ceiling(max(dat.sum$age_year))
    p1 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
      geom_line(data=model.out, aes(x=age, y=seroprevalence), color = "royalblue") + #scale_colour_manual(values=c("dodgerblue", "tomato")) +# guide= guide_legend(override.aes = list(linetype = c("blank", "blank", "solid", "solid"), shape = c(NA, NA, 1,1)))) +
      theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
      coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
      theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
    print(p1)
    
    #and also plot a version where the model gets binned over .5-year increments
    model.out$age_year <- floor(model.out$age)
    new.mod <- ddply(model.out, .(age_year), summarize, seroprevalence= mean(seroprevalence))
    #new.mod <- bin.seroprev.mod(model.out)
    p2 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
      geom_line(data=new.mod, aes(x=age_year, y=seroprevalence), color = "royalblue") + #scale_colour_manual(values=c("dodgerblue", "tomato")) +# guide= guide_legend(override.aes = list(linetype = c("blank", "blank", "solid", "solid"), shape = c(NA, NA, 1,1)))) +
      theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
      coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
      theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
    print(p2)
    
    
    if(do.save==TRUE){
      ggsave(file = filename,
             units="mm",  
             width=80, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
    
  }
  
  #return your model stats, including estimated par
  
  
  #we penalize only for the estimated parameters in AIC
  k = as.numeric(length(par))
  AIC = 2*k + 2*(as.numeric(llik))
  
  #check AICc for small sample size
  n = length(dat.tmp$age) #should be 109 so should be fine but good to check
  AICc = 2*(k)*(n/(n-k-1)) + 2*(as.numeric(llik))
  
  
  fit.stat <- data.frame(cbind(species1, type1, mort, mort_juv, adult_fec,  recov, sigma, sigma2, wane,  boost, sero, mu.sick, llik, convergence))
  names(fit.stat) <- c("species", "type", "mort", "juv_mort", "adult_fec", "recov", "sigma", "sigma2", "wane",  "boost", "sero", "mu_sick","neg_llik", "convergence")
  fit.stat$AIC = AIC
  fit.stat$AICc = AICc
  fit.stat$k = k 
  fit.stat$n = n
  fit.stat$cutoff = cutoff1
  fit.stat$model=model
  rownames(fit.stat) = c()
  
  #and return your confidence intervals, as estimated by the hessian
  
  
  #and the chain of betas separately
  return(list(fit.stat, beta, ConfidenceInterval))
}


wrap.msir.fit.wane.age.CI.method2 <- function(model, dat1, cutoff1, species1, type1, vis_split, s, sim_pop, yrs, ntyr, recov, sero, wane, sigma, sigma2, boost, mort, mort_juv, N_stat,
                                              adult_fec, mu.sick, add.inf.mort, beta, age.brk, slope.wane, burnin, do.plot, do.save, obs.proc, filename){
  
  #firstz pull out your data...
  dat.tmp <- subset(dat1, species==species1 & type==type1)
  dat.tmp <- subset(dat.tmp, !is.na(age))
  #dat.tmp <- subset(dat.tmp, fit==1). do you use all the data or just the longitudinal resamples?
  
  #identify your parameter(s) for fitting based on your model of choice
  if(model=="MSIR" | model=="MSI"){
    par = c(wane, beta)  
  }else if(model=="MSIRS" | model=="MSIRN"){
    par = c(wane, beta, sigma)  
  }else if (model=="MSRIR"){
    par = c(wane, beta, sero)
  }else if (model=="MSIRNR"){
    par = c(wane, beta, sigma, boost)
  } else if (model=="MSIRNP" | model=="MSIRNR2"){
    par = c(wane, beta, sigma, boost, sigma2)
  }
  print(par)
  #estimate parameters
  
  out.lik <- optim(par=log(par),fn=log.lik.wane.age, method="Nelder-Mead", data=dat.tmp, model=model,  burnin=burnin, sigma=sigma, boost=boost,
                   cutoff=cutoff1, s=s, sim_pop=sim_pop, ntyr=ntyr, yrs=yrs, wane=wane, recov=recov, slope.wane=slope.wane, sero=sero, age.brk=age.brk,
                   mort=mort, mort_juv=mort_juv, adult_fec=adult_fec, mu.sick=mu.sick, add.inf.mort=add.inf.mort, N_stat=N_stat,
                   control=list(maxit=1000), hessian = TRUE)
  
  #save the raw CI values from the hessian (not including variance)
  hess <- solve(out.lik$hessian)
  prop_sigma <-sqrt(diag(hess))
  upper<-exp(out.lik$par)+1.96*prop_sigma
  lower<-exp(out.lik$par)-1.96*prop_sigma
  ConfidenceInterval<-data.frame(value=exp(out.lik$par), upper=upper, lower=lower)
  
  ConfidenceInterval[  ConfidenceInterval<0] <- 0
  
  ConfidenceInterval$par <- NA
  #assign par names par below
  
  
  
  #write over parameters with new optimized version
  if(model=="MSIR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:length(par)])
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:length(par)] <- "beta"
    
  }else if(model=="MSIRS" | model=="MSIRN"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-1)])
    sigma = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-1)] <- "beta"
    ConfidenceInterval$par[length(par)] <- "sigma"
    
  }else if (model=="MSRIR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-1)])
    sero = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-1)] <- "beta"
    ConfidenceInterval$par[length(par)] <- "sero"
    
  }else if (model=="MSIRNR"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-2)])
    sigma = exp(out.lik$par[length(par)-1])
    boost = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-2)] <- "beta"
    ConfidenceInterval$par[length(par)-1] <- "sigma"
    ConfidenceInterval$par[length(par)] <- "boost"
    
  } else if (model=="MSIRNP"| model=="MSIRNR2"){
    wane = exp(out.lik$par[1])
    beta = exp(out.lik$par[2:(length(par)-3)])
    sigma = exp(out.lik$par[length(par)-2])
    boost = exp(out.lik$par[length(par)]-1)
    sigma2 = exp(out.lik$par[length(par)])
    
    ConfidenceInterval$par[1] <- "wane"
    ConfidenceInterval$par[2:(length(par)-3)] <- "beta"
    ConfidenceInterval$par[length(par)-2] <- "sigma"
    ConfidenceInterval$par[length(par)-1] <- "boost"
    ConfidenceInterval$par[length(par)] <- "sigma2"
    
  } 
  
  #save you likelihood and convergence history
  llik=out.lik$value
  convergence = out.lik$convergence
  
  
  #solve the hessian to get CIs using mvnorm
  #sigma is the covariance matrix which is the inverse of the hessian when you are minimizing the neg log liklihood in your function
  par.dist <- rmvnorm(n=500, exp(out.lik$par), sigma=(solve(out.lik$hessian)))
  par.dist[ par.dist < 0] <- 0
  
  #we can't have negative pars, so replace those with 0s
  
  #now have 500 guesses for each of our 3 pars, all calculated in tandem.
  
  ## get predicts across this span - including variance covariance
  ## span should account for biwk of year and age across 1 year and the full age range of the model
  ## now have a matrix to hold output for all ages and biweeks
  age.span = rep(0:(s-1), each = ntyr)
  biwk.span = rep(1:ntyr, each = s)
  
  #drop the last term in each
  age.span <- age.span[1:(length(age.span)-1)]
  biwk.span <- biwk.span[1:(length(biwk.span)-1)]
  
  preds <- matrix(NA,nrow=length(biwk.span),ncol=length(par.dist[,1]))
  
  
  
  # for each, run the model to get seroprev with the fitted version, 
  # then run it 500 times to get CIs on seroprev from the variance
  
  #make your CI terms (fill them in later)
  beta_lci <- beta_uci <- wane_lci <- wane_uci <- sigma_lci <- sigma_uci <- sigma2_lci <- sigma2_uci <- sero_lci <- sero_uci <- boost_lci <- boost_uci <- 0
  
  if (model=="MSIR"){
    
    beta_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="beta"]
    beta_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="beta"]
    
    wane_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="wane"]
    wane_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="wane"]
    
    #run it to get the mean
    model.out  = sim.met.MSIR.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=beta,
                                  age.brk = age.brk,
                                  sero=sero,
                                  recov=recov,
                                  mort=mort, 
                                  mort_juv=mort_juv, 
                                  adult_fec=adult_fec, 
                                  wane=wane, 
                                  slope.wane=slope.wane,
                                  sigma=sigma,
                                  model=model,
                                  mu.sick = mu.sick,
                                  add.inf.mort = add.inf.mort)  
    
    #get CIs one way
    for (j in 1:length(par.dist[,1])){
      
      #and the CIs
      tmp <- sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.dist[j,2],
                              age.brk = age.brk,
                              sero=sero,
                              recov=recov,
                              mort=mort, #annual
                              mort_juv=mort_juv, #annual
                              adult_fec=adult_fec, #annual
                              wane=par.dist[j,1], #in biweeks
                              slope.wane=slope.wane,
                              sigma=sigma,#biwks
                              model=model,
                              mu.sick = mu.sick,
                              add.inf.mort = add.inf.mort)
      
      #each column is a new guess. each row is output by age and biweek
      preds[,j] <- tmp$seroprevalence
      
    }
    
    
    
  }else if (model=="MSRIR"){
    
    beta_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="beta"]
    beta_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="beta"]
    
    wane_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="wane"]
    wane_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="wane"]
    
    sero_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sero"]
    sero_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sero"]
    
    
    model.out  = sim.met.MSIR.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=beta,
                                  age.brk = age.brk,
                                  sero=sero,
                                  recov=recov,
                                  mort=mort, 
                                  mort_juv=mort_juv, 
                                  adult_fec=adult_fec, 
                                  wane=wane, 
                                  slope.wane=slope.wane,
                                  sigma=sigma,
                                  model=model,
                                  mu.sick = mu.sick,
                                  add.inf.mort = add.inf.mort) 
    
    
    for (j in 1:length(par.dist[,1])){
      
      tmp <- sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.dist[j,2],
                              age.brk = age.brk,
                              sero=par.dist[j,3],
                              recov=recov,
                              mort=mort, #annual
                              mort_juv=mort_juv, #annual
                              adult_fec=adult_fec, #annual
                              wane=par.dist[j,1], #in biweeks
                              slope.wane=slope.wane,
                              sigma=sigma,#biwks
                              model=model,
                              mu.sick = mu.sick,
                              add.inf.mort = add.inf.mort)
      
      #each column is a new guess. each row is output by age and biweek
      preds[,j] <- tmp$seroprevalence
      
    }
    
    
    
    
    
  }else if (model=="MSIRS"){
    
    beta_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="beta"]
    beta_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="beta"]
    
    wane_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="wane"]
    wane_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="wane"]
    
    sigma_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sigma"]
    sigma_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sigma"]
    
    
    model.out  = sim.met.MSIR.age(burnin=burnin,
                                  sim_pop=sim_pop,
                                  yrs=yrs,
                                  ntyr=ntyr,
                                  s=s,
                                  beta=beta,
                                  age.brk = age.brk,
                                  sero=sero,
                                  recov=recov,
                                  mort=mort, 
                                  mort_juv=mort_juv, 
                                  adult_fec=adult_fec, 
                                  wane=wane, 
                                  slope.wane=slope.wane,
                                  sigma=sigma,
                                  model=model,
                                  mu.sick = mu.sick,
                                  add.inf.mort = add.inf.mort) 
    
    for (j in 1:length(par.dist[,1])) {
      
      tmp <- sim.met.MSIR.age(burnin=burnin,
                              sim_pop=sim_pop,
                              yrs=yrs,
                              ntyr=ntyr,
                              s=s,
                              beta=par.dist[j,2],
                              age.brk = age.brk,
                              sero=sero,
                              recov=recov,
                              mort=mort, #annual
                              mort_juv=mort_juv, #annual
                              adult_fec=adult_fec, #annual
                              wane=par.dist[j,1], #in biweeks
                              slope.wane=slope.wane,
                              sigma=par.dist[j,3],#biwks
                              model=model,
                              mu.sick = mu.sick,
                              add.inf.mort = add.inf.mort)
      
      #each column is a new guess. each row is output by age and biweek
      preds[,j] <- tmp$seroprevalence
    }
    
    
    
  }else if (model =="MSIRN"){
    
    beta_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="beta"]
    beta_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="beta"]
    
    wane_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="wane"]
    wane_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="wane"]
    
    sigma_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sigma"]
    sigma_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sigma"]
    
    
    model.out  = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=beta,
                                   age.brk = age.brk,
                                   boost = boost,
                                   recov=recov,
                                   mort=mort, 
                                   mort_juv=mort_juv,
                                   adult_fec=adult_fec, 
                                   wane=wane,
                                   slope.wane=slope.wane,
                                   sigma=sigma,
                                   model=model,
                                   mu.sick = mu.sick,
                                   add.inf.mort = add.inf.mort,
                                   N_stat=N_stat) 
    
    
    for (j in 1:length(par.dist[,1])) {
      
      tmp  = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,2],
                               age.brk = age.brk,
                               boost = boost,
                               recov=recov,
                               mort=mort, 
                               mort_juv=mort_juv,
                               adult_fec=adult_fec,
                               wane=par.dist[j,1], 
                               slope.wane=slope.wane,
                               sigma=par.dist[j,3],
                               model=model,
                               mu.sick = mu.sick,
                               add.inf.mort = add.inf.mort,
                               N_stat=N_stat) 
      
      preds[,j] <- tmp$seroprevalence
    }
    
    
    
  }else if (model =="MSIRNR"){
    
    beta_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="beta"]
    beta_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="beta"]
    
    wane_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="wane"]
    wane_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="wane"]
    
    sigma_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sigma"]
    sigma_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sigma"]
    
    boost_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="boost"]
    boost_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="boost"]
    
    model.out  = sim.met.MSIRN.age(burnin=burnin,
                                   sim_pop=sim_pop,
                                   yrs=yrs,
                                   ntyr=ntyr,
                                   s=s,
                                   beta=beta,
                                   age.brk = age.brk,
                                   boost = boost,
                                   recov=recov,
                                   mort=mort, 
                                   mort_juv=mort_juv,
                                   adult_fec=adult_fec, 
                                   wane=wane,
                                   slope.wane=slope.wane,
                                   sigma=sigma,
                                   model=model,
                                   mu.sick = mu.sick,
                                   add.inf.mort = add.inf.mort,
                                   N_stat=N_stat) 
    
    for (j in 1:length(par.dist[,1])) {
      
      tmp  = sim.met.MSIRN.age(burnin=burnin,
                               sim_pop=sim_pop,
                               yrs=yrs,
                               ntyr=ntyr,
                               s=s,
                               beta=par.dist[j,2],
                               age.brk = age.brk,
                               boost = par.dist[j,4],
                               recov=recov,
                               mort=mort, #annual
                               mort_juv=mort_juv, #annual
                               adult_fec=adult_fec, #annual
                               wane=par.dist[j,1], #in biweeks
                               slope.wane=slope.wane,
                               sigma=par.dist[j,3],#biwks
                               model=model,
                               mu.sick = mu.sick,
                               add.inf.mort = add.inf.mort,
                               N_stat=N_stat) 
      
      preds[,j] <- tmp$seroprevalence
    }
    
    
    
  }else if (model =="MSIRNP"){
    
    beta_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="beta"]
    beta_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="beta"]
    
    wane_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="wane"]
    wane_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="wane"]
    
    sigma_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sigma"]
    sigma_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sigma"]
    
    boost_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="boost"]
    boost_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="boost"]
    
    sigma2_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sigma2"]
    sigma2_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sigma2"]
    
    model.out  = sim.met.MSIRNP.age(burnin=burnin,
                                    sim_pop=sim_pop,
                                    yrs=yrs,
                                    ntyr=ntyr,
                                    s=s,
                                    beta=beta,
                                    age.brk = age.brk,
                                    boost = boost,
                                    recov=recov,
                                    mort=mort, 
                                    mort_juv=mort_juv, 
                                    adult_fec=adult_fec, 
                                    wane=wane, 
                                    slope.wane=slope.wane,
                                    sigma=sigma,
                                    sigma2=sigma2,
                                    model=model,
                                    mu.sick = mu.sick,
                                    add.inf.mort = add.inf.mort,
                                    N_stat=N_stat) 
    
    for (j in 1:length(par.dist[,1])) {
      
      tmp  = sim.met.MSIRNP.age(burnin=burnin,
                                sim_pop=sim_pop,
                                yrs=yrs,
                                ntyr=ntyr,
                                s=s,
                                beta=par.dist[j,2],
                                age.brk = age.brk,
                                boost = par.dist[j,4],
                                recov=recov,
                                mort=mort, #annual
                                mort_juv=mort_juv, #annual
                                adult_fec=adult_fec, #annual
                                wane=par.dist[j,1], #in biweeks
                                slope.wane=slope.wane,
                                sigma=par.dist[j,3],#biwks
                                sigma2=par.dist[j,5],
                                model=model,
                                mu.sick = mu.sick,
                                add.inf.mort = add.inf.mort,
                                N_stat=N_stat) 
      
      preds[,j] <- tmp$seroprevalence
    }
    
    
  }else if (model =="MSIRNR2"){
    
    beta_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="beta"]
    beta_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="beta"]
    
    wane_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="wane"]
    wane_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="wane"]
    
    sigma_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sigma"]
    sigma_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sigma"]
    
    boost_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="boost"]
    boost_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="boost"]
    
    sigma2_lci <- ConfidenceInterval$lower[ConfidenceInterval$par=="sigma2"]
    sigma2_uci <- ConfidenceInterval$upper[ConfidenceInterval$par=="sigma2"]
    
    
    model.out  = sim.met.MSIRNR2.age(burnin=burnin,
                                     sim_pop=sim_pop,
                                     yrs=yrs,
                                     ntyr=ntyr,
                                     s=s,
                                     beta=beta,
                                     age.brk = age.brk,
                                     boost = boost,
                                     recov=recov,
                                     mort=mort, 
                                     mort_juv=mort_juv, 
                                     adult_fec=adult_fec, 
                                     wane=wane, 
                                     slope.wane=slope.wane,
                                     sigma=sigma,
                                     sigma2=sigma2,
                                     model=model,
                                     mu.sick = mu.sick,
                                     add.inf.mort = add.inf.mort,
                                     N_stat=N_stat) 
    
    for (j in 1:length(par.dist[,1])) {
      
      tmp  = sim.met.MSIRNR2.age(burnin=burnin,
                                 sim_pop=sim_pop,
                                 yrs=yrs,
                                 ntyr=ntyr,
                                 s=s,
                                 beta=par.dist[j,2],
                                 age.brk = age.brk,
                                 boost = par.dist[j,4],
                                 recov=recov,
                                 mort=mort, #annual
                                 mort_juv=mort_juv, #annual
                                 adult_fec=adult_fec, #annual
                                 wane=par.dist[j,1], #in biweeks
                                 slope.wane=slope.wane,
                                 sigma=par.dist[j,3],#biwks
                                 sigma2=par.dist[j,5],
                                 model=model,
                                 mu.sick = mu.sick,
                                 add.inf.mort = add.inf.mort,
                                 N_stat=N_stat) 
      
      preds[,j] <- tmp$seroprevalence
    }
    
    
  }
  
  #then, we process the cumulative fits to get the CI
  #replace all negative vals with 0
  
  preds[preds<0] <- 0
  
  pred.qs  <- apply(X=preds,MARGIN=1,FUN=quantile,probs=c(0.025,0.975), na.rm=TRUE) 
  
  #and combine to produce the time series:
  model.out$seroprev_lci = pred.qs[1,]
  model.out$seroprev_uci = pred.qs[2,]
  
  
  #now plot your model output with the data
  
  #get data in age-seroprev form
  dat.sum = get.seroprev.dat(data=dat.tmp[!is.na(dat.tmp$age),], vis_split = vis_split, cutoff = cutoff1)
  
  
  #and plot
  if(do.plot==TRUE){
    #plot with data
    max.a = ceiling(max(dat.sum$age_year))
    p1 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
      geom_line(data=model.out, aes(x=age, y=seroprevalence), color = "royalblue") + 
      geom_ribbon(data=model.out, aes(x=age, ymin=seroprev_lci, ymax=seroprev_uci), fill = "royalblue", alpha =.3) +
      theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
      coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
      theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
    print(p1)
    
    #and also plot a version where the model gets binned over .5-year increments
    model.out$age_year <- floor(model.out$age)
    new.mod <- ddply(model.out, .(age_year), summarize, seroprevalence= mean(seroprevalence), seroprev_lci=mean(seroprev_lci), seroprev_uci=mean(seroprev_uci), seroprev_lci2=mean(seroprev_lci2), seroprev_uci2=mean(seroprev_uci2))
    
    p2 <- ggplot() + geom_point(data=dat.sum, aes(x=age_year, y=prevalence, size=count), color ="black", pch=1, alpha=.9) + theme_bw() +
      geom_line(data=new.mod, aes(x=age_year, y=seroprevalence), color = "royalblue") + 
      geom_ribbon(data=new.mod, aes(x=age_year, ymin=seroprev_lci, ymax=seroprev_uci), fill = "royalblue", alpha =.3) +
      theme(legend.position = c(.27,.9),  legend.title=element_blank(), legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=10), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
      coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
      theme(axis.title = element_text(size=11), axis.text = element_text(size=10), panel.grid = element_blank())
    print(p2)
    
    
    if(do.save==TRUE){
      ggsave(file = filename,
             units="mm",  
             width=80, 
             height=60, 
             scale=3, 
             dpi=300)
    }
    
    
  }
  
  #return your model stats, including estimated pars, the lci/uci and the output seroprev with CIs
  
  
  #we penalize only for the estimated parameters in AIC
  k = as.numeric(length(par))
  AIC = 2*k + 2*(as.numeric(llik))
  
  #check AICc for small sample size
  n = length(dat.tmp$age) #should be 109 so should be fine but good to check
  AICc = 2*(k)*(n/(n-k-1)) + 2*(as.numeric(llik))
  
  
  
  fit.stat <- data.frame(cbind(species1, type1, mort, mort_juv, adult_fec,  recov, beta, beta_lci, beta_uci, sigma, sigma_lci, sigma_uci, sigma2, sigma2_lci, sigma2_uci, wane, wane_lci, wane_uci, boost, boost_lci, boost_uci, sero, sero_lci, sero_uci, mu.sick, llik, convergence))
  names(fit.stat) <- c("species", "type", "mort", "juv_mort", "adult_fec", "recov", "beta", "beta_lci", "beta_uci", "sigma", "sigma_lci", "sigma_uci", "sigma2", "sigma2_lci", "sigma2_uci", "wane", "wane_lci", "wane_uci", "boost", "boost_lci", "boost_uci", "sero", "sero_lci", "sero_uci", "mu_sick","neg_llik", "convergence")
  fit.stat$AIC = AIC
  fit.stat$AICc = AICc
  fit.stat$k = k 
  fit.stat$n = n
  fit.stat$cutoff = cutoff1
  fit.stat$model=model
  rownames(fit.stat) = c()
  
  #attach identifying info to model output
  
  model.out$model = model
  model.out$cutoff = cutoff1
  model.out$species = species1
  model.out$type = type1
  
  
  #then, if need be, add any identifying information to the model output as well, 
  #and return them both
  
  
  
  return(list(fit.stat, model.out))
}


#simulate model

out.sim <- sim.met.MSIRN.age(burnin=20,
                             sim_pop=10000,
                             yrs=50,
                             ntyr=26,
                             s=20,
                             beta = 2.203968,
                             age.brk = 20,
                             #sero=0,
                             recov=1, 
                             mort=.207,
                             mort_juv=0.456,
                             adult_fec=.48,
                             wane=0.0850142487956748,
                             slope.wane=1,
                             sigma=0.00748049787194744,
                             #sigma2= 0,
                             boost=0, 
                             mu.sick = 1,
                             add.inf.mort = FALSE, 
                             N_stat="matAB",
                             model = "MSIRN")


#eidolon parameters:
# mort = .207; mort_juv= .456; adult_fec = .48
#pteropus parameters:
# mort = .489; mort_juv= .456; adult_fec = .48

MSIRN.EID.mean.NiV <- wrap.msir.fit.wane.age.CI.method2(model="MSIRN",
                                                                   dat1=dat,
                                                                   species1 = "Eidolon dupreanum",
                                                                   type1="NiV",
                                                                   cutoff1 = "mean",
                                                                   vis_split = 3,
                                                                   burnin=20,
                                                                   sim_pop=10000,
                                                                   yrs=50,
                                                                   ntyr=26,
                                                                   s=20,
                                                                   beta = 2.203968,
                                                                   age.brk = 20,
                                                                   sero=0,
                                                                   recov=1, 
                                                                   mort=.207,
                                                                   mort_juv=0.456,
                                                                   adult_fec=.48,
                                                                   wane=0.0850142487956748,
                                                                   slope.wane=1,
                                                                   sigma=0.00748049787194744,
                                                                   sigma2= 0,
                                                                   boost=0, 
                                                                   mu.sick = 1,
                                                                   add.inf.mort = FALSE, 
                                                                   N_stat="matAB",
                                                                   do.plot=FALSE,
                                                                   do.save=FALSE,
                                                                   filename = NA)

save(MSIRN.EID.mean.NiV, file=paste0(homewd, "/data/MSIRN.EID.mean.NiV.Rdata"))
