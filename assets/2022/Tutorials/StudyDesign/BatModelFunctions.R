
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
sim.model <- function(par, model){
  
  #write down the par as inputs into each model
  if(model=="MSIRN"){
    model.sim <- sim.met.MSIRN.age(burnin=par$burnin,
                                   sim_pop=par$sim_pop,
                                   yrs=par$yrs,
                                   ntyr=par$ntyr,
                                   s=par$s,
                                   beta = par$beta,
                                   age.brk = par$age.brk,
                                   #sero=0,
                                   recov=par$recov, 
                                   mort=par$mort,
                                   mort_juv=par$mort_juv,
                                   adult_fec=par$adult_fec,
                                   wane=par$wane,
                                   slope.wane=par$slope.wane,
                                   sigma=par$sigma,
                                   #sigma2= 0,
                                   boost=par$boost, 
                                   mu.sick = par$mu.sick,
                                   add.inf.mort = par$add.inf.mort, 
                                   N_stat=par$N_stat,
                                   model = par$model)
    
  }else if(model=="MSIRS"){
    model.sim  = sim.met.MSIR.age(burnin=par$burnin,
                                  sim_pop=par$sim_pop,
                                  yrs=par$yrs,
                                  ntyr=par$ntyr,
                                  s=par$s,
                                  beta=par$beta,
                                  age.brk = par$age.brk,
                                  sero=par$sero,
                                  recov=par$recov,
                                  mort=par$mort, 
                                  mort_juv=par$mort_juv, 
                                  adult_fec=par$adult_fec, 
                                  wane=par$wane, 
                                  slope.wane=par$slope.wane,
                                  sigma=par$sigma,
                                  model=par$model,
                                  mu.sick = par$mu.sick,
                                  add.inf.mort = par$add.inf.mort) 
    
  }
  
  return(model.sim)
}
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
  #print(paste0("lambda = ", lambda)) #we need the pop to replace itself or grow slightly - must be 1 or greater
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
  
  
  
  
  dat.tot = data.frame(rbind(cbind(dat.tot[,1],dat.tot[,2]), cbind(dat.tot[,1], dat.tot[,3]), cbind(dat.tot[,1], dat.tot[,4]), cbind(dat.tot[,1], dat.tot[,5])))
  names(dat.tot) = c("time", "count")
  dat.tot$class = rep(c("M", "S", "I", "R"), each= length(times))
  dat.tot$class = factor(dat.tot$class, levels=c("M", "S", "I", "R"))
  
  
  dat.tot$proportion = dat.tot$count/sim_pop
  
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
  
  #now, for simplicity, collapse to just age-seroprevalence cross-sectionally
  out.mat$age <- trunc(out.mat$age)
  #and because you can't have partial bats
  out.mat$seropos <- trunc(out.mat$seropos)
  out.mat$N <- trunc(out.mat$N)
  age.ser.dat <- ddply(out.mat, .(age), summarise, seropos = sum(seropos), N= sum(N))
  age.ser.dat$seroprevalence <- age.ser.dat$seropos/age.ser.dat$N
  
  ggplot(age.ser.dat) + geom_point(aes(x=age, y=seroprevalence, size=N)) +
    geom_line(aes(x=age, y=seroprevalence))                             
  
  
  
  #return both seasonal and non-seasonal data
  #return(out.mat)#, seas.prev))
  return(age.ser.dat)
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
  
  prop.tot = dat.tot
  prop.tot$M = prop.tot$M/prop.tot$pop
  prop.tot$S = prop.tot$S/prop.tot$pop
  prop.tot$I = prop.tot$I/prop.tot$pop
  prop.tot$R = prop.tot$R/prop.tot$pop
  prop.tot$N = prop.tot$N/prop.tot$pop
  
  
  
  
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
  
  #eq.dat = subset(age.dat.tot, times >=(5) & times <=(6))
  eq.dat = subset(age.dat.tot, times >=(yrs-1))
  
  
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
  
  
  
  #return both seasonal and non-seasonal data
  #return(out.mat)#, seas.prev))
  return(age.ser.dat)
}
fit.mod.dat <- function(model, dat1, par.dat, do.plot, do.save, obs.proc, filename){
  
  #build data
  age.split <- dlply(dat1, .(age))
  
  age.df <- lapply(age.split, build.data)
  
  age.df <- data.table::rbindlist(age.df)
  
  #identify your parameter(s) for fitting based on your model of choice
  
  # if(model=="MSIR" | model=="MSI"){
  #   par = c(wane, beta)  
  # }else if(model=="MSIRS" | model=="MSIRN"){
  #   par = c(beta)  
  # }else if (model=="MSRIR"){
  #   par = c(wane, beta, sero)
  # }else if (model=="MSIRNR"){
  #   par = c(wane, beta, sigma, boost)
  # } else if (model=="MSIRNP" | model=="MSIRNR2"){
  #   par = c(wane, beta, sigma, boost, sigma2)
  # }
  # print(par)
  #estimate parameters
  
  #for MSIRS and MSIRN, par here is just beta
  par = par.dat$beta
  
  out.lik <- optim(par=log(par),fn=log.lik.wane.age, method="Nelder-Mead", model= model,
                   data=age.df, par.dat=par.dat,
                   control=list(maxit=1000), hessian = TRUE)
  
  
  #now run your model output with these new data
  par.dat$beta <- exp(out.lik$par)
  
  
  
  model.sim <- sim.model(par = par.dat,
                         model = model)
  
  
  
  #and plot with the data
  
  #we removed the confidence interval section here to speed things up
  
  if(do.plot==TRUE){
    #plot with data
    max.a = ceiling(max(dat1$age))
    p1 <- ggplot() + geom_point(data=dat1, aes(x=age, y=seroprevalence, size=N), color ="black", alpha=.9) + theme_bw() +
      geom_line(data=model.sim, aes(x=age, y=seroprevalence), linewidth= 1, color = "royalblue") + 
      theme(  legend.spacing = unit(.03,"cm"),  legend.text=element_text(size=12), legend.key.size = unit(.4, "cm"), legend.box = "horizontal") +
      coord_cartesian(ylim=c(0,1), xlim=c(0, max.a)) + ylab("seroprevalence") + xlab("age(yrs)")  + guides(colour = guide_legend(order = 2), size = guide_legend(order = 1)) +
      theme(axis.title = element_text(size=18), axis.text = element_text(size=14), panel.grid = element_blank())
    print(p1)
    
    
    
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
  AIC = 2*k + 2*(as.numeric(out.lik$value))
  
  #check AICc for small sample size
  n = length(age.df$age) #should be field_sample
  AICc = 2*(k)*(n/(n-k-1)) + 2*(as.numeric(out.lik$value))
  
  
  #fit.stat <- data.frame(cbind(species, type, mort, mort_juv, adult_fec,  recov, beta, beta_lci, beta_uci, sigma, sigma_lci, sigma_uci, sigma2, sigma2_lci, sigma2_uci, wane, wane_lci, wane_uci, boost, boost_lci, boost_uci, sero, sero_lci, sero_uci, mu.sick, llik, convergence))
  #names(fit.stat) <- c("species", "type", "mort", "juv_mort", "adult_fec", "recov", "beta", "beta_lci", "beta_uci", "sigma", "sigma_lci", "sigma_uci", "sigma2", "sigma2_lci", "sigma2_uci", "wane", "wane_lci", "wane_uci", "boost", "boost_lci", "boost_uci", "sero", "sero_lci", "sero_uci", "mu_sick","neg_llik", "convergence")
  fit.stat <- par.dat
  fit.stat$beta <- exp(out.lik$par)
  fit.stat$value <- out.lik$value
  fit.stat$convergence <- out.lik$convergence
  
  
  fit.stat$AIC = AIC
  fit.stat$AICc = AICc
  fit.stat$k = k 
  fit.stat$n = n
  #fit.stat$cutoff = cutoff1
  fit.stat$model=model
  rownames(fit.stat) = c()
  
  #attach identifying info to model output
  
  model.sim$model = model
  
  
  #then, if need be, add any identifying information to the model output as well, 
  #and return them both
  #return(list(fit.stat, model.out))
  
  #return only your comparative pars
  return(fit.stat)
}
build.data <- function(df){
  
  dat = cbind.data.frame(age = rep(df$age, df$N), sero = 0)
  dat$index = rownames(dat)
  #randomly change a subset of seronegs to seropos
  if(df$seropos>0){
    seropos.ind <- sample(dat$index, size = df$seropos, replace = F)  
    dat[seropos.ind,]$sero <- 1
  }
  
  
  
  return(dat)
}
sub.sample <- function(sim.dat, field_sample){
  #now convert to an actual data table, where each bat is given an age, and a serostatus
  age.split <- dlply(sim.dat, .(age))
  
  
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
  return(sub.df)
  
}
log.lik.wane.age <- function(par, data, model, par.dat){
  
  #run your model of choice
  if (model =="MSIRS"){
    out = sim.met.MSIR.age(burnin=par.dat$burnin,
                           sim_pop=par.dat$sim_pop,
                           yrs=par.dat$yrs,
                           ntyr=par.dat$ntyr,
                           s=par.dat$s,
                           beta=exp(par),
                           sero=par.dat$sero,
                           recov=par.dat$recov,
                           age.brk=par.dat$age.brk,
                           mort=par.dat$mort, #annual
                           mort_juv=par.dat$mort_juv, #annual
                           adult_fec=par.dat$adult_fec, #annual
                           wane=par.dat$wane,
                           slope.wane=par.dat$slope.wane,
                           sigma=par.dat$sigma,
                           model=par.dat$model,
                           mu.sick = par.dat$mu.sick,
                           add.inf.mort = par.dat$add.inf.mort)
    
    
  }else if (model =="MSIRN"){
    out = sim.met.MSIRN.age(burnin=par.dat$burnin,
                            sim_pop=par.dat$sim_pop,
                            yrs=par.dat$yrs,
                            ntyr=par.dat$ntyr,
                            s=par.dat$s,
                            beta=exp(par),
                            boost=par.dat$boost,
                            recov=par.dat$recov,
                            mort=par.dat$mort, #annual
                            mort_juv=par.dat$mort_juv, #annual
                            adult_fec=par.dat$adult_fec, #annual
                            wane=par.dat$wane,
                            age.brk=par.dat$age.brk,
                            slope.wane=par.dat$slope.wane,
                            sigma=par.dat$sigma,
                            model=par.dat$model,
                            mu.sick = par.dat$mu.sick,
                            add.inf.mort = par.dat$add.inf.mort,
                            N_stat=par.dat$N_stat)
  }
  #first, run your model at the desired parameters for estimation (par)
  #sometimes, they don't produce anything at all - in which case, just return a giant ll here and move on
  
  #print the guess par
  #print(exp(par))
  
  if(is.na(unique(out$seropos)[1])){
    ll=-1*10^9#write over with real number that is huge
    print(-ll)
    return(-ll)
  }else{
    
    
    
    data = arrange(data, age)
    
    
    ll=0 #sets initial log-likelihood to 0
    
    #then, compare with data by age and biweek
    for(a in 1:length(data$age)){ 
      #now run through data and compare with model
      #for-loops over every individual in the dataset, starting above the lower bound - 
      
      
      pprev <- out$seroprevalence[out$age==data$age[a]]  
      
      #print(pprev)
      
      ll=ll+dbinom(data$sero[a],1,pprev,log=T)    
      
      
    }
    print(-ll)
    #calculate and return the negative log-likelihood 
    return(-ll)
  }}
