#remove background
library(plyr)
library(dplyr)
library(ggplot2)

rm(list-ls())

#read in csv
dat = read.csv(file="WN.csv", header = T, stringsAsFactors = F)
head(dat)

#subset to age and infection status, by species
dat.sub <- dplyr::select(dat, Genre.et.especes, A, SX, Res)
head(dat.sub)
names(dat.sub) <- c("species", "class", "sex", "res")
unique(dat.sub$sex) #what is I?
unique(dat.sub$class)
#sum by species and age
dat.sum <-ddply(dat.sub, .(species, class),summarize, N_pos=sum(res), N=length(res))
dat.sum$seroprevalence = dat.sum$N_pos/dat.sum$N

#and plot
ggplot(data=dat.sum) + geom_point(aes(x=species, y=seroprevalence, color=class, size=N))

dat.slim = subset(dat.sum, N>9)

#too many species to see. just look at those species where you caught at least 10
ggplot(data=dat.slim) + geom_point(aes(x=species, y=seroprevalence, color=class, size=N))

#look at age-seroprev across all species summed
dat.tot= ddply(dat.sum, .(class), summarize, N_pos=sum(N_pos), N=sum(N))
dat.tot$seroprevalence = dat.tot$N_pos/dat.tot$N
  
ggplot(data=dat.tot) + geom_point(aes(x=class, y=seroprevalence, color=class, size=N))
