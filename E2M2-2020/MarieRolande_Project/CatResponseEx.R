rm(list=ls())

library(plyr)
library(dplyr)
library(MASS) #allows you to do ordinal logistic regression

dat <- read.csv(file="dat1.csv", header = T, stringsAsFactors = F)
head(dat)

#narrow to columns of interest
dat2 <- dplyr::select(dat, land_use_type, village_code, CS,EN,crown_height)
head(dat2)

# first, replace "not endemic" with introduced
dat2$EN[dat2$EN=="not endemic"] <- "introduced"
#do an ordinal logistic regression
#response variable has a natural order
#if there is no natural order to the response variable, you use multinomial logistic regression.

#see here for multinomial regression:
# https://www.analyticsvidhya.com/blog/2016/02/multinomial-ordinal-logistic-regression/

#first, make your response variable into an ordered categorical variable
unique(dat2$CS)

#subset the data to those with known threat levels

dat3 <- subset(dat2, CS!="unknown" & CS!="uk" & CS!="" & CS!="DD")
head(dat3)
unique(dat3$CS)

#now order these variables logically
dat3$CS = factor(dat3$CS, levels=c("LC", "NT", "VU", "EN", "CR"))
unique(dat3$CS)


#make sure that the other variables are factors too
dat3$land_use_type = as.factor(dat3$land_use_type)
dat3$village_code = as.factor(dat3$village_code)
dat3$crown_height = as.factor(dat3$crown_height)

#reorder EN so that you use something other than endemic as a reference
dat3$EN = factor(dat3$EN, levels= c("native but not endemic", "endemic", "introduced", "naturalized"))

#write the response variable as a numeric rank
dat3$CS_rank <- as.factor(as.numeric(dat3$CS))
#check that it worked
dat3[dat3$CS=="EN",]

#the function polr is in the package MASS and alkows you to model these interactions with no random effects
#you can see the example of this here: https://stats.idre.ucla.edu/r/dae/ordinal-logistic-regression/
m1 <- polr(CS_rank ~ EN + crown_height + land_use_type, data = dat3, Hess = TRUE)
summary(m1)

#to add random effects, you'll need a new package
#try the clmm or clmm2 function in the package 'ordinal'
install.packages('ordinal')
library(ordinal)


# this does not converge because I haev eliminated plants with unknown or data deficitient threats
# and that means that some of the remaining villages do not have all land use types or 

#this function fails to converge
m2 <- clmm(CS_rank ~  crown_height + EN + land_use_type + (1|village_code), data = dat3, Hess = T)
#this uses a new algorithm and converges burt not for EN == naturalized because the 
#response variable is ALWAYS the same (LC)
m2 <- clmm2(CS_rank ~  crown_height + EN + land_use_type, random=village_code, data = dat3, Hess = 1)
summary(m2)

#so you either accept the partial convergence or you drop this term
m3 <- clmm(CS_rank ~ crown_height + land_use_type + (1|village_code), data = dat3, Hess = T)
summary(m3) # this converges when you don't use EN as a predictor

#up to you whether to use m1, m2, or m3
#fine to use m2 but you need to explain why it maybe didn't converge (because all naturalized plants have a repsonse of LC)

