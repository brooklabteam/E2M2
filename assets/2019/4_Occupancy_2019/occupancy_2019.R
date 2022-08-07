#############################################################################
########--------------Introduction to occupancy modeling-----------##########
#############################################################################



####load library unmarked
library(unmarked)
###import data frame called fakedata.csv
data<-read.csv("fakedata.csv")
####data narrative and details
#researchers are interested in studying the prevalence
#of a parasite in a population. However, this pathogen was recently described in this species
#and no reliable test exist (no gold standard). 

#data: the researcher looked at blood smear slides three times and noted whether an animal was positive or negative for a parasite
#Pos, Pos2, Pos 3
#he/she also measured weight, age and sex and scaled the values so that mean=0, (x-mean/sd)

data$scaleW<-scale(data$Weight)
data$scaleA<-scale(data$Age)
data$scaleS<-data$Sex
data$scaleF<-scale(data$Forearm)
summary(data)


####columns 7 to 9 in the dataframe represent the different surveys 
####and we can consider each individual as a site. columns 13 to 16 are the site specific covariates 
y<-data[,7:9]
para.site<-data[,13:16]

######Unmarked uses a specific type of data frame (S4)
para <- unmarkedFrameOccu(y = y, siteCovs = para.site)
summary(para)

##################################################################################################
#CREATING MODELS
#simple occupancy model no covariates
fm1<-occu(~1 ~1,para)
fm1# Estimate for occupancy and detection are on logit scale and need transformation

#back transformations using the function backTransform 
#This can only work when there are no covariates
backTransform(fm1,'det')
backTransform(fm1,"state")
#probability of detection given the parasite is present, is only 44%
#probability that toxo is present is 12.5%

#Create some more occupancy models
#constant detection, occupancy predicted by Sex 
fm2<-occu(~1 ~scaleS,para)
#occupancy varies between Weight
fm3<-occu(~1~scaleW,para)
#occupancy predicted by Age
fm4<-occu(~1 ~scaleA,para)
#occupancy predicted by Forearm
fm5<-occu(~1~scaleF,para)

########Which of these models is the best??????
#built-in model selection function modSel
#labeled list of models
fmlist<-fitList("null"=fm1,"psi_Sex"=fm2,"psi_Weight"=fm3, "psi_Age"=fm4, "psi_Forearm"=fm5)
#model selection -- ranked list of models and AIC
modSel(fmlist)

####### Goodness of fit test (Chi square) to assess whether our model fits with the data
####ie. is our model valid?
###

chisq<-function(fm) {
  observed<-getY(fm@data)
  expected<-fitted(fm)
  sum((observed-expected)^2/expected)
}
(pb<-parboot(fm4,statistic=chisq,nsim=200, parallel=FALSE))

######predict to backtransform estimates when there are covariates in the model
####e.g. want to transform the occupancy estimates obtained from the model fm4
backTransform(fm4,'det')
#again, backTransform works only when there are no covariates to explain the detection or state process
#when there are covariates, need to build a 
lc<-linearComb(fm4,c(1,0),type="state") #predict the value of psi when scaledAge is at 0
# = when age is at its mean
backTransform(lc)
#we can also create a new data frame that spans the range of "Age-scaled",
#which will allo us to predict specific values of psi (Probability of occupancy) for a value of age scale
# Create a variable spanning the range of scaleAge, -2 to 2
plotData <- data.frame(scaleA = seq(from=-2, to=2, length=1200))
head(plotData)
# Compute predicted values of occupancy probability for age using our best model
AgePred <- predict(fm4, newdata=plotData, type="state",appendData=T)

# Plot probability of occupancy versus age including upper and lower bounds of 95% CI
plot(plotData$scaleA, AgePred[,"Predicted"], type="l",
     xlab="Age (scaled)", ylab="Probability of occupancy (95% CI)",ylim=c(0,1))
lines(plotData$scaleA, AgePred[,"Predicted"]-1.96*AgePred[,"SE"], lty=2)
lines(plotData$scaleA, AgePred[,"Predicted"]+1.96*AgePred[,"SE"], lty=2)
#another way to do that graph using ggplot2
library(ggplot2)
qplot(scaleA, Predicted,data=AgePred, geom = "line",
      xlab = "Scaled age", ylab = "Probability of occupancy", ylim=c(0,1))+
  geom_ribbon(aes(x = scaleA, ymin = lower, ymax = upper),alpha=0.1)+
  theme_classic()


######################################################################################
#####################
#occupancy models also allows us to estimate abundance, here the number of parasites in animals
##########
######################################################################################
#first we need identify in our dataframe the columns that specify
counts<-data[,10:12]

#load site level (individual) covariates 
para_count.site<-data[,13:16]

#######Unmarked uses a specific type of data frame (S4)
#but this time we use unmarkedFramePcount
para_count <- unmarkedFramePCount(y = counts, siteCovs = para_count.site)
summary(para_count)

#we build occupancy numbers
pc1<-pcount(~1~1,para_count,mixture="NB")
pc1
pc2<-pcount(~1~scaleA,para_count,mixture="NB")
pc2
coef(pc2)
confint(pc2, type = "state", level = 0.95)

#convert estimates from link-scale to original scale 
#"backTransform" works when no covariates
backTransform(pc1,"state")
backTransform(pc1,"det")

#when there are covariates we need to use linearcomb
lc<-linearComb(pc2,c(1,0),type="state")
backTransform(lc)
#we can also create a new data frame that spans the range of "Age-scaled",
#which will allo us to predict specific values of psi (Probability of occupancy) for a value of age scale
# Create a variable spanning the range of scaleAge, -2 to 2
plotData <- data.frame(scaleA = seq(from=-2, to=2, length=1200))
head(plotData)
# Compute predicted values of occupancy probability for age using our best model
AgePredAbund <- predict(pc2, newdata=plotData, type="state",appendData=T)

# Plot probability of occupancy versus age including upper and lower bounds of 95% CI
plot(plotData$scaleA, AgePredAbund[,"Predicted"], type="l",
     xlab="Age (scaled)", ylab="Probability of occupancy (95% CI)",ylim=c(0,100))
lines(plotData$scaleA, AgePredAbund[,"Predicted"]-1.96*AgePredAbund[,"SE"], lty=2)
lines(plotData$scaleA, AgePredAbund[,"Predicted"]+1.96*AgePredAbund[,"SE"], lty=2)


#another way of doing the same using ggplot package
qplot(plotData$scaleA, ylim=c(0,100),Predicted,data=AgePredAbund, geom = "line",
      xlab = "Scaled age", ylab = "Estimated abundance of parasites")+
  geom_ribbon(aes(x = plotData$scaleA, ymin = lower, ymax = upper),alpha=0.1)+
  theme_classic()
