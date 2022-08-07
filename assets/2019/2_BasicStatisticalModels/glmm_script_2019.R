# In this tutorial we will learn to use statistical methods to study dynamical systems
# The database we will use are utilization rates from 4 health centers (CSB) during three years.
# 2 CSBs (A,B) have benefited from a user-fee removal intervention since Feb2014, and 2 CSBs (C,D) have
# benefited from 2 user-fee removal interventions, since Feb2014 and Oct2014. We will study the impact
# of each intervention on utilization rates. 
# We include some covariates that could affect utilization rates and bias our estimations: 
# number of staff per month (staff) and number of referals per month (ref), time-dependent variables

# Load data # 
# Write the path of the folder where your database is, and specify database name and extension
setwd('C:/Users/Garchito/Documents/Teaching and courses/E2M2 2018/Statistics')
csb.data=read.csv('csb database.csv')
csb.data$Date=as.Date(csb.data$Date)

# Explore the database #
head(csb.data)
dim(csb.data)
names(csb.data)
summary(csb.data)
str(csb.data)

#------------------------------------------------------
# 1. Data exploration: Plots and summary stadistics
#------------------------------------------------------

# Distribution of outcome variable
hist(csb.data$outpatient, col='grey', main='',xlab='Number of outpatient visits per month')

# Outpatient visits for each health center (exploration of association for categorical variables)
boxplot(csb.data$outpatient~csb.data$csb, ylab='Number of outpatient visits per month')

# Outpatient visits after interventions were in place
boxplot(csb.data$outpatient~csb.data$int1, xlab="Intervention 1 in place", ylab='Number of outpatient visits per month')
boxplot(csb.data$outpatient~csb.data$int2, xlab="Intervention 2 in place", ylab='Number of outpatient visits per month')

# Correlation plots (exploration of association for quantitative variables)
plot(csb.data$staff,csb.data$outpatient, ylab='Number of outpatient visits per month')
abline(lm(outpatient~staff, data=csb.data))
plot(csb.data$ref,csb.data$outpatient, ylab='Number of outpatient visits per month')
abline(lm(outpatient~ref, data=csb.data))

#-------------------------------------------------
# 2. Exploratory analyses and graphics over time
#-------------------------------------------------

# Plot trends over time 
plot(csb.data$Date, csb.data$outpatient)
abline(lm(outpatient~Date, data=csb.data))
ts.avg=tapply(csb.data$outpatient, csb.data$Date, mean, na.rm=T)
lines(unique(csb.data$Date), ts.avg, type='l', col='red', lwd=2)

# Create a seasonal variable to better assess trends in our model (otherwise our estimates will be biased)
  # Estimate average number of consultations per month and plot that 
  monthly.avg=tapply(csb.data$outpatient, csb.data$month, mean, na.rm=T)
  plot(1:12, monthly.avg, type='b')
  # It's hard to see the seasonal trend ploting just one year 
  # This would be the equivalent over three years:
  x=rep(1:36) ; y=rep(monthly.avg,3)
  plot(x, y,type='b')
  
  # It seems that there's a clear seasonal effect. We'll construct several seasonal effects 
  # and see which one seems more appropriate. A seasonal effect looks as follows: 
  shift=1
  period=12
  season=sin(2*pi*(csb.data$month-shift)/period)
  # where "shift" is the horizontal shift of the sine function in a X axis and "period"
  # is the number of months it takes for the sine function to be repeated. For instance:
  season1=sin(2*pi*(x-0)/36)
  season2=sin(2*pi*(x-0)/24)
  season3=sin(2*pi*(x-0)/12)
  season4=sin(2*pi*(x-1)/12)
  season5=sin(2*pi*(x-2)/12)
  season6=sin(2*pi*(x-3)/12)
  
  # Plot of seasonal effects with different periods
  par(mfrow=c(2,1))
  plot(x, y,type='b', main="Observations")
  plot(x,season1, type='l', col='red', xlab='Month', ylab='Value',lwd=2, main="Seasonal effect")
  lines(x,season2,col='blue', lwd=2)
  lines(x,season3,col='blue', lwd=2, lty=3)
  
  # Plot of seasonal effects with different horizontal shifts
  par(mfrow=c(2,1))
  plot(x, y,type='b', main="Observations")
  plot(x,season3, type='l', col='red', xlab='Month', ylab='Value',lwd=2, main="Seasonal effect")
  lines(x,season4,col='blue', lwd=2)
  lines(x,season5,col='blue', lwd=2, lty=2)
  lines(x,season6,col='blue', lwd=2, lty=3)
  
  # Based on visual inspection, we select the one that best fits the data:
  par(mfrow=c(2,1))
  plot(x, y,type='b', main="Observations")
  plot(x,sin(2*pi*(x-11)/12),col='blue', lwd=2, type='b', main="Seasonal effect")
  csb.data$season=sin(2*pi*(unique(csb.data$month-11))/12)

#---------------------------
# 4. Fit a linear model 
#---------------------------
  
# Before fitting a linear model we need to think about our outcome variable in order to assess
# which model family to use. Our utilization rates have discrete, non-zero values and look non-normal
par(mfrow=c(1,1))
hist(csb.data$outpatient, col='grey', main='',xlab='Number of outpatient visits per month')

# We thus could think of a poisson model, but let's see if it fulfills the assumptions:
mean(csb.data$outpatient, na.rm=T)
var(csb.data$outpatient, na.rm=T)

# When the variance is higher than the mean, we have overdispersion. A model that deals fine with this
# is called "Negative Binomial".

# The second thing we need to think about is whether all our observations are independent. In the case 
# where we have repeated observations of the same individuals/sites/healthcenters, observations are not independent.
# The values in a site of a certain month (x) are correlated to the values of that site the month before (x-1)
# To deal with this we will use "mixed-effects models". The best package for doing this is called lme4:
library(lme4)

# A mixed-effects model includes 2 components: a "fixed effect", which is common for all individuals/sites (similar to a linear model)
# and a random effect, which is specific of the individual/site. This random effect can be a random intercept or a random slope.

# Let's first construct a univariate model (one variable) to see the basic structure and outcomes of the model
m1=glmer.nb(outpatient~int1+ (1 | csb), data=csb.data)
summary(m1)
# We can look at the random intercept (see how the intercept changes for each CSB)
ranef(m1)

# We'll check a multivariate model that includes number of medicall staff and all time-dependent variables 
# (we don't consider a full model for simplicity and lack of sufficient observations). 
m1=glmer.nb(outpatient~ staff+ int1+ int2+ season+ (1 | csb), data=csb.data )
summary(m1)

#--------------------------------------
# 5. Back-transformation of estimates 
#--------------------------------------
# Remember that each family has a certain link to associate the linear predictor with the response variable
# To be able to interpret the results given in the model output, we need to back-transform them
# Let's look at one coefficient and it's confidence intervals
m1tab=summary(m1)$coefficients ; m1tab

intervention1=exp(m1tab[3,1]) ; intervention1
lowerCI=exp(m1tab[3,1]-1.96*m1tab[3,2]) ; lowerCI
upperCI=exp(m1tab[3,1]+1.96*m1tab[3,2]) ; upperCI

intervention2=exp(m1tab[4,1]) ; intervention2
lowerCI=exp(m1tab[4,1]-1.96*m1tab[4,2]) ; lowerCI
upperCI=exp(m1tab[4,1]+1.96*m1tab[4,2]) ; upperCI

#-----------------------
# 6. Check model fit 
#-----------------------
# Let's see how our model predicts utilization patterns
# Overall trends
predictions=exp(predict(m1,csb.data))
ts.pred=tapply(predictions, csb.data$Date, mean, na.rm=T)
plot(csb.data$Date, csb.data$outpatient)
lines(unique(csb.data$Date), ts.avg, col='red', lwd=2)
lines(unique(csb.data$Date), ts.pred, col='blue', lwd=2, lty='dashed')

# For each health center
time.temp=csb.data$Date[!is.na(csb.data$outpatient)]
outcome.var.temp=csb.data$outpatient[!is.na(csb.data$outpatient)]
par(mfrow=c(2,2))
for (i in 1:4){
  x=time.temp[m1@ frame$csb==unique(m1@ frame$csb)[i]]
  y1=outcome.var.temp[m1@ frame$csb==unique(m1@ frame$csb)[i]]
  y2=fitted(m1)[m1@ frame$csb==unique(m1@ frame$csb)[i]]
  
  plot(as.numeric(x),y1,xlab="Time",ylab='Utilization',cex.lab=1.3,cex.axis=1.3,main=unique(m1@ frame$csb)[i], pch=16, col='blue')
  lines(as.numeric(x),y2,lty=2, col='red',lwd=2)
}

#----------------------------
# 6. Check model assumptions 
#----------------------------

# a. Normality of the residuals
par(mfrow=c(1,2))
qqnorm(resid(m1))
qqline(resid(m1))
hist(resid(m1))
shapiro.test(resid(m1))

# b. Predicted vs Residuals (homogeneity)
par(mfrow=c(1,1))
plot(fitted(m1),resid(m1))
abline(lm(resid(m1)~fitted(m1)), lwd=2, col="red")

# c. Residuals vs Variables (independence)
par(mfrow=c(2,2))
varList=names(m1@ frame)[2:(length(names(m1@frame))-1)]
for (i in 1:length(varList)){
  var=m1@ frame[,1+i]
  plot(var,resid(m1),main=varList[i])
  abline(lm(resid(m1)~var), lwd=2, col="red")
}

# d. Temporal correlation of residuals (independence)
# The autocorrelation function (ACF) evaluates the correlation between the values of time x with the values 
# of time x-1, x-2, x-3, etc. This is called lag, and the residuals shouldn't show correlations at any lags

# We need to prepare a little the database in order to add NAs into the residual structure
E= residuals(m1, type="pearson")
I1=!is.na(csb.data$outpatient)
Efull=NA
Efull[I1]=E   # This is done to add NAs in the right place
temp=data.frame(csb.data$csb,Efull) ; colnames(temp)[1]="csb"
# Plot 
par(mfrow=c(2,2))
for (i in 1:4){
  acf(temp[temp$csb==levels(temp$csb)[i],2],main=levels(temp$csb)[i], na.action=na.pass)
}

# After validation of fit and assumptions, we can see that we haven't done yet a great job at predicting
# utilization patterns. The fit is not super accurate, especially for some CSBs, and we can see temporal 
# autocorrelation of residuals among other issues. A way forward would be to include lagged utilization 
# into the model, test for progressive changes in addition to abrupt ones, include other covariates that 
# we may think could bias our estimates, or explore a model with random slopes.  