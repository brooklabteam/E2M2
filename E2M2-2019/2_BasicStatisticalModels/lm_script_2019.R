# In this tutorial we will learn to use basic statistical methods, from exploratory analyses to linear models
# The database we will use are a group  of 100 lemurs for which a number of measurements have been carried out in the field 
# These include height (in cm), age (in months), gender, number of gastrointestinal parasites found in their faeces, and malaria infection
# We will explore the relationship between these different variables, and will develop statistical models to quantify their association

# Load data # 
# Write the path of the folder where your database is, and specify database name and extension
lemur.data=read.table('C:/Users/garchito/Documents/Teaching and courses/E2M2 2019/Statistics/Lemur data.txt', head=T)

#------------------------------------------------------
# 1. Data exploration: Plots and summary stadistics
#------------------------------------------------------
# Explore the database #
head(lemur.data)
dim(lemur.data)
names(lemur.data)
summary(lemur.data)
str(lemur.data)
attach(lemur.data)

# Distribution of outcome variable
hist(taille, col='grey')
boxplot(taille, ylab='Height (cm)')

# Distribution of explanatory variables
par(mfrow=c(2,2))
hist(age, col='grey') ; plot(sexe, col='grey', main="Gender") ; hist(GIparasites, col='grey') ; plot(malaria, col='grey', main="Malaria") ;
pairs(lemur.data)

#-------------------------------------------------------------
# 2. Explore individual associations: Plots and linear models
#-------------------------------------------------------------
par(mfrow=c(1,1))
# Relationship between height and age (univariate analyses)
plot(age, taille, pch=19)
m1=lm(taille~age, data=lemur.data)
abline(m1)
summary(m1)

# Relationship between height and GI parasites (univariate analyses)
plot(GIparasites, taille, pch=19)
m2=lm(taille~GIparasites, data=lemur.data)
abline(m2)
summary(m2)

# Relationship between height and gender (univariate and bivariate analyses)
boxplot(taille~sexe)
plot(age, taille, col=c('blue','red')[sexe],pch=19)
abline(lm(taille[sexe=='Female']~age[sexe=='Female']), col='blue',lwd=2)
abline(lm(taille[sexe=='Male']~age[sexe=='Male']), col='red',lwd=2)

m3=lm(taille~age+sexe, data=lemur.data)
summary(m3)

#---------------------------------------
# 3. Create a multivariate linear model
#---------------------------------------
# Include all relevant variables
m4=lm(taille~age+sexe+GIparasites+malaria, data=lemur.data)
summary(m4)


# Stepwise selection to exclude variables that do not help explain the outcome variable
step(m4)   # One can also do forward selection (add one variable at a time and see how the model improves) ; do backwards selection (the opposite) ; or multimodel selection (a combination of methods)
m5=step(m4)
summary(m5)

# Simple verification of model assumptions
plot(m5)
shapiro.test(resid(m5))

#--------------------------------------
# 4. Create a generalized linear model 
#--------------------------------------
# Instead of height, let's explore determinants of parasite numbers and of infectious status
# For this, we need to understand first what type of outcome variable we're dealing with and the most appropriate model
hist(GIparasites, col='grey')
plot(sickGIparasites, col="grey", main="Infection status")

# The variable "Number of parasites" is count data and it's poisson distributed. 
# This type of variable is typically modelled with poisson models 
m6=glm(GIparasites~age+sexe++malaria, family='poisson', data=lemur.data)
summary(m6)

m7=step(m6)
summary(m7)

# The variable "infectious status" is a binary variable 
# This type of variable is typically modelled with binomial models (also known as logistic regression)
m8=glm(sickGIparasites~age+sexe++malaria, family='binomial', data=lemur.data)
summary(m8)

m9=step(m8)
summary(m9)

#----------------------------------------------
# 5. Back-transformation of estimates for GLMs
#----------------------------------------------
# Remember that each family has a certain link to associate the linear predictor with the response variable
# To be able to interpret the results given in the model output, we need to back-transform them
# Let's look at the intercept of our poisson model of GI parasites and one coefficient 
m7tab=summary(m7)$coefficients

intercept=exp(m7tab[1,1]) ; intercept
effet.sexe=exp(m7tab[2,1]) ; effet.sexe
