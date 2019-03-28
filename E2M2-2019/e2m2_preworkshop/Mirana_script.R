rm(list=ls())

library(reshape2)
library(ggplot2)
library(dplyr)
library(lme4)

dat = read.csv(file ="mirana-data.csv", header = T, stringsAsFactors = F)
head(dat)
names(dat)
names(dat) <- c("yr", "gdp_growth", "ODA", "export", "import", "debt", "broad_money", "nat_rent", "inflation")
head(dat)


#separate and rejoin using 'melt' function in reshape2

dat.long = melt(dat,
                id.vars = "yr",
                variable.name = "indicator")
                
head(dat.long)
names(dat.long)[3] <- "quantity"


#look at indicator names
unique(dat.long$indicator)

crisis_yr = c(1972, 1991, 2002, 2009)
crisis_names = c(1972, 1991, 2002, 2009)
crisis.dat = cbind.data.frame(crisis_yr, crisis_names)
crisis.dat2 = crisis.dat
crisis.dat2$indicator = "broad_money"
head(crisis.dat)

#as points
ggplot(dat = dat.long) + geom_point(aes(x=yr, y=quantity)) + facet_grid(indicator~., scales = "free_y") + 
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 3) + geom_label(data=crisis.dat2, aes(x=crisis_yr, y=25, label=crisis_names))

#as lines - much clearer
ggplot(dat = dat.long) + geom_line(aes(x=yr, y=quantity, color=indicator)) + facet_grid(indicator~., scales = "free_y") + 
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_label(data=crisis.dat2, aes(x=crisis_yr, y=25, label=crisis_names))

ggsave(file = "Mirana_cumulative.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)


#then replace GDP growth with GDP
dat.long <- subset(dat.long, indicator!="gdp_growth")

#and attach the GDP data belwo

#regression on export
head(dat.long)
dat.long$crisis_yr = "N"
dat.long$crisis_yr[dat.long$yr==1972] = "Y"
dat.long$crisis_yr[dat.long$yr==1991] = "Y"
dat.long$crisis_yr[dat.long$yr==2002] = "Y"
dat.long$crisis_yr[dat.long$yr==2009] = "Y"

dat.long$post_crisis_yr = "N"
dat.long$post_crisis_yr[dat.long$yr==1973] = "Y"
dat.long$post_crisis_yr[dat.long$yr==1992] = "Y"
dat.long$post_crisis_yr[dat.long$yr==2003] = "Y"
dat.long$post_crisis_yr[dat.long$yr==2010] = "Y"






#now try loading the raw GDP data
dat3 <- read.csv(file="GDP.csv", header=T, stringsAsFactors = F)

head(dat3)
names(dat3) = dat3[3,]
dat3 <- dat3[4:nrow(dat3),]

head(dat3)

dat3 <- dat3[dat3$`Country Name`=="Madagascar",]

gdp.dat <- melt(dat3, 
                 id.vars = c("Country Name", "Country Code", "Indicator Name", "Indicator Code"),
                 variable.name = "Year")

head(gdp.dat)

#and plot
gdp.dat$Year = as.numeric(as.character(gdp.dat$Year))
gdp.dat$value = as.numeric(as.character(gdp.dat$value))
ggplot(data=gdp.dat) + geom_line(aes(x=Year, y=value)) + ylab("GDP (US$)") +
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_label(data=crisis.dat, aes(x=crisis_yr, y=10, label=crisis_names))

head(gdp.dat)
gdp.dat <- dplyr::select(gdp.dat, -('Country Code'), -('Indicator Name'), -('Indicator Code'))

head(gdp.dat)


#and join
head(dat.long)

gdp.dat$crisis_yr <- "N"
gdp.dat$crisis_yr[gdp.dat$Year==1972 | gdp.dat$Year==1991 | gdp.dat$Year==2002 | gdp.dat$Year==2009] <- "Y"


gdp.dat$post_crisis_yr <- "N"
gdp.dat$post_crisis_yr[gdp.dat$Year==1973 | gdp.dat$Year==1992 | gdp.dat$Year==2003 | gdp.dat$Year==2010] <- "Y"

head(gdp.dat)
gdp.dat <- dplyr::select(gdp.dat, Year, value, crisis_yr, post_crisis_yr)
gdp.dat$indicator <- "GDP"
names(gdp.dat)[1:2] <-c("yr", "quantity")

gdp.dat <- dplyr::select(gdp.dat, names(dat.long))

dat.long <- rbind(gdp.dat, dat.long)


#as lines - much clearer
ggplot(dat = dat.long) + geom_line(aes(x=yr, y=quantity, color=indicator)) + facet_grid(indicator~., scales = "free_y") + 
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_label(data=crisis.dat2, aes(x=crisis_yr, y=25, label=crisis_names))

ggsave(file = "Mirana_cumulative.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)


## now that you have consistent fitted time series, go ahead and fit a GAM model to each of these time series to look at what years represent statistical anomalies
## So, for each variable try the following structure:

#try fitting a better model - gam
library(mgcv)

#fit one for every indicator
unique(dat.long$indicator)

#eliminate gaps in the data
dat.long <- dat.long[complete.cases(dat.long),]

gam1 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="GDP"))
gam2 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="ODA"))
gam3 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="export"))
gam4 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="import"))
gam5 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="debt"))
gam6 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="broad_money"))
gam7 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="nat_rent"))
gam8 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="inflation"))

dat.long$gam_predict <- dat.long$gam_residuals <- NA
dat.long$gam_predict[dat.long$indicator=="GDP"] <- predict.gam(gam1)
dat.long$gam_residuals[dat.long$indicator=="GDP"] <- residuals(gam1)

dat.long$gam_predict[dat.long$indicator=="ODA"] <- predict.gam(gam2)
dat.long$gam_residuals[dat.long$indicator=="ODA"] <- residuals(gam2)

dat.long$gam_predict[dat.long$indicator=="export"] <- predict.gam(gam3)
dat.long$gam_residuals[dat.long$indicator=="export"] <- residuals(gam3)

dat.long$gam_predict[dat.long$indicator=="import"] <- predict.gam(gam4)
dat.long$gam_residuals[dat.long$indicator=="import"] <- residuals(gam4)

dat.long$gam_predict[dat.long$indicator=="debt"] <- predict.gam(gam5)
dat.long$gam_residuals[dat.long$indicator=="debt"] <- residuals(gam5)

dat.long$gam_predict[dat.long$indicator=="broad_money"] <- predict.gam(gam6)
dat.long$gam_residuals[dat.long$indicator=="broad_money"] <- residuals(gam6)

dat.long$gam_predict[dat.long$indicator=="nat_rent"] <- predict.gam(gam7)
dat.long$gam_residuals[dat.long$indicator=="nat_rent"] <- residuals(gam7)

dat.long$gam_predict[dat.long$indicator=="inflation"] <- predict.gam(gam8)
dat.long$gam_residuals[dat.long$indicator=="inflation"] <- residuals(gam8)

crisis.dat2 = crisis.dat
crisis.dat2$indicator <- "broad_money"

ggplot(data=dat.long) + geom_line(aes(x=yr, y=quantity)) + facet_grid(indicator~., scales="free_y") +
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_line(aes(x=yr, y= gam_predict), col="blue", size = 1) +
  geom_label(data=crisis.dat2, aes(x=crisis_yr, y=25, label=crisis_names)) 

ggsave(file = "Mirana_cumulative_wGAM.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)



#then plot residuals

ggplot(data=dat.long)  + geom_hline(yintercept = 0)  + geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + ylab("residuals from predictions") +
   geom_point(aes(x=yr, y=gam_residuals), size=4) + geom_label(data=crisis.dat2, aes(x=crisis_yr, y=6, label=crisis_names)) +  facet_grid(indicator~., scales="free_y")





ggsave(file = "Mirana_cumulative_wGAM_residuals.pdf",
       units="mm",  
       width=60, 
       height=90, 
       scale=3, 
       dpi=300)


#then add in autocorrelation







datWEO <- read.csv(file="WEO-Data-FMI-Mada.csv", header=T, stringsAsFactors = F)
head(datWEO)
names(datWEO)
datWEO <- cbind.data.frame(datWEO[,2], datWEO[,5], datWEO[,7], datWEO[,10:48])  
names(datWEO) <- c("country", "variable", "units", seq(1980, 2018,1))

datWEO.long <- melt(datWEO, id.vars = c("country", "variable", "units"),
                    variable.name = "yr")
head(datWEO.long)

unique(datWEO.long$variable)
unique(datWEO.long$units)

#select only those with USD
datWEO.long = subset(datWEO.long, units == "U.S. dollars")
head(datWEO.long)
unique(datWEO.long$variable)

#plot GDP
WEO.gdp <- subset(datWEO.long, variable=="Gross domestic product, current prices")
head(WEO.gdp)
WEO.gdp$value <- as.numeric(WEO.gdp$value)
WEO.gdp$yr <- as.numeric(as.character(WEO.gdp$yr))

ggplot(data=WEO.gdp) + geom_line(aes(x=yr, y=value)) + ylab("GDP (US$)") +
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_label(data=crisis.dat, aes(x=crisis_yr, y=10, label=crisis_names)) 
  #geom_line(aes(x=Year, y= gam_predict), col="blue", size = 2)

ggsave(file = "WEO_gdp.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)


## looks the same, which is good.


## take a look at the crisis data that Mirana compiled and do some AIC model comparison
## 
crisis.dat <- read.csv2(file = "crisesA.csv", header=T, stringsAsFactors = F)
head(crisis.dat)
crisis.dat$Election <-  gsub(pattern = ",", replacement = ".", x=crisis.dat$Election)
crisis.dat$Election <-  gsub(pattern = "%", replacement = "", x=crisis.dat$Election)
crisis.dat$Election <- as.numeric(crisis.dat$Election)






################################
################################

## Old analyses - can mostly ignore
#add a column for change in GDP 1 year, 5 years, 10 years later

## compare predictors of crisis severity
## start with most complex crisis

m1 <- lm(Death~Transition + Pre.crisis + Post.crisis + Crisis, data=crisis.dat)
summary(m1)


m2 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="ODA"))
summary(m2)

m3 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="export"))
summary(m3)

m4 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="import"))
summary(m4)

m5 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="debt"))
summary(m5)

m6 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="broad_money"))
summary(m6)

m7 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="nat_rent"))
summary(m7)

m8 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="inflation"))
summary(m8)


#look at whether crisis matters each post-crisis yr
m1 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="gdp_growth"))
summary(m1)


m2 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="ODA"))
summary(m2)


m3 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="export"))
summary(m3)
##this sig - want to know the duration of effect...

m4 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="import"))
summary(m4)

m5 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="debt"))
summary(m5)

m6 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="broad_money"))
summary(m6)

m7 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="nat_rent"))
summary(m7)

m8 <- glm(quantity~yr +  post_crisis_yr, data=subset(dat.long, indicator=="inflation"))
summary(m8)




#look at whether crisis matters each year
m1 <- glm(quantity~yr +  crisis_yr, data=subset(dat.long, indicator=="gdp_growth"))
summary(m1)




gdp.dat = subset(dat.long, indicator=="gdp_growth")
gdp.dat = gdp.dat[complete.cases(gdp.dat),]
gdp.dat$prediction = predict(m1, type="response")
#this is sig

#try plotting with projection

#as lines - much clearer
ggplot(dat = gdp.dat) + geom_line(aes(x=yr, y=quantity), color="blue") + facet_grid(indicator~., scales = "free_y") + geom_line(aes(x=yr, y=prediction), col="black", size=2) +
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_label(data=crisis.dat, aes(x=crisis_yr, y=10, label=crisis_names))

#find the residuals
gdp.dat$residuals = resid(m1)

ggplot(data=gdp.dat)  + geom_hline(yintercept = 0)  + geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + ylab("residuals from predictions, GDP growth") +
  geom_label(data=crisis.dat, aes(x=crisis_yr, y=6, label=crisis_names)) + geom_point(aes(x=yr, y=residuals), size=4)

ggsave(file = "Mirana_GDP_growth_residuals.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)






m2 <- glm(value~Year + crisis_yr, data=gdp.dat)
summary(m2) 

m3 <- glm(value~Year, data=gdp.dat)
summary(m3) 

gdp.dat$lm_predict = predict(m3)
gdp.dat$lm_residuals = residuals(m3)

ggplot(data=gdp.dat) + geom_line(aes(x=Year, y=value)) + ylab("GDP (US$)") +
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_label(data=crisis.dat, aes(x=crisis_yr, y=10, label=crisis_names)) +
  geom_line(aes(x=Year, y= lm_predict), col="blue")


ggplot(data=gdp.dat)  + geom_hline(yintercept = 0)  + geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + ylab("residuals from predictions, GDP") +
  geom_point(aes(x=Year, y=lm_residuals), size=4) + geom_label(data=crisis.dat, aes(x=crisis_yr, y=6, label=crisis_names))

ggplot(data=gdp.dat)  + geom_hline(yintercept = 0)  + geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + ylab("residuals from predictions, GDP") +
  geom_point(aes(x=Year, y=gam_residuals), size=4) + geom_label(data=crisis.dat, aes(x=crisis_yr, y=6, label=crisis_names))




