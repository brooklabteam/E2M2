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


#then replace GDP growth with GDP
dat.long <- subset(dat.long, indicator!="gdp_growth")

#and attach the GDP data below

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

ggsave(file = "Fig1_Mirana_cumulative.pdf",
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
#and look at autocorrelation in the residuals:
acf(resid(gam1), lag.max = 36, main = "ACF") #no clear pattern so not an issue
pacf(resid(gam1), lag.max = 36, main = "pACF")

gam2 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="ODA"))
acf(resid(gam2), lag.max = 36, main = "ACF")
pacf(resid(gam2), lag.max = 36, main = "pACF")

gam3 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="export"))
acf(resid(gam3), lag.max = 36, main = "ACF")
pacf(resid(gam3), lag.max = 36, main = "pACF")

gam4 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="import"))
acf(resid(gam4), lag.max = 36, main = "ACF")
pacf(resid(gam4), lag.max = 36, main = "pACF")

gam5 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="debt"))
acf(resid(gam5), lag.max = 36, main = "ACF")
pacf(resid(gam5), lag.max = 36, main = "pACF")

gam6 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="broad_money"))
acf(resid(gam6), lag.max = 36, main = "ACF")
pacf(resid(gam6), lag.max = 36, main = "pACF")

gam7 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="nat_rent"))
acf(resid(gam7), lag.max = 36, main = "ACF")
pacf(resid(gam7), lag.max = 36, main = "pACF")

gam8 <- gam(quantity~s(yr), data = subset(dat.long, indicator=="inflation"))
acf(resid(gam8), lag.max = 36, main = "ACF")
pacf(resid(gam8), lag.max = 36, main = "pACF")

#looks like autocorrelation is not a huge issue here... we'll keep the GAMS as is

#here are the gam summaries you might want to report
summary(gam1)
summary(gam2)
summary(gam3)
summary(gam4)
summary(gam5)
summary(gam6)
summary(gam7)
summary(gam8)


#and add both the predictions (with upper and lower confidence limits) to the original data.
#do the same for the residuals
dat.long$gam_predict <- dat.long$gam_predict_lci <- dat.long$gam_predict_uci <- dat.long$gam_residuals <- dat.long$gam_residuals_lci <- dat.long$gam_residuals_uci <- NA
dat.long$gam_predict[dat.long$indicator=="GDP"] <- predict.gam(gam1, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="GDP"] <- predict.gam(gam1, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam1, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="GDP"] <- predict.gam(gam1, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam1, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="GDP"] <- residuals(gam1)
dat.long$gam_residuals_lci[dat.long$indicator=="GDP"] <- dat.long$gam_predict[dat.long$indicator=="GDP"] - dat.long$gam_predict_lci[dat.long$indicator=="GDP"]
dat.long$gam_residuals_uci[dat.long$indicator=="GDP"] <- dat.long$gam_predict[dat.long$indicator=="GDP"] - dat.long$gam_predict_uci[dat.long$indicator=="GDP"]


dat.long$gam_predict[dat.long$indicator=="ODA"] <- predict.gam(gam2, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="ODA"] <- predict.gam(gam2, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam2, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="ODA"] <- predict.gam(gam2, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam2, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="ODA"] <- residuals(gam2)
dat.long$gam_residuals_lci[dat.long$indicator=="ODA"] <- dat.long$gam_predict[dat.long$indicator=="ODA"] - dat.long$gam_predict_lci[dat.long$indicator=="ODA"]
dat.long$gam_residuals_uci[dat.long$indicator=="ODA"] <- dat.long$gam_predict[dat.long$indicator=="ODA"] - dat.long$gam_predict_uci[dat.long$indicator=="ODA"]



dat.long$gam_predict[dat.long$indicator=="export"] <- predict.gam(gam3, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="export"] <- predict.gam(gam3, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam3, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="export"] <- predict.gam(gam3, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam3, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="export"] <- residuals(gam3)
dat.long$gam_residuals_lci[dat.long$indicator=="export"] <- dat.long$gam_predict[dat.long$indicator=="export"] - dat.long$gam_predict_lci[dat.long$indicator=="export"]
dat.long$gam_residuals_uci[dat.long$indicator=="export"] <- dat.long$gam_predict[dat.long$indicator=="export"] - dat.long$gam_predict_uci[dat.long$indicator=="export"]


dat.long$gam_predict[dat.long$indicator=="import"] <- predict.gam(gam4, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="import"] <- predict.gam(gam4, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam4, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="import"] <- predict.gam(gam4, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam4, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="import"] <- residuals(gam4)
dat.long$gam_residuals_lci[dat.long$indicator=="import"] <- dat.long$gam_predict[dat.long$indicator=="import"] - dat.long$gam_predict_lci[dat.long$indicator=="import"]
dat.long$gam_residuals_uci[dat.long$indicator=="import"] <- dat.long$gam_predict[dat.long$indicator=="import"] - dat.long$gam_predict_uci[dat.long$indicator=="import"]


dat.long$gam_predict[dat.long$indicator=="debt"] <- predict.gam(gam5, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="debt"] <- predict.gam(gam5, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam5, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="debt"] <- predict.gam(gam5, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam5, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="debt"] <- residuals(gam5)
dat.long$gam_residuals_lci[dat.long$indicator=="debt"] <- dat.long$gam_predict[dat.long$indicator=="debt"] - dat.long$gam_predict_lci[dat.long$indicator=="debt"]
dat.long$gam_residuals_uci[dat.long$indicator=="debt"] <- dat.long$gam_predict[dat.long$indicator=="debt"] - dat.long$gam_predict_uci[dat.long$indicator=="debt"]


dat.long$gam_predict[dat.long$indicator=="broad_money"] <- predict.gam(gam6, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="broad_money"] <- predict.gam(gam6, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam6, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="broad_money"] <- predict.gam(gam6, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam6, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="broad_money"] <- residuals(gam6)
dat.long$gam_residuals_lci[dat.long$indicator=="broad_money"] <- dat.long$gam_predict[dat.long$indicator=="broad_money"] - dat.long$gam_predict_lci[dat.long$indicator=="broad_money"]
dat.long$gam_residuals_uci[dat.long$indicator=="broad_money"] <- dat.long$gam_predict[dat.long$indicator=="broad_money"] - dat.long$gam_predict_uci[dat.long$indicator=="broad_money"]


dat.long$gam_predict[dat.long$indicator=="nat_rent"] <- predict.gam(gam7, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="nat_rent"] <- predict.gam(gam7, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam7, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="nat_rent"] <- predict.gam(gam7, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam7, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="nat_rent"] <- residuals(gam7)
dat.long$gam_residuals_lci[dat.long$indicator=="nat_rent"] <- dat.long$gam_predict[dat.long$indicator=="nat_rent"] - dat.long$gam_predict_lci[dat.long$indicator=="nat_rent"]
dat.long$gam_residuals_uci[dat.long$indicator=="nat_rent"] <- dat.long$gam_predict[dat.long$indicator=="nat_rent"] - dat.long$gam_predict_uci[dat.long$indicator=="nat_rent"]


dat.long$gam_predict[dat.long$indicator=="inflation"] <- predict.gam(gam8, type = "response", se.fit=T)$fit
dat.long$gam_predict_lci[dat.long$indicator=="inflation"] <- predict.gam(gam8, type = "response", se.fit=T)$fit - 1.96*predict.gam(gam8, type = "response", se.fit=T)$se.fit
dat.long$gam_predict_uci[dat.long$indicator=="inflation"] <- predict.gam(gam8, type = "response", se.fit=T)$fit + 1.96*predict.gam(gam8, type = "response", se.fit=T)$se.fit
dat.long$gam_residuals[dat.long$indicator=="inflation"] <- residuals(gam8)
dat.long$gam_residuals_lci[dat.long$indicator=="inflation"] <- dat.long$gam_predict[dat.long$indicator=="inflation"] - dat.long$gam_predict_lci[dat.long$indicator=="inflation"]
dat.long$gam_residuals_uci[dat.long$indicator=="inflation"] <- dat.long$gam_predict[dat.long$indicator=="inflation"] - dat.long$gam_predict_uci[dat.long$indicator=="inflation"]

dat.long$strip_label = dat.long$indicator
dat.long$strip_label[dat.long$strip_label=="broad_money"] = "broad money"
dat.long$strip_label[dat.long$strip_label=="nat_rent"] = "nat rent"

crisis.dat2 = crisis.dat
crisis.dat2$strip_label <- "broad money"

ggplot(data=dat.long) + geom_line(aes(x=yr, y=quantity)) + facet_grid(strip_label~., scales="free_y") + theme_bw() + 
  theme(panel.grid = element_blank(), strip.background = element_blank()) + ylab("") + xlab("") +
  geom_ribbon(aes(x=yr, ymin=gam_predict_lci, ymax=gam_predict_uci), alpha=.3, fill="blue") +
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + geom_line(aes(x=yr, y= gam_predict), col="blue", size = 1) +
  geom_label(data=crisis.dat2, aes(x=crisis_yr, y=25, label=crisis_names)) 

ggsave(file = "Fig2_Mirana_cumulative_wGAM.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)



#then plot residuals

ggplot(data=dat.long)  + geom_hline(yintercept = 0, col="blue")  + 
  geom_ribbon(aes(x=yr, ymin=gam_residuals_lci, ymax=gam_residuals_uci), fill = "blue", alpha=.3) +
  geom_vline(data=crisis.dat, aes(xintercept=crisis_yr), col="red", size = 1) + xlab("") +
  ylab("residuals from predictions") + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank()) +
  geom_point(aes(x=yr, y=gam_residuals), size=3) + geom_label(data=crisis.dat2, aes(x=crisis_yr, y=6, label=crisis_names)) +  
  facet_grid(strip_label~., scales="free_y")


ggsave(file = "Fig3_Mirana_cumulative_wGAM_residuals.pdf",
       units="mm",  
       width=60, 
       height=90, 
       scale=3, 
       dpi=300)


#now select those years that are outside the bounds of the confidence interval (blue band) for multiple economic indicators
#and assess what is common across the years with indicators that are consistent outliers.
