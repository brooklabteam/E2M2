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




