rm(list=ls())

library(ggplot2)
library(plyr)
library(dplyr)
library(lubridate)
library(mgcv)


#read in and examin data
dat = read.csv(file ="tbdata_Mada.csv", header = T, stringsAsFactors = F)
head(dat)
unique(dat$tb_type)
unique(dat$tb_type_cat)

#make full date column for your datasets until you get the raw data
dat$day = 01
dat$date = paste( dat$year, dat$month, dat$day, sep="-")
dat$date = as.Date(dat$date, format = "%Y-%m-%d")
head(dat)


#here, load a dataset you created that gives the start and end date of the winter and summer season for all years in your dataset
#the ymx value gives the height of the bar plot as appropriate to the max number of observed cases for each type of TB
load("season.dat.Rdata")


#let's first look at whether tb is seasonally by type across the full lenth of our time series
#summarize cases by month and type across all years
dat.sum = ddply(dat, .(year, month, date, tb_type_cat), summarize, tot_case = length(ID) )
head(dat.sum)
#remove those rows with NA values:
dat.sum <- dat.sum[complete.cases(dat.sum),]

dat.sum.all.1 <- subset(dat.sum, tb_type_cat==1)
dat.sum.all.2 <- subset(dat.sum, tb_type_cat==2)
dat.sum.all.3 <- subset(dat.sum, tb_type_cat==3)
dat.sum.all.4 <- subset(dat.sum, tb_type_cat==4)

#do a separate GAM for each of the 4 types of TB. fix smoothing terms at K as recommended by package author
gam1 <- gam(tot_case ~ s(as.numeric(date), k=7,  bs="cr"), family="poisson", data=dat.sum.all.1)
gam2 <- gam(tot_case ~ s(as.numeric(date), k=7,  bs="cr"), family="poisson", data=dat.sum.all.2)
gam3 <- gam(tot_case ~ s(as.numeric(date), k=7,  bs="cr"), family="poisson", data=dat.sum.all.3)
gam4 <- gam(tot_case ~ s(as.numeric(date), k=7,  bs="cr"), family="poisson", data=dat.sum.all.4)

#and look at output of GAMs
summary(gam1) #sig s() value shows that there is significant seasonality for type 1 TB
summary(gam2) #sig s() value shows that there is significant seasonality for type 2 TB
summary(gam3) #sig s() value shows that there is significant seasonality for type 3 TB
summary(gam4) #p-value is lower (.01) so this type of TB is not really significantly seasonal

#now add the GAM predictions to the datasets so that you can plot it
dat.sum.all.1$gam_pred <- predict(gam1, type="response", se.fit=T)$fit
dat.sum.all.2$gam_pred <- predict(gam2, type="response", se.fit=T)$fit
dat.sum.all.3$gam_pred <- predict(gam3, type="response", se.fit=T)$fit
dat.sum.all.4$gam_pred <- predict(gam4, type="response", se.fit=T)$fit

#and inlcude lower 95% confidence intervals
dat.sum.all.1$gam_pred_lci <- predict(gam1, type="response", se.fit=T)$fit - 1.96*predict(gam1, type="response", se.fit=T)$se.fit
dat.sum.all.2$gam_pred_lci <- predict(gam2, type="response", se.fit=T)$fit - 1.96*predict(gam2, type="response", se.fit=T)$se.fit
dat.sum.all.3$gam_pred_lci <- predict(gam3, type="response", se.fit=T)$fit - 1.96*predict(gam3, type="response", se.fit=T)$se.fit
dat.sum.all.4$gam_pred_lci <- predict(gam4, type="response", se.fit=T)$fit - 1.96*predict(gam4, type="response", se.fit=T)$se.fit

#anbd upper 95% confidence intervals
dat.sum.all.1$gam_pred_uci <- predict(gam1, type="response", se.fit=T)$fit + 1.96*predict(gam1, type="response", se.fit=T)$se.fit
dat.sum.all.2$gam_pred_uci <- predict(gam2, type="response", se.fit=T)$fit + 1.96*predict(gam2, type="response", se.fit=T)$se.fit
dat.sum.all.3$gam_pred_uci <- predict(gam3, type="response", se.fit=T)$fit + 1.96*predict(gam3, type="response", se.fit=T)$se.fit
dat.sum.all.4$gam_pred_uci <- predict(gam4, type="response", se.fit=T)$fit + 1.96*predict(gam4, type="response", se.fit=T)$se.fit

#rejoin your data across all the types of TB:
dat.sum.all <- rbind(dat.sum.all.1, dat.sum.all.2, dat.sum.all.3, dat.sum.all.4)

#give a new label by which you will factor your data
dat.sum.all$type_label <- paste0("TB type-", dat.sum.all$tb_type_cat)

#then plot the results of this model with the data across all 4 types of TB:
#first specify the colors you wnat for your winter/summer bars:
filz1 = c( "summer" = "yellow", "winter" = "seagreen")

pTB_all_long <- ggplot(data=dat.sum.all) + geom_point(aes(x=date, y=tot_case)) + facet_grid(type_label~., scales="free_y") + xlab("") +
   scale_fill_manual(values=filz1) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  geom_line(aes(x=date, y=gam_pred)) + geom_ribbon(aes(x=date, ymin=gam_pred_lci, ymax=gam_pred_uci), alpha=.3, show.legend = F) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=ymx, fill=season), alpha=.3) #+ coord_cartesian(ylim = c(0,55))
print(pTB_all_long)

ggsave(file = "Fig1_all_TB_type_longitudinal.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)


#Now, to get a better sense of seasonality, try looking at the total counts by month, where your GAM is constrained within a calendar year
#you can use a cyclic cubic spline ("cc") that means that the output at the end of December is the same as the input at the beginning of January.
#Also can include random effects of year.

#First, summarize TB type counts by month:
month.new.sum = ddply(dat, .(tb_type_cat, year, month), summarize, monthly_cases=length(ID))
head(month.new.sum)
month.new.sum <- arrange(month.new.sum, month, tb_type_cat)

#get rid of NAs:
month.new.sum <- month.new.sum[complete.cases(month.new.sum),]
month.sum.1 <- subset(month.new.sum, tb_type_cat==1)
month.sum.2 <- subset(month.new.sum, tb_type_cat==2)
month.sum.3 <- subset(month.new.sum, tb_type_cat==3)
month.sum.4 <- subset(month.new.sum, tb_type_cat==4)


#There are only 4 seasons in a year (roughly), so k=4 here: (we include a random effect of year)
gam1.mon.re <- gam(monthly_cases ~ s(as.numeric(month), k=4, bs="cc") + s(year, bs="re"), family="poisson", data=month.sum.1)
gam2.mon.re <- gam(monthly_cases ~ s(as.numeric(month), k=4, bs="cc")+ s(year, bs="re"), family="poisson", data=month.sum.2)
gam3.mon.re <- gam(monthly_cases ~ s(as.numeric(month), k=4, bs="cc")+ s(year, bs="re"), family="poisson", data=month.sum.3)
gam4.mon.re <- gam(monthly_cases ~ s(as.numeric(month), k=4,  bs="cc")+ s(year, bs="re"), family="poisson", data=month.sum.4)

#These are the outputs we want to report (the p-values on the s()):
summary(gam1.mon.re) #monthly smoothing term is significant (meaning there is seasonality in the data) and year is somewhat sig (p=.06) as a random effect, menaing that are data have different patterns across the years
summary(gam2.mon.re) #monthly smoothing term is sig and random effect of year is also sig
summary(gam3.mon.re) #neither term significant
summary(gam4.mon.re) #neither term significant

#It can be messy to plot random effects, so rerun the model here without them to creat predictions for the plot
gam1.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4,  bs="cc") , family="poisson", data=month.sum.1)
gam2.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4,  bs="cc"), family="poisson", data=month.sum.2)
gam3.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4,  bs="cc"), family="poisson", data=month.sum.3)
gam4.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4,  bs="cc"), family="poisson", data=month.sum.4)

summary(gam1.mon) #yes sig seasonality
summary(gam2.mon) # yes sig seasonality
summary(gam3.mon) # yes sug seasonality
summary(gam4.mon) # no sig seasonality

#add predictions to monthly data so that you can plot
month.sum.1$gam_pred <- predict(gam1.mon, type="response", se.fit=T)$fit
month.sum.2$gam_pred <- predict(gam2.mon, type="response", se.fit=T)$fit
month.sum.3$gam_pred <- predict(gam3.mon, type="response", se.fit=T)$fit
month.sum.4$gam_pred <- predict(gam4.mon, type="response", se.fit=T)$fit

#and lower 95% confidence intervals
month.sum.1$gam_pred_lci <- predict(gam1.mon, type="response", se.fit=T)$fit - 1.96*predict(gam1.mon, type="response", se.fit=T)$se.fit
month.sum.2$gam_pred_lci <- predict(gam2.mon, type="response", se.fit=T)$fit - 1.96*predict(gam2.mon, type="response", se.fit=T)$se.fit
month.sum.3$gam_pred_lci <- predict(gam3.mon, type="response", se.fit=T)$fit - 1.96*predict(gam3.mon, type="response", se.fit=T)$se.fit
month.sum.4$gam_pred_lci <- predict(gam4.mon, type="response", se.fit=T)$fit - 1.96*predict(gam4.mon, type="response", se.fit=T)$se.fit

#and higher 95% confidence intervals
month.sum.1$gam_pred_uci <- predict(gam1.mon, type="response", se.fit=T)$fit + 1.96*predict(gam1.mon, type="response", se.fit=T)$se.fit
month.sum.2$gam_pred_uci <- predict(gam2.mon, type="response", se.fit=T)$fit + 1.96*predict(gam2.mon, type="response", se.fit=T)$se.fit
month.sum.3$gam_pred_uci <- predict(gam3.mon, type="response", se.fit=T)$fit + 1.96*predict(gam3.mon, type="response", se.fit=T)$se.fit
month.sum.4$gam_pred_uci <- predict(gam4.mon, type="response", se.fit=T)$fit + 1.96*predict(gam4.mon, type="response", se.fit=T)$se.fit

month.new.sum <- rbind(month.sum.1,month.sum.2, month.sum.3, month.sum.4)
month.new.sum <- arrange(month.new.sum, month, year)

month.new.sum$type_label = paste0("TB type-", month.new.sum$tb_type_cat)

#and plot with data
load("month.season.Rdata")

pTB1_yr <-ggplot(data=month.new.sum) + geom_point(aes(x=month, y=monthly_cases)) + facet_grid(type_label~., scales = "free_y") + xlab("") +
  scale_fill_manual(values=filz1) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  geom_line(aes(x=month, y=gam_pred)) + geom_ribbon(aes(x=month, ymin=gam_pred_lci, ymax=gam_pred_uci), alpha=.3, show.legend = F) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = month.season, aes(xmin=month_start, xmax=month_end, ymin=0, ymax=ymx, fill=season), alpha=.3) #+ coord_cartesian(ylim = c(0,55))

print(pTB1_yr)

ggsave(file = "Fig2_all_TB_type_1year.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)



#show that the seasonal patterns are mimicked across men and women, so that it is valid to look at seasonality in the cumulative dataset:
#summarize counts by type and by sex:
dat.sum.type = ddply(dat, .( tb_type_cat, sex, date), summarize, tot_case = length(ID) )
head(dat.sum.type)
dat.sum.type <- dat.sum.type[complete.cases(dat.sum.type),]
dat.sum.type$tb_type_cat = as.factor(dat.sum.type$tb_type_cat)
dat.sum.type$sex
dat.sum.type$sex[dat.sum.type$sex==0] = "F"
dat.sum.type$sex[dat.sum.type$sex==1] = "M"
dat.sum.type$sex = as.factor(dat.sum.type$sex)
dat.sum.type <- arrange(dat.sum.type, date, sex)

dat.sum.type$type_label = paste0("TB type-", dat.sum.type$tb_type_cat)

colz = c('F' = "magenta", 'M' = "blue")
filz = c('F' = "magenta", 'M' = "blue", "summer" = "yellow", "winter" = "seagreen")


pTB_all_long_data <- ggplot(data=dat.sum.type) + geom_line(aes(x=date, y=tot_case,  color=sex)) + facet_grid(type_label~., scales="free_y") +
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=ymx, fill=season), alpha=.3) 
print(pTB_all_long_data)

#We see that the data curves are mimicked in the first 3 seasonal datasets across men and women, so it is okay that we grouped them together as above.

ggsave(file = "Fig3all_TB_type_data_bysex.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)
