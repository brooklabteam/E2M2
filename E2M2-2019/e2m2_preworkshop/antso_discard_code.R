#make season dat
season=rep(c("winter", "summer"), each=5)
date_start = c("2010-04-01", "2011-04-01", "2012-04-01", "2013-04-01", "2014-04-01", "2010-10-01", "2011-10-01", "2012-10-01", "2013-10-01", "2014-10-01")
date_end = c("2010-09-01", "2011-09-01", "2012-09-01", "2013-09-01", "2014-09-01", "2011-03-01", "2012-03-01", "2013-03-01", "2014-03-01", "2015-03-01")
season.dat = cbind.data.frame(season, date_start, date_end)
season.dat$date_start = as.Date(as.character(season.dat$date_start))
season.dat$date_end = as.Date(as.character(season.dat$date_end))


season.dat$type_label = "TB type-1"
season.dat$ymx = 55
season.dat7 <- season.dat
season.dat8 <- season.dat
season.dat9 <- season.dat
season.dat7$type_label = "TB type-2"
season.dat7$ymx = 25
season.dat8$type_label = "TB type-3"
season.dat8$ymx = 25
season.dat9$type_label = "TB type-4"
season.dat9$ymx = 10


season.dat <- rbind(season.dat, season.dat7, season.dat8, season.dat9)


season.dat2 = season.dat
season.dat2 = season.dat2[5:7,]
season.dat2$month_start = month(season.dat2$date_start)
season.dat2$month_end = month(season.dat2$date_end)
season.dat2$month_end[2] = 12
season.dat2$month_start[3] = 1


season.dat2$type_label = "TB type-1"
season.dat2$ymx = 55
season.dat3 <- season.dat2
season.dat4 <- season.dat2
season.dat5 <- season.dat2
season.dat3$type_label = "TB type-2"
season.dat3$ymx = 25
season.dat4$type_label = "TB type-3"
season.dat4$ymx = 25
season.dat5$type_label = "TB type-4"
season.dat5$ymx = 10

season.dat2 <- rbind(season.dat2, season.dat3, season.dat4,season.dat5)



#summarize case counts by month
dat.sum = ddply(dat, .(year, month, date), summarize, tot_case = length(ID) )
head(dat.sum)





ggplot(data=dat.sum) + geom_line(aes(x=date, y=tot_case)) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=150, fill=season), alpha=.3) + ylim(c(0,150))



#GAM modeling by sex within a year

dat.sum.type <- dat.sum.type[complete.cases(dat.sum.type),]
dat.sum.1 <- subset(dat.sum.type, tb_type_cat==1)
dat.sum.2 <- subset(dat.sum.type, tb_type_cat==2)
dat.sum.3 <- subset(dat.sum.type, tb_type_cat==3)
dat.sum.4 <- subset(dat.sum.type, tb_type_cat==4)

#do a separate gam for each of the 4 types of TB. separate smoothers by sex
gam1 <- gam(tot_case ~ s(as.numeric(date), k=7, by = sex, bs="cr"), family="poisson", data=dat.sum.1)
gam2 <- gam(tot_case ~ s(as.numeric(date), k=7, by = sex, bs="cr"), family="poisson", data=dat.sum.2)
gam3 <- gam(tot_case ~ s(as.numeric(date), k=7, by = sex, bs="cr"), family="poisson", data=dat.sum.3)
gam4 <- gam(tot_case ~ s(as.numeric(date), k=7, by = sex, bs="cr"), family="poisson", data=dat.sum.4)
summary(gam1)
summary(gam2)
summary(gam3)
summary(gam4) # this seasonal smoother (tb type) is not significant, meaning type 4 is not significant

#now add the predictions to the datasets
dat.sum.1$gam_pred <- predict(gam1, type="response", se.fit=T)$fit
dat.sum.2$gam_pred <- predict(gam2, type="response", se.fit=T)$fit
dat.sum.3$gam_pred <- predict(gam3, type="response", se.fit=T)$fit
dat.sum.4$gam_pred <- predict(gam4, type="response", se.fit=T)$fit


dat.sum.1$gam_pred_lci <- predict(gam1, type="response", se.fit=T)$fit - 1.96*predict(gam1, type="response", se.fit=T)$se.fit
dat.sum.2$gam_pred_lci <- predict(gam2, type="response", se.fit=T)$fit - 1.96*predict(gam2, type="response", se.fit=T)$se.fit
dat.sum.3$gam_pred_lci <- predict(gam3, type="response", se.fit=T)$fit - 1.96*predict(gam3, type="response", se.fit=T)$se.fit
dat.sum.4$gam_pred_lci <- predict(gam4, type="response", se.fit=T)$fit - 1.96*predict(gam4, type="response", se.fit=T)$se.fit


dat.sum.1$gam_pred_uci <- predict(gam1, type="response", se.fit=T)$fit + 1.96*predict(gam1, type="response", se.fit=T)$se.fit
dat.sum.2$gam_pred_uci <- predict(gam2, type="response", se.fit=T)$fit + 1.96*predict(gam2, type="response", se.fit=T)$se.fit
dat.sum.3$gam_pred_uci <- predict(gam3, type="response", se.fit=T)$fit + 1.96*predict(gam3, type="response", se.fit=T)$se.fit
dat.sum.4$gam_pred_uci <- predict(gam4, type="response", se.fit=T)$fit + 1.96*predict(gam4, type="response", se.fit=T)$se.fit




#and plot pred with data
tot.dat.sum <- rbind(dat.sum.1, dat.sum.2, dat.sum.3, dat.sum.4)
tot.dat.sum$tb_type_cat = as.factor(tot.dat.sum$tb_type_cat)
tot.dat.sum$sex[tot.dat.sum$sex==0] = "F"
tot.dat.sum$sex[tot.dat.sum$sex==1] = "M"
tot.dat.sum$sex = as.factor(tot.dat.sum$sex)
head(tot.dat.sum)

tot.dat.sum$type_label = paste0("TB type-", tot.dat.sum$tb_type_cat)





#and plot these all together - both longitudinally and within a year

pTB_all_long <- ggplot(data=tot.dat.sum) + geom_line(aes(x=date, y=tot_case,  color=sex)) + facet_grid(type_label~., scales="free_y") +
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  #geom_line(aes(x=date, y=gam_pred,  color=sex)) + geom_ribbon(aes(x=date, ymin=gam_pred_lci, ymax=gam_pred_uci, fill=sex), alpha=.3, show.legend = F) + 
  ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=ymx, fill=season), alpha=.3) #+ coord_cartesian(ylim = c(0,55))
print(pTB_all_long)

ggsave(file = "all_TB_type_longitudinal_bysex.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)











month.dat.sum = ddply(tot.dat.sum, .(tb_type_cat, sex, year, month), summarize, monthly_cases=sum(tot_case))
head(month.dat.sum)
month.dat.sum <- arrange(month.dat.sum, month, tb_type_cat)
month.sum.1 <- subset(month.dat.sum, tb_type_cat==1)
month.sum.2 <- subset(month.dat.sum, tb_type_cat==2)
month.sum.3 <- subset(month.dat.sum, tb_type_cat==3)
month.sum.4 <- subset(month.dat.sum, tb_type_cat==4)


gam1.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc") + s(year, bs="re"), family="poisson", data=month.sum.1)
gam2.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc")+ s(year, bs="re"), family="poisson", data=month.sum.2)
gam3.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc")+ s(year, bs="re"), family="poisson", data=month.sum.3)
gam4.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc")+ s(year, bs="re"), family="poisson", data=month.sum.4)

#no re on year (for plotting purposes only)
gam1.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc") , family="poisson", data=month.sum.1)
gam2.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc"), family="poisson", data=month.sum.2)
gam3.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc"), family="poisson", data=month.sum.3)
gam4.mon <- gam(monthly_cases ~ s(as.numeric(month), k=4, by = sex, bs="cc"), family="poisson", data=month.sum.4)

summary(gam1.mon) #men are seasonal; women marginally so; sig random effect of year
summary(gam2.mon) #women are seasonal; men not so; sig random effect of year
summary(gam3.mon) #women are seasonal; men not so; sig random effect of year
summary(gam4.mon) # no sig seasonality

month.sum.1$gam_pred <- predict(gam1.mon, type="response", se.fit=T)$fit
month.sum.2$gam_pred <- predict(gam2.mon, type="response", se.fit=T)$fit
month.sum.3$gam_pred <- predict(gam3.mon, type="response", se.fit=T)$fit
month.sum.4$gam_pred <- predict(gam4.mon, type="response", se.fit=T)$fit


month.sum.1$gam_pred_lci <- predict(gam1.mon, type="response", se.fit=T)$fit - 1.96*predict(gam1.mon, type="response", se.fit=T)$se.fit
month.sum.2$gam_pred_lci <- predict(gam2.mon, type="response", se.fit=T)$fit - 1.96*predict(gam2.mon, type="response", se.fit=T)$se.fit
month.sum.3$gam_pred_lci <- predict(gam3.mon, type="response", se.fit=T)$fit - 1.96*predict(gam3.mon, type="response", se.fit=T)$se.fit
month.sum.4$gam_pred_lci <- predict(gam4.mon, type="response", se.fit=T)$fit - 1.96*predict(gam4.mon, type="response", se.fit=T)$se.fit


month.sum.1$gam_pred_uci <- predict(gam1.mon, type="response", se.fit=T)$fit + 1.96*predict(gam1.mon, type="response", se.fit=T)$se.fit
month.sum.2$gam_pred_uci <- predict(gam2.mon, type="response", se.fit=T)$fit + 1.96*predict(gam2.mon, type="response", se.fit=T)$se.fit
month.sum.3$gam_pred_uci <- predict(gam3.mon, type="response", se.fit=T)$fit + 1.96*predict(gam3.mon, type="response", se.fit=T)$se.fit
month.sum.4$gam_pred_uci <- predict(gam4.mon, type="response", se.fit=T)$fit + 1.96*predict(gam4.mon, type="response", se.fit=T)$se.fit

month.dat.sum <- rbind(month.sum.1,month.sum.2, month.sum.3, month.sum.4)

month.dat.sum$type_label = paste0("TB type-", month.dat.sum$tb_type_cat)




#and plot these all together - both longitudinally and within a year
pTB_all_long <- ggplot(data=tot.dat.sum) + geom_point(aes(x=date, y=tot_case,  color=sex)) + facet_grid(type_label~., scales="free_y") +
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  geom_line(aes(x=date, y=gam_pred,  color=sex)) + geom_ribbon(aes(x=date, ymin=gam_pred_lci, ymax=gam_pred_uci, fill=sex), alpha=.3, show.legend = F) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=ymx, fill=season), alpha=.3) #+ coord_cartesian(ylim = c(0,55))
print(pTB_all_long)

ggsave(file = "all_TB_type_longitudinal_bysex.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)


#and within a year
#and look within a single year
pTB1_yr <-ggplot(data=month.dat.sum) + geom_point(aes(x=month, y=monthly_cases,  color=sex)) + facet_grid(type_label~., scales = "free_y") +
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  geom_line(aes(x=month, y=gam_pred,  color=sex)) + geom_ribbon(aes(x=month, ymin=gam_pred_lci, ymax=gam_pred_uci, fill=sex), alpha=.3, show.legend = F) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat2, aes(xmin=month_start, xmax=month_end, ymin=0, ymax=ymx, fill=season), alpha=.3) #+ coord_cartesian(ylim = c(0,55))

print(pTB1_yr)

ggsave(file = "all_TB_type_1year_by_sex.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)


#and try not separated by gender, since seasonality appears to correspond:




pTB1_yr <-ggplot(data=month.dat.sum) + geom_point(aes(x=month, y=monthly_cases)) + facet_grid(type_label~., scales = "free_y") +
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  geom_line(aes(x=month, y=gam_pred)) + geom_ribbon(aes(x=month, ymin=gam_pred_lci, ymax=gam_pred_uci), alpha=.3, show.legend = F) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat2, aes(xmin=month_start, xmax=month_end, ymin=0, ymax=ymx, fill=season), alpha=.3) #+ coord_cartesian(ylim = c(0,55))

print(pTB1_yr)

ggsave(file = "all_TB_type_1year_by_sex.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)







##TB type 1

#and add the model to each plot, one-by-one
pTB1_long <- ggplot(data=tot.dat.sum) + geom_point(aes(x=date, y=tot_case,  color=sex)) + facet_grid(type_label~.) +
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.title = element_blank())+
  geom_line(aes(x=date, y=gam_pred,  color=sex)) + geom_ribbon(aes(x=date, ymin=gam_pred_lci, ymax=gam_pred_uci, fill=sex), alpha=.3) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=ymx, fill=season), alpha=.3) + coord_cartesian(ylim = c(0,55))

print(pTB1_long)

##TB type 2
pTB2_long <- ggplot(data=subset(tot.dat.sum, tb_type_cat==2)) + geom_point(aes(x=date, y=tot_case,  color=sex)) + 
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank())+
  geom_line(aes(x=date, y=gam_pred,  color=sex)) + geom_ribbon(aes(x=date, ymin=gam_pred_lci, ymax=gam_pred_uci, fill=sex), alpha=.3) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=60, fill=season), alpha=.3) + coord_cartesian(ylim = c(0,55))

print(pTB2_long)


#and look within a single year
pTB1_yr <-ggplot(data=subset(month.dat.sum, tb_type_cat==1)) + geom_point(aes(x=month, y=monthly_cases,  color=sex)) + 
  scale_color_manual(values=colz) + scale_fill_manual(values=filz) + theme_bw() + theme(panel.grid = element_blank())+
  geom_line(aes(x=month, y=gam_pred,  color=sex)) + geom_ribbon(aes(x=month, ymin=gam_pred_lci, ymax=gam_pred_uci, fill=sex), alpha=.3) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat2, aes(xmin=month_start, xmax=month_end, ymin=0, ymax=150, fill=season), alpha=.3) + coord_cartesian(ylim = c(0,55))

print(pTB1_yr)











#and look by gender
dat.sum.sex = ddply(dat, .(year, month, date, sex), summarize, tot_case = length(ID) )
head(dat.sum.sex)
dat.sum.sex = dat.sum.sex[complete.cases(dat.sum.sex),]
dat.sum.sex$sex[dat.sum.sex$sex==0] = "F"
dat.sum.sex$sex[dat.sum.sex$sex==1] = "M"

ggplot(data=dat.sum.sex) + geom_line(aes(x=date, y=tot_case, color=sex)) + ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=150, fill=season), alpha=.3) + ylim(c(0,150))


ggsave(file = "Antso_seasonal_TB_by_sex.pdf",
       units="mm",  
       width=80, 
       height=60, 
       scale=3, 
       dpi=300)


#now take a look at age data
head(dat)
unique(dat$age)
#round to 1
dat$age = round(dat$age, 0)

dat$sex = as.numeric(as.character(dat$sex))
dat$sex[dat$sex==0] = "F"
dat$sex[dat$sex==1] = "M"
dat <- dplyr::select(dat, -(age_month))
dat <- dat[complete.cases(dat),]

#summarize by age
dat.age = ddply(dat, .(age, sex), summarize, tot_case = length(ID))
ggplot(data=dat.age) + geom_point(aes(x=age, y=tot_case, col=sex), size=3)

ggsave(file = "Antso_TB_by_age_sex.pdf",
       units="mm",  
       width=40, 
       height=40, 
       scale=3, 
       dpi=300)


#try seasonality by type of TB
dat.TB.seas = ddply(dat, .(year, month, date, tb_type_cat), summarize, tot_case = length(ID) )
head(dat.TB.seas)
dat.TB.seas$tb_type_cat = as.factor(dat.TB.seas$tb_type_cat)

ggplot(data=dat.TB.seas) + geom_line(aes(x=date, y=tot_case, color=tb_type_cat)) + #ylab("total TB cases - one hospital Tana") +
  geom_rect(data = season.dat, aes(xmin=date_start, xmax=date_end, ymin=0, ymax=90, fill=season), alpha=.3) + ylim(c(0,90))

ggsave(file = "TB_type_seasonality.pdf",
       units="mm",  
       width=40, 
       height=40, 
       scale=3, 
       dpi=300)

#try a gam to look at seasonality - start just with pos microscopy cases
pos.mic <- subset(dat.TB.seas, tb_type_cat==1)
pos.mic$month <- month(pos.mic$date)

gam.month <- gam(tot_case ~ s(month, bs="cc") + s(year, bs="re"), data=pos.mic)
summary(gam.month)
pos.mic$predict_case = predict(gam.month)

#and plot as month
ggplot(data=pos.mic) + geom_point(aes(x=month, y=tot_case), pch=1, size=3) + geom_line(aes(x=month, y=predict_case), col="red", size=3) +
  geom_rect(data = season.dat2, aes(xmin=month_start, xmax=month_end, ymin=0, ymax=70, fill=season), alpha=.3) + ylab("positive microscopy TB cases - one hospital Tana") +
  xlab("month of year") + coord_cartesian(xlim=c(1,12), expand = F)

ggsave(file = "seasonal_gam.pdf",
       units="mm",  
       width=40, 
       height=40, 
       scale=3, 
       dpi=300)
