#remove any background in your environment
rm(list=ls())

#load a few helpful packages
library(ggplot2)
library(cowplot)
library(plyr)
library(dplyr)
library(lubridate)

#check wd is this folder
getwd()

#load data
dat <- read.csv(file = "E2M2_angonoka_pressure_fire_05 01 2020.csv", header = T, stringsAsFactors = F)
#look at data
head(dat)

#rename columns
names(dat) <- c("site", "date", "number_patrollers", "count_tortoise", "count_indirect_sign", "count_fire")

#look again
head(dat)

#rename sites with spaces
unique(dat$site)
dat$site[dat$site=="Bloc de l'ouest"] <- "Ouest"

#replace NAs with 0
dat$count_tortoise[is.na(dat$count_tortoise)] = 0
dat$count_indirect_sign[is.na(dat$count_indirect_sign)] = 0
dat$count_fire[is.na(dat$count_fire)] = 0

#reformat date to default in R
dat$date = as.Date(dat$date, format = "%d/%m/%Y")

#and look at data
#tortoise
p1 <- ggplot(data=dat) + geom_point(aes(x=date, y=count_tortoise, color=site))
print(p1)
#most observations appear to be 1 or 2 tortoises

#indirect sign
p2 <- ggplot(data=dat) + geom_point(aes(x=date, y=count_indirect_sign, color=site))
print(p2)

#fire
p3 <- ggplot(data=dat) + geom_point(aes(x=date, y=count_fire, color=site))
print(p3)

#plot together
cowplot::plot_grid(p1,p2,p3, nrow=3)
#not particularly informative

#calculate the mean number of tortoises/indirect sign/fire/ per patrol per month per site
#first, add a month column
dat$month = month(dat$date)
head(dat)
sum.dat <- ddply(dat, .(site, month), summarize, mean_tortoise =mean(count_tortoise), mean_sign = mean(count_indirect_sign), mean_fire=mean(count_fire), N_patrol = length(unique(date)))
head(sum.dat)

#make your data so that you can facet by count
tort.dat <- dplyr::select(sum.dat, -(mean_sign), -(mean_fire))
tort.dat$type = "tortoise"

sign.dat <- dplyr::select(sum.dat, -(mean_tortoise), -(mean_fire))
sign.dat$type = "indirect_sign"

fire.dat <- dplyr::select(sum.dat, -(mean_tortoise), -(mean_sign))
fire.dat$type = "fire"

#synchronize names
names(tort.dat) <- names(sign.dat) <- names(fire.dat) <- c("site", "month", "count_per_patrol", "N_patrol", "type")

#and bind
sum.dat.long <- rbind(tort.dat, sign.dat, fire.dat)

head(sum.dat.long)

#and redo the same plots as above - but faceted
# also make the point size correspond to the sample size (number of patrols)
p4 <- ggplot(data=sum.dat.long) + geom_point(aes(x=month, y=count_per_patrol, color=site, size = N_patrol)) + facet_grid(type~.)
print(p4)
#add lines
p5 <- p4 + geom_line(aes(x=month, y=count_per_patrol, color=site))
print(p5)

#try aggregrating across all sites
tot.sum <- ddply(sum.dat.long, .(month, type), summarize, count_per_patrol = sum(count_per_patrol), N_patrol=sum(N_patrol))
head(tot.sum)

#and plot
p6 <- ggplot(data=tot.sum) + geom_point(aes(x=month, y=count_per_patrol, size = N_patrol)) + facet_grid(type~.)
print(p6)

#add a line
p7 <- p5 + geom_line(aes(x=month, y=count_per_patrol))
print(p7)
#much easier to visualize! probably best to use the aggregated data
#edit the axes
p8 <- p7 + ylab("count per patrol") + 
      scale_x_discrete(limits=c(1,2,3,4,5,6,7,8,9,10,11,12), labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sept", "Oct", "Nov", "Dec")) +
      theme_bw() + theme(panel.grid = element_blank(), axis.title.x = element_blank())
print(p8)

#save this plot - will save autmoatically to your working directory
ggsave(file = "Liantsoa_cumulative_plot.pdf",
       units="mm",  
       width=40, 
       height=50, 
       scale=3, 
       dpi=300)

