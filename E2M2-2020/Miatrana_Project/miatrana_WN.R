#remove background
library(plyr)
library(dplyr)
library(ggplot2)
library(lubridate)

rm(list-ls())

#read in csv
dat = read.csv(file="WN.csv", header = T, stringsAsFactors = F)
head(dat)


#subset to variables of interest
dat.sub <- dplyr::select(dat, Ordre, Genre.et.especes, Region, Date.prlvs, A, SX, Res)
head(dat.sub)
names(dat.sub) <- c("order", "species", "region", "date",  "class", "sex", "res")
unique(dat.sub$sex) 
unique(dat.sub$class)
unique(dat.sub$order)
dat.sub$date = as.Date(dat.sub$date, format = "%m/%d/%y")

#check out 'risk factors' for testing seropos to WNV
dat.sub$season = as.factor(dat.sub$season)
dat.sub$species = as.factor(dat.sub$species)
dat.sub$class = as.factor(dat.sub$class)
dat.sub$sex = as.factor(dat.sub$sex)
mod1 <- glm(res~season+species+sex+class+sex, data=dat.sub, family = "binomial")
summary(mod1)

#no difference by season. One species implicated only. Hypsipetes madagascariensis - how many captured?
length(dat.sub$species[dat.sub$species=="Hypsipetes madagascariensis"]) # 14
length(dat.sub$species[dat.sub$species=="Hypsipetes madagascariensis" & dat.sub$res==1]) # 10
hyp.sub = subset(dat.sub, species=="Hypsipetes madagascariensis")

#try season by  order
dat.seas.order <- ddply(dat.sub, .(season, order), summarize, N_pos = sum(res), N_tot=length(res))
dat.seas.order$seroprevalence <- dat.seas.order$N_pos/dat.seas.order$N_tot
ggplot(data=dat.seas.order) + geom_point(aes(x=season, y=seroprevalence, color=order, size=N_tot))


#now summarize based on species
dat.spp <- ddply(dat.sub, .(order, species), summarize, N_pos=sum(res), N_tot=length(res), seroprevalence=N_pos/N_tot)
#and rank and plot by species
dat.spp <- arrange(dat.spp, desc(seroprevalence))
dat.spp$species = factor(dat.spp$species, levels=unique(dat.spp$species))

#and plot
ggplot(data=dat.spp) + geom_point(aes(x=species, y=seroprevalence, size=N_tot))

#or try as a bar plot - with N=X as the label for sample size
dat.spp$label = paste0("N=", dat.spp$N_tot) 
dat.spp$outline = 0
dat.spp$outline[dat.spp$species=="Hypsipetes madagascariensis"] = 1
dat.spp$outline = factor(dat.spp$outline)
colz= c('0'="black", '1'="red")


p1<-  ggplot(data=dat.spp) + geom_bar(aes(x=species, y=seroprevalence, fill=seroprevalence, color=outline), stat = "identity", size = 1.1, show.legend = F) +
                        geom_label(aes(x=species, y=seroprevalence+.05, label=label),label.size=0, size=1.9, fill=NA) + theme_bw() + scale_color_manual(values=colz) +
                        theme(axis.text.x = element_text(angle=90, face="italic"), axis.title.x = element_blank(), panel.grid = element_blank()) +
                        scale_fill_continuous(low="white", high="darkblue") + coord_cartesian(ylim=c(0,1.1), expand = F)
print(p1)  
  #and save
  ggsave(file = "seroprev_by_spp.pdf",
         units="mm",  
         width=70, 
         height=35, 
         scale=3, 
         dpi=300)


  
###################################################################################################
###################################################################################################
  
#other exploratory analyses that are not significant
  
### AGE  
#sum by species and age
dat.sum <-ddply(dat.sub, .(species, class),summarize, N_pos=sum(res), N=length(res))
dat.sum$seroprevalence = dat.sum$N_pos/dat.sum$N

#and plot
ggplot(data=dat.sum) + geom_point(aes(x=species, y=seroprevalence, color=class, size=N))

dat.slim = subset(dat.sum, N>9)

#too many species to see. just look at those species where you caught at least 10
ggplot(data=dat.slim) + geom_point(aes(x=species, y=seroprevalence, color=class, size=N))

#look at age-seroprev across all species summed
dat.tot= ddply(dat.sum, .(class), summarize, N_pos=sum(N_pos), N=sum(N))
dat.tot$seroprevalence = dat.tot$N_pos/dat.tot$N
  
ggplot(data=dat.tot) + geom_point(aes(x=class, y=seroprevalence, color=class, size=N))


### SEASON


#get month and year from date
dat.sub$year = year(dat.sub$date)
dat.sub$month = month(dat.sub$date)

#look at seroprevalence in wet vs. dry
dat.sub$season <- NA
dat.sub$season[dat.sub$month==1|dat.sub$month==11|dat.sub$month==12] <- "wet"
dat.sub$season[dat.sub$month==5] <- "dry"

dat.seas <- ddply(dat.sub, .(season, order, species), summarize, N_pos = sum(res), N_tot=length(res))
dat.seas$seroprevalence <- dat.seas$N_pos/dat.seas$N_tot
ggplot(data=dat.seas) + geom_point(aes(x=season, y=seroprevalence, color=species, size=N_tot))
#it looks like the seroprevalence is actually higher in the dry season
#try a boxplot
boxplot(seroprevalence~season, data=dat.seas, ylab = "seroprevalence by species")

