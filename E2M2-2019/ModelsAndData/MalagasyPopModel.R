#make data frame
library(ggplot2)

mad.pop <- cbind.data.frame(c(1970, 1980, 1990, 2000, 2010), c(6500000, 8900000, 11500000, 16100000, 21500000))
names(mad.pop) <- c("year", "population")
mad.pop

p1 <- ggplot(data=mad.pop) + geom_point(aes(x=year, y=population), size=4) + xlab("") + ylab("millions") +
      theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill=NA), axis.text = element_text(size=14), axis.title = element_text(size=18))+
      scale_y_discrete(limits=c(0,5000000,10000000,15000000,20000000,25000000), labels=c(0,5,10,15,20,25))+
      coord_cartesian(ylim=c(0,25000000))
print(p1)

ggsave(file = "mad_pop.pdf",
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)

#then, with a linear regression:

m1 <- lm(population~year, data=mad.pop)
mad.pop$linear_pred <- predict(m1)
summary(m1)

p2 <- ggplot(data=mad.pop) + geom_point(aes(x=year, y=population), size=4) + xlab("") + ylab("millions") +
  geom_line(aes(x=year, y= linear_pred), col="red", size=1) +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill=NA), axis.text = element_text(size=14), axis.title = element_text(size=18))+
  scale_y_discrete(limits=c(0,5000000,10000000,15000000,20000000,25000000), labels=c(0,5,10,15,20,25))+
  coord_cartesian(ylim=c(0,25000000))
print(p2)

ggsave(file = "mad_pop_linear.pdf",
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)

#and a non-linear regression (exponential)
m2 <- lm(log(mad.pop$population)~mad.pop$year)
mad.pop$exp_pred <- exp(predict(m2))
summary(m2)

p3 <- ggplot(data=mad.pop) + geom_point(aes(x=year, y=population), size=4) + xlab("") + ylab("millions") +
  geom_line(aes(x=year, y= exp_pred), col="blue", size=1) +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill=NA), axis.text = element_text(size=14), axis.title = element_text(size=18))+
  scale_y_discrete(limits=c(0,5000000,10000000,15000000,20000000,25000000), labels=c(0,5,10,15,20,25))+
  coord_cartesian(ylim=c(0,25000000))
print(p3)


ggsave(file = "mad_pop_exponential.pdf",
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)

#and make your mechanistic model:
pop.grow <- function(dat0, times, r){
  dat.proj = rep(NA, length(times))
  dat.proj[1] = dat0
  for (i in 1:(length(times)-1)){
    dat.proj[i+1] = dat.proj[i] + r*dat.proj[i]
  }
  return(dat.proj)
}
lst.sq <- function(par, real.dat, times){
  real.dat$predicted <- pop.grow(dat0=real.dat$population[1], times= real.dat$year, r=par)
  real.dat$diff_sq <- real.dat$population-real.dat$predicted
  real.dat$diff_sq <- (real.dat$diff_sq)^2
  sm.sqs = sum(real.dat$diff_sq)
  return(sm.sqs)
  
}
wrap.fit <- function(mad.dat, par.guess){
  out <- optim(par = par.guess, fn=lst.sq, real.dat=mad.dat, times=mad.dat$year)
  out.fit <- c(out$par, out$value)
  names(out.fit) <- c("r", "sum_sqs")
  return(out.fit)
  
}

#and run it
out.fit <- wrap.fit(mad.dat = mad.pop, par.guess = .3)

#and run and plot:
mad.pop$mech_pred <- pop.grow(dat0=6500000, times=mad.pop$year, r=out.fit['r'])

p4 <- ggplot(data=mad.pop) + geom_point(aes(x=year, y=population), size=4) + xlab("") + ylab("millions") +
  geom_line(aes(x=year, y= mech_pred), col="green", size=1) +
  theme_bw() + theme(panel.grid = element_blank(), axis.text = element_text(size=14), axis.title = element_text(size=18))+
  scale_y_discrete(limits=c(0,5000000,10000000,15000000,20000000,25000000), labels=c(0,5,10,15,20,25))+
  coord_cartesian(ylim=c(0,25000000))
print(p4)

ggsave(file = "mad_pop_mechanistic.pdf",
       units="mm",  
       width=50, 
       height=40, 
       scale=3, 
       dpi=300)






#blank background
p4 <- ggplot(data=mad.pop) + geom_point(aes(x=year, y=population), size=4, color="white") + xlab("") + ylab("millions") +
  geom_line(aes(x=year, y= mech_pred), col="green", size=1) +
  theme_bw() + theme(panel.grid = element_blank(),panel.border = element_rect(color="white"), axis.ticks = element_line(color="white"), panel.background = element_rect(fill=NA), axis.text = element_text(size=14, color="white"), axis.title = element_text(size=18, color="white"), plot.background = NULL)+
  scale_y_discrete(limits=c(0,5000000,10000000,15000000,20000000,25000000), labels=c(0,5,10,15,20,25))+
  coord_cartesian(ylim=c(0,25000000))
print(p4)