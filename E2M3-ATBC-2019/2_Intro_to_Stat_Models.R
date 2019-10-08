rm(list = ls())

## set wd
setwd("/Users/caraebrook/Documents/R/R_repositories/E2M2/E2M3-ATBC-2019/")
## Import your data.

library(ggplot2)

mad.pop <- cbind.data.frame(c(1970, 1980, 1990, 2000, 2010), c(6500000, 8900000, 11500000, 16100000, 21500000))
names(mad.pop) <- c("year", "population")
mad.pop

p1 <- ggplot(data=mad.pop) + geom_point(aes(x=year, y=population), size=4) + xlab("") + ylab("millions") +
  theme_bw() + theme(panel.grid = element_blank(), panel.background = element_rect(fill=NA), axis.text = element_text(size=14), axis.title = element_text(size=18))+
  scale_y_discrete(limits=c(0,5000000,10000000,15000000,20000000,25000000), labels=c(0,5,10,15,20,25))+
  coord_cartesian(ylim=c(0,25000000))
print(p1)

linear.regress= function(m,x,b){
  y = m*x + b
  return(y)
}


data= mad.pop
names(mad.pop) = c("x", "y")

sum.of.sqs = function(par,data){
  m = par[1] 
  b = par[2] 
  
  y.out <- linear.regress(m=m, x =data$x, b=b)
  
  sm.sq <- sum((data$y-y.out)^2)
  
  return(sm.sq)
}


parameter.out = optim(par = c(500000, 100), fn=sum.of.sqs, data=mad.pop)

y.out = linear.regress(m= 6469.822, x = mad.pop$x, b= 127597.578)

out.line = cbind.data.frame(mad.pop$x, y.out)
names(out.line) = c("x", "y")

p2 <- p1 + geom_line(data=out.line, aes(x,y), color="red")
print(p2)

R.linear.regression = lm(y~x, data=mad.pop)
R.predict= predict(R.linear.regression)

p3 = p2 + geom_line(aes(x= mad.pop$x, y = R.predict), color="green")
print(p3)


R.linear.regression.exp = lm(log(y)~(x), data=mad.pop)
R.predict.exp= exp(predict(R.linear.regression.exp))

p4 = p3 + geom_line(aes(x= mad.pop$x, y = R.predict.exp), color="blue")
print(p4)



get.sum.sq = function(par, data){
  m = exp(par[1])
  b = exp(par[2])
  out.mod = linear.regress(m=m, x=dat$x, b=b)
  
  sum.sq = sum((dat$y-out.mod)^2)
  return(sum.sq)
  
}

out.stat = optim(par=log(c(5500000,100)), fn=get.sum.sq, data=mad.pop)





y.out = linear.regress(m=500000, x=c(1970, 1980, 1990, 2000, 2010), b=-1000)

p1 + geom_line(aes(x=c(1970, 1980, 1990, 2000, 2010), y= y.out), col="red") + ylim(0,25)




stat.model = function(m,x,b){
  y = m*x+b
  return(y)
}


#and minimmize
out.stat = optim(par=log(c(1,100)), fn=get.sum.sq, data=dat)
exp(out.stat$par)
#run
out.mod = stat.model(m=exp(out.stat$par[1]), x=dat$x, b=exp(out.stat$par[2]))

#plot
plot(x=dat$x, y=dat$y, type="b")
lines(x=dat$x, y=out.mod, col="red")

#test with internal R script
lm(y~x, data=dat)
