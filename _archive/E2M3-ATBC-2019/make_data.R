
count.host = c(100,500,700,600,900)
count.flea= c(1000,5000,7000,6000,9000)

dat = cbind.data.frame(count.host, count.flea)
names(dat)= c("host", "parasite")
head(dat)
dat$time = seq(1:length(dat$host))

with(dat, plot(time, parasite, col="red", type="l", ylim=c(0,10000), ylab="count"))
with(dat, lines(time, host, col="blue"))


pop.grow <- function(dat0, times, r){
  dat.host = rep(NA, length(times))
  dat.para = rep(NA, length(times))
  dat.host[1] = dat0
  dat.para[1] = dat0
  for (i in 1:(length(times)-1)){
    dat.host[i+1] = dat.host[i] + r*dat.host[i]
    dat.para[i+1] = dat.para[i] + r*dat.para[i]*dat.host[i]
  }
  
  
  return(list(dat.host, dat.para))
}

out <- pop.grow(dat0=6500000, times=seq(1,10,1), r=.0001)

out.dat = cbind.data.frame(out[[1]], out[[2]])
out.dat$time = seq(1,10,1)
names(out.dat) = c("host", "parasite", "time")

with(out.dat, plot(x=out.dat$time, y=out.dat$parasite, col="red", type="l", ylim=c(0,10000000000)))
with(out.dat, lines(x=out.dat$time, y=out.dat$host, col="blue"))












list.host <- c("bat", "bird", "fossa", "lemur")
list.host = rep(list.host, length = 100)

list.ecto <- c("flea", "mite", "tick")
list.ecto = rep(list.ecto, length = 100)

rand.ecto = sample(list.ecto, 100)
rand.host = sample(list.host, 100)

dat.rand = cbind.data.frame(rand.host, rand.ecto)
names(dat.rand) = c("host", "ecto")
head(dat.rand)
dat.rand$host = as.character(dat.rand$host)
dat.rand$ecto = as.character(dat.rand$ecto)
dat.rand$ecto[dat.rand$host=="bat"] <- "bat fly"

head(dat.rand)
dat.rand$season = 1
dat.rand1 = dat.rand

dat.rand$season = 2
dat.rand2 = dat.rand

dat.rand$season = 3
dat.rand$season=4

dat.rand3 = dat.rand
dat.rand4 = dat.rand

dat.rand = rbind(dat.rand1, dat.rand2, dat.rand3, dat.rand4)
write.csv(dat.rand, file = "E2M3_data.csv")
