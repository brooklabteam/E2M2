# we will use the diamonds dataset in the ggplot package to demonstrate for loops, if-else statements and creating function

#load ggplot2 package

library(ggplot2)

data<-diamonds

# let's say we need to check through this dataset to see if an entry meets certain condition

head(data)

#we want to check through he first 30 elements of this dataset we want the carat [weight], the price and the difference in price with the following diamond
#of the diamond and if a diamond has carats higher than 0.8, we say it is high quality
#we can check it manually and print the values and make the calculations ourselves in our dataset

for (i in 1:30){ ## This line creates a variable "i" that is 30 lines long
  print(paste(data$carat[i], data$price[i], (data$price[i+1]-data$price[i])))
  #Here we loop through and print each carat value of "i" 
  #as specified in the initial loop
  }



#############################################################################
## If-Else Statements 

## R also has a cool feature that allows you to specify certain actions if certain qualifications are met. These If-Else
## Statements take on the following general form:

## if (condition is TRUE) {
# do something
# } else {
# do different thing
# }

qc<-c()
for (i in 1:length(data$carat)){
  if(data$carat[i]>=0.8){
    qc<-append(qc,"high")}else{
      qc<-append(qc,"low-medium")}
}


#there is a function called ifelse() that does just the same thing
data$quality<-ifelse(data$carat>=0.8,"high","low-medium")
summary(data$quality)

###we want to characterize diamonds that have carat of 0.8 or more AND Premium has high quality, others will be low-medium
## to accomplish our precedent task, we can also use a for loop along with if-else statements
qc<-c()
for (i in 1:length(data$carat)){
  if(data$carat[i]>=0.8&data$cut[i]=="Premium"){
    qc<-append(qc,"high")}else{
      qc<-append(qc,"low-medium")}
}



summary(as.factor(qc))
#now what if we want to automate all of this and make it so that we only provide the "carat" and "cut" and it spits out those numbers?

q.check<-function(car,cu){
  qual<-c()
  for (i in 1:length(data$carat)){
    if(data$carat[i]>=car&data$cut[i]==cu){
      qual<-append(qual,"high")}else{
        qual<-append(qual,"low-medium")}
  }
  return(summary(as.factor(qual)))
}
q.check(3.0,"Premium")

######
#enough about diamonds let's continue writing the function about the coin toss

coin<-function(n,p)# we created a function called coin() where, we would need to specify the number of toss n and the probability of getting tail
{
  #the function will generate n numbers (either 0 or 1) with a p prob of getting tails.
  Tail<-rbinom(n,1,p)
  #count number of Tails
  numTail<-sum(Tail)
  return(numTail)
}

#let's try it out with a various number of toss and probabilities
coin(1000,0.5)

#### now we recruit x friends to help with an experiment
reps<-replicate(5000,coin(10,0.5)) #5000 of my best friends each toss the coin 10 times
hist(reps)
#what is the pobability of getting exactly 5 tails 
sum(reps==5)

for (i in 1:10) {print (i)}

##can you write a for loop to have the probabilities of getting exactly 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 tails

for(i in 1:10){
  print(sum(reps==i)/5000)}

reps0.3<-replicate(5000,coin(10,0.3)) #5000 of my best friends each toss the coin 10 times but this time their coin is weighted
hist(reps0.3, add=T)
