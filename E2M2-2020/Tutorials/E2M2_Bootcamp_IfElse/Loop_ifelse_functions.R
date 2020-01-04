# we will use the diamonds dataset in the ggplot package to demonstrate for loops, if-else statements and creating function

#load ggplot2 package

library(ggplot2)

data<-diamonds

# let's say we need to check through this dataset to see if an entry meets certain condition

head(data)

#we want to check through this dataset and if a diamond has cara higher than 0.8, we say it is high quality
#we can check it manually and print the carat of each diamond in our dataset

for (i in 1:30){ ## This line creates a variable "i" that is 20 lines long
  print(data$carat[i])#Here we loop through and print each carat value of "i" as specified in the initial loop
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




#now what if we want to automate all of this and make it so that we only provide the carat and cut and it spits out those numbers?

q.check<-function(car,cu){
  for (i in 1:length(data$carat)){
    if(data$carat[i]>=car&data$cut[i]==cu){
      qc<-append(qc,"high")}else{
        qc<-append(qc,"low-medium")}
  }
  return(qc)
}
q.check(2.0,"Premium")

#### finally we can write our own function to add this to our initial data frame

add.df<-function(df,vect){
  new.df<-cbind(df,vect)
  return(new.df)
}
data<-add.df(data,qc)
summary(data$vect)
