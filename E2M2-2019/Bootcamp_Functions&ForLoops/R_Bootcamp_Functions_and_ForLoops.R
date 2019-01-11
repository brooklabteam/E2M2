library(plyr)
library(dplyr)
library(ggplot2)
### Writing Functions, For-Loops, and If-Else Statements
#############################################################################
##   Ecological and Epidemiological Modeling in Madagascar (E2M2)
##   Jan 12-21, 2019
##   Ranomafana, Fianarantsoa, Madaascar
#############################################################################
##   Cara E. Brook, 2019
##   UC Berkeley
#############################################################################

##  In this tutorial, we will learn how to write "for-loops" and set up our 
##  own functions in R. Both for-loops and functions offer a powerful way
##  to perform the same activity across many iterations of the data.
##  Let's start by loading the same dataset that Christian showed us earlier.

##  Do you remember how to read in a .csv file?

dat.e2m2 <- read.csv(file = "e2m2.csv", header=T, stringsAsFactors = F)


## Let's double-check that this is the right dataset
head(dat.e2m2)

## Notice that it doesn't matter that we named our dataset differently 
## than Christian - it still looks the same!

## Now, what it we would like to print out the Forearm of every bat in our
## dataset? Christian showed you that you can simply ask the dataset to print
## out a column, like this:
dat.e2m2$Forearm

## Or, equally:
dat.e2m2["Forearm"]

## But we can also use a for-loop and a function to do the same. R
## has many built-in functions, including one called "print" which will 
## print out data for you.

#############################################################################
## Writing for-loops:

## To loop through our data, we use the following syntax:

for (i in 1:length(dat.e2m2$Id)){ ## This line creates a variable "i" that is as long as the Id column in our dataset
  print(dat.e2m2$Forearm[i]) ## Here we loop through and print each value of "i" as specified in the initial loop
}

## You could also specify values for i to print only the first 10 numbers...

## For example:
for (i in 1:10){ ## This line creates a variable "i" that is as long as the Id column in our dataset
  print(dat.e2m2$Forearm[i]) ## Here we loop through and print each value of "i" as specified in the initial loop
}

## Or the second 10 numbers:
for (i in 10:20){ ## This line creates a variable "i" that is as long as the Id column in our dataset
  print(dat.e2m2$Forearm[i]) ## Here we loop through and print each value of "i" as specified in the initial loop
}

## The function "print" does not demonstrate the true power of a for-loop 
## because you could just as easily read out the data column. But for-loops 
## allow you to perform more complex tasks on the data. What if you wanted to 
## take the mean of every set of three forearms in your dataset?

## Try:
for (i in 1:length(dat.e2m2$Id)){ ## This line creates a variable "i" that is as long as the Id column in our dataset
  mean(dat.e2m2$Forearm[i], dat.e2m2$Forearm[i+1], dat.e2m2$Forearm[i+2]) ## Here we loop through and print each value of "i" as specified in the initial loop
}

## I get an error message! Why is this?
## Notice that the third Forearm you listed is indexed as dat.e2m2$Forearm[i+2]
## Since you specified i as 1:length(dat.e2m2$Id), this means that you are
## asking for the Forearm at the index of length(dat.e2m2$Id) + 2, which is:
length(dat.e2m2$Forearm) + 2 # 1202
## Out data columns are only actually 1200 long, so this index does not exist.

## Try, instead, to constrain your call to the for-loop:
for (i in 1:(length(dat.e2m2$Id)-2)){ ## Here, we constrain our largest number to the index which exists after 2 is added
  mean(dat.e2m2$Forearm[i], dat.e2m2$Forearm[i+1], dat.e2m2$Forearm[i+2]) 
} 
## The for-loop runs now, but we see no output. 

## Try printing the means instead:
for (i in 1:(length(dat.e2m2$Id)-2)){ ## Here, we constrain our largest number to the index which exists after 2 is added
  print(mean(dat.e2m2$Forearm[i], dat.e2m2$Forearm[i+1], dat.e2m2$Forearm[i+2]))
} 

## You can also save the means, but you need to initialize an empty list
## to hold them first:
out.mean <- list()
for (i in 1:(length(dat.e2m2$Id)-2)){ ## Here, we constrain our largest number to the index which exists after 2 is added
  out.mean[[i]] <- mean(dat.e2m2$Forearm[i], dat.e2m2$Forearm[i+1], dat.e2m2$Forearm[i+2])
} 
## And print the output list to see what it saved.
print(out.mean)

#############################################################################
## If-Else Statements

## R also has a cool feature that allows you to specify certain actions if certain qualifications are met. This can be
## paired powerfully with a for-loop. For instance, can you print the forearm length only of bats that are positive
## ("Pos=1") for the pathogen?
## Alternatively, if not positive, have it print the phrase: "Bat is pathogen negative."

for (i in 1:length(dat.e2m2$Id)){
  if(dat.e2m2$Pos[i]==1){
    print(dat.e2m2$Forearm)   
  }else{
    print("Bat is pathogen negative.")
  }
}

#############################################################################
## Writing Functions:

## In the above section, we learned how to for-loop over our data and apply useful built-in functions in R (i.e. 'mean').
## Fidy showed you a list of some of these in his lecture this morning. Try to build a vector and use some of these
## functions to accomplish tasks in R. Calculate the min, max, and mean of your vector and tell us how long it is.

example.vector= c(1,2,5,8,3,5,6,9)
min(example.vector) #1
max(example.vector) #9
mean(example.vector) #4.875

## In reality, R is simply hiding a lot of directional code 'under the hood' such that you can't see what is happening 
## but can simply call the function from it's name. In the case of mean, for example, the built-in R function first 
## adds all the components of the vector, then divides by the length of elements.

## All functions in R take on the following general syntax:
#function_name <- function(arguments) {
 # body
#}

## You can name them whatever you like and specify as many arguments as are needed.
## Can you write an alternative function to the default in R that calculates the mean of your vector? Try it.
## Remember that you will first need to highlight and run your function to save it into your working directory, and
## you'll then be able to use it later.

## Load your new function:
alternative.mean <- function(vector){
  sum.of.elements <- sum(vector)
  length.of.elements <- length(vector)
  alt.mean <- sum.of.elements/length.of.elements
  return(alt.mean)
}

## Run your new function on your vector:
alternative.mean(vector=example.vector) #4.875

## CHALLENGE: Can you write a third function to get the mean of the vector, but without using any built-in 
## functions in R, specifically sum() and length()? (Hint: Try using a for-loop)

alternative.mean2 <- function(vector){
  new.vector <- vector
  for (i in 2:(length(vector))){
    #sum the two
    tot.sum <- new.vector[i] + new.vector[i-1]
    
    #write over the current entry for the next step
    new.vector[i] <- tot.sum
    }    
    sum.of.elements <- new.vector[length(vector)]
    
    length.of.elements <- length(vector)
    
    alt.mean2 <- sum.of.elements/length.of.elements
    
    return(alt.mean2)}

alternative.mean2(vector=example.vector) #4.875


#############################################################################
## The Power of R: Apply Functions

## Many of the built-in functions in R are secretly running C-code (ususally)
## for-loops under the hood, but translating them to English so that you 
## don't actually have to think about it. One of the most powerful
## tools you can use in R, is the "apply" family of functions, which can be used
## to perform a general function over a list of sub-datasets in a linearized fashion.

## For example, what if we want to calculate the mean of female bat forearms and 
## the mean of male bat forearms separately?

## Christian showed you that you can do this with indexing:
mean(dat.e2m2$Forearm[dat.e2m2$Sex == "f"]) # females = 70.81329
mean(dat.e2m2$Forearm[dat.e2m2$Sex == "m"]) # males = 78.09067

## How might you determine this using a for-loop?
list.females = list()
list.males = list()
for(i in 1:length(dat.e2m2$Id)){
  if(dat.e2m2$Sex[i]=="f"){
    list.females[[i]] = dat.e2m2$Forearm[i] 
  }else{
    list.males[[i]] = dat.e2m2$Forearm[i] 
  }
}

teach.fem = sum(c(unlist(list.females)))/ length(list.females)
teach.male = sum(c(unlist(list.males)))/ length(list.males)

## Or you could use an apply function.
## First, separate data into a list of two: m and f
sep.dat.list = dlply(dat.e2m2, .(Sex))

## You now have a list with f and m data. Apply the mean function 
## over the list. 

out.MF <- lapply(sep.dat.list, get.mean.forearm) 

## new to write a function
get.mean.forearm <- function(data){
  data2 = mean(data$Forearm)
  return(data2)
}



