

## Let's work together to re-analyze Christian's dataset:


##############################################
# Load data and visualize
##############################################
## First, read in the .csv file.

## What do we do if we want to print out the Forearm length of every bat in the dataset?

##############################################
# For-loops
##############################################

## Can we use a for-loop to do this? Demonstrate.

## Remember the syntax:

## for (variable in vector) {
# do something
# }


##############################################
# If-Else Statements
##############################################

## What if we only want to print the forearm for the bats that are positive for infection?

## Remember the syntax:

## if (condition is TRUE) {
# do something
# } else {
# do different thing
# }

##############################################
# Functions
##############################################


## In the above section, we learned how to for-loop over our data and apply useful built-in functions in R (i.e. 'mean').
## Fidy showed you a list of some of these fucntions in his lecture this morning. Try to build a vector and use some of these
## functions to accomplish tasks in R. Calculate the min, max, and mean of the bat forearms in Christian's data, and
## tell us how many there are.





## Next, can you write an alternative function to the default in R that calculates the mean of your vector? Try it.
## Remember that you will first need to highlight and run your function to save it into your working directory, and
## you'll then be able to use it later.

## Remember the syntax:

## All functions in R take on the following general syntax:
#function_name <- function(arguments) {
# body
#}

## And, finally, can you write a third function to get the mean of the vector, but without using any built-in 
## functions in R, specifically sum()? (Hint: Try using a for-loop)





#############################################################################
## Apply Functions

## Many of the built-in functions in R are secretly running C-code (ususally)
## for-loops under the hood, but translating them to English so that you 
## don't actually have to think about it. One of the most powerful
## tools you can use in R, is the "apply" family of functions, which can be used
## to perform a general function over a list of sub-datasets in a linearized fashion.

## For example, what if we want to calculate the mean of female bat forearms and 
## the mean of male bat forearms separately?

## Try it with indexing.

## Try it with a for-loop.

## Try it by applying a function.




