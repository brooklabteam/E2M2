################################################################################
#-------------------------------------------------------- R BASICS          ###
################################################################################

################################################################################
### MODULE 1: A Sample R Session                                             ### 
################################################################################
## OBJECTIVE:
## The goal of this exercise is to become accustomed to R and the way it 
## responds to line commands. 


## PART 1 ##
# Here you will be prelimilarly introduced to several important concepts that we 
# will cover in more detail during the rest of the workshop


# To access on-line manuals, references, and other material, you can type:
help.start()


### A. WHAT IS AN OBJECT AND HOW TO CREATE IT ##################################

# Simply, an object in R is a piece of memory with a name. Objects contains some
# information, which is often used for figures or analyses.

# Objects can have multiple types of data

#TYPES OF DATA WITHIN OBJECTS ############################################

# There are 4 types of data commonly used in R:
# 1. Numerical - Numbers - e.g. 5, 0.4, 3.885
# 2. Characters - text - e.g. "a", "workshop", "R is awesome"
# 3. Logical - i.e. TRUE, FALSE
# 4. Special Values - e.g. NA, Inf, NULL, NaN

# The easiest way to create an object is using the arrow operator *<-*

x <- 10 # This creates an object called *x* that contains the value 10

# Typing the object name returns the information stored inside the object 
# by printing it in the console
x

y <- 5.3 # This created another object called *y* that contains the value 5.3
y


# As you might have already noticed, R does not run lines that start with '#'.

# This is used to create commentaries. 

# # You can store the ouput of a function in an object. Then, these output can be 
# used for other purposes: (more on functions below)
seq(from = 4, to = 10, by = 2) # This just prints the result in the screen

x <- seq(from=4, to=10, by=2) # This RE-WRITES the object *x* with the sequence

x 
# Now the object *x* contains the sequence from 4 to 10 by 2.





# You now can do things with the values stored in 'x'. 
#For example you can do basic arithmetics: 
# the arithmetic operators for sum '+', substraction '-', multiplication '*', and
# division '/'.
#we can multiply the values in 'x' by a constant:
x*2
#add pi to each value of x
x+pi

#other arithmetic function: "is x equal to, greater than, less than  ...?" 
x
x==6
x>6
x>=6


#or calculate other statistics
# function 'mean' to calculate the mean of the values in 'x':
mean(x)

# You can calculate other statistics, like the standard deviation, or just create
# a summary of the values in x:
sd(x)
summary(x)

# You can change the order of the values:
sort(x)
sort(x,decreasing=T)


# We can also write a single line of code that performs multiple actions and 
# saves the output, for exmple
y <- rnorm(4)*2

# This creates 4 random values from a normal distribution, then multiplies
# each value by 2 and finally stores the results in an object named 'y'. 



#To see the values in 'y', just type:
y

plot(x,y, xlim=c(0,12))

############################################################################################
#---------------------------------------- PART 2-------------------------------------------#
############################################################################################

# In this second part, you will continue to play with various elements of R.
# and get accustomed to writing commands.
# It would be best if you try to type the code into the console rather than just copy-paste:

### A. WRITING IN R IS SIMILAR TO WRITING IN ENGLISH ###########################

# A command in English: 

# "Step three times forward"

# 1. *Step* is the word that defines the requested action (verb)
# 2. *three times* and *forward* are two modifiers of that action


# Another command in English:

# "Generate a sequence from 5 to 20 with values every 0.5"

# 1. *Generate a sequence* is what defines the action to be taken
# 2. *From 5* modifies the action defining where to start
# 3. *to 20* defines the end of the sequence
# 4. *with values every 0.5 * defines the spacing of values in the sequence


# How do you give the above command to a computer using R?

seq(from = 5, to = 20, by = 0.5)

# 1. *seq* is the name of a function, and defines an action to be taken in R.
# 2. *from*: is the name of an argument that defines the beginning of the 
#    sequence, In this case the value given to this argument is *5*.
# 3. *to* is a second argument that defines the end of the sequence. In this
#    example this argument takes a value of *20*
# 4. *by* is a third argument that defines the separation between values in
#    the sequence. The argument *by* in this example takes the value of *0.5*


## FUNCTION: a function is an element in R that requests an action from your
## computer. Functions contain algorithms to perform particular tasks ##

## ARGUMENT: an argument specifies or modifies how a function works. Arguments
## are given between parenthesis after a function name ##

### B. SOME "RULES" ABOUT WRITING COMMANDS IN R ##################################

# RULE 1. Each function has its own arguments.

rep(x = 8, times = 3) # *rep* has *x*, *times*
rnorm(n = 10, mean = 2, sd = 1) # *rnorm* has *n*, *mean*, *sd*


# RULE 2. The description of each function and its arguments can be found in the
# help page of the function. To access this information, use the function *help*.

help(topic="seq")


# RULE 3. Arguments USUALLY have names (e.g.: from, to, by), and values are 
# passed to each argument using "=".

seq(from = 5, to = 20, by = 0.5)

# RULE 4. Names of arguments can be eliminated if values are given *in the
# pre-determined order.* For example, these two commands are equivalent:

seq(from = 5, to = 20, by = 0.5) # Most explicit (all with names)
seq(5, 20, 0.5) # Quickest, no argument names

# But they are different from this:
seq(0.5, 5, 20)


# RULE 5. Each function has a pre-determined order for its arguments. 
# For example, for the function *seq* the order is *from* first, then *to*, 
# and then *by*.

help(topic="seq") # See this info in the help page of the function *seq*


# RULE 6. The order of the arguments can be changed ONLY if you use their names. 
# For example, these commands are equivalent:

seq(from = 5, to = 20, by = 0.5)
seq(by = 0.5, from = 5, to = 20)


# RULE 7. Some arguments have pre-determined values!

rnorm(n=10)

# What are the pre-determined values for the argumens in function *rnorm*?
#mean=
#sd=

# These arguments with pre-determined values don't need to be specified for the
# function to work, but one has to be careful. Make sure the defaults are what
# you want!


# RULE 8. When *...* appears in the help file of a function, it frequently means
# multiple arguments with no names. For example, in the function *c*, ... means
# many values that will be concatenated:

help(c) # Also could have been written as help(topic="c")

z<-c(9, 5, 3, 5)

z
# RULE 9. R is case sensitive, so the function *seq* exist, but the 
# function *Seq* does not:
seq(5,20,0.5)
Seq(from = 5, to = 20, by = 0.5)


# RULE 10. In R, white space is meaningless:

seq(from=5, to=20, by=0.5)

seq(from     =     5,      to     =     20,      by     =     0.5)

seq(from=5, 
    to=20, 
    by=
      0.5)



### D. SOME ADDITIONAL EXAMPLES OF SIMPLE COMMANDS (FUNCTIONS AND ARGUMENTS) ###

# *rep*
rep(x = "R", times = 10)

rep(times = 10, x = "R")

rep("R", 10)

rep(10, "R") # Why doesn't this work?

