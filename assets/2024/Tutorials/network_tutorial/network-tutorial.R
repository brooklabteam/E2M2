# E2M2 Network Tutorial 2022
# Project: 
# Authors: 

# Michelle V Evans
# Github: mvevans89
# Institut de Recherche pour le Developpement
# Email: evans.michelle@ird.fr

# Script originated November 2022

# Description of script and Instructions #####################

#' This script contains the tutorial on using social networks in
#' epidemiology that was provided as part of the E2M2 workshop 
#' in December 2022.
#' 

# Packages & Options #########################################

options(stringsAsFactors = F, scipen = 999, digits = 4)

#networks
library(igraph)
library(asnipe)

library(tidyr)
library(dplyr)

# Loading the datasets ###################################


#' Public health officials are concerned about an outbreak of an unknown parasite in a community near the great republic of E2M2
#' they want to evaluate the risks of disease spreading in a small community
#' and have collected interactions data among a group of 79 scientists from 3 institutions. 
#' Each person recorded whether they have spent some time with each other during the past 3 days.

# import the dataset
edgelist <- read.csv("edgelist-e2m2.csv")
head(edgelist)
str(edgelist)

# reformat into a matrix
network <- pivot_wider(edgelist, names_from = Node2, values_from = Link)
#move name of node1 to the rowname and turn into a matrix
network.mat <- select(network, -Node1) %>%
  as.matrix()
rownames(network.mat) <- network$Node1
#set diagonal (contact with oneself) to 0
diag(network.mat) <- 0

#investigate small portion of this matrix to see what it looks like
network.mat[1:6,1:6]

#load the metadata
nodeinfo <- read.csv("nodeinfo_all.csv",header = T)
str(nodeinfo)
head(nodeinfo)

# Visualize the network ###########################

#simple plot
net <- graph.adjacency(network.mat, mode = "undirected", diag = FALSE, weighted = TRUE)
plot(net)

#use metadata to change node attributes
#first check that the order of nodes in the nodeinfo matches that in the network
nodeinfo$CitizenID == V(net)$name #all should be TRUE

#if it is not true, then you need to reorder the nodeinfo based on order of names in network
nodeinfo <- nodeinfo[match(V(net)$name,nodeinfo$CitizenID),]

# shape is dependent on sex
V(net)$shape <- ifelse(nodeinfo$Sex == "F", "circle", "square")

#color is dependent on institution
#create a named color palette
color.pal  <- c("instructor" = "blue", "mentor" = "yellow", "student" = "red")

#assign color by matching institution to this palette
V(net)$color = color.pal[match(nodeinfo$Institution, names(color.pal))]

#plot with legend
fr <- layout_nicely(net, dim=2) #uses a specific algorithm to draw the network
plot(net, 
     vertex.label=NA,
     vertex.size=4,
     layout=fr, 
     edge.lty=1,edge.width=1,
     margin=c(-0.1,-0.1,0,0))
legend(x=0, y=-0.9, c("Male","Female"), 
       pch=c(0,1), pt.cex=1.5, cex=.8, bty="n", ncol=1,y.intersp = 1.2)
legend(x=-0.9, y=-0.9, c("student", "mentor", "instructor"), 
       pch=16, col=c("red","yellow","blue"), bg=c("red","yellow","blue"), 
       pt.cex=1.5, cex=.8, bty="n", ncol=1, y.intersp= 1.2)

# Calculate network characteristics ##################

## Node-level characteristics ########################

#degree
nodeinfo$deg <- degree(net)
#betweeness
nodeinfo$bet <- betweenness(net, directed=F)
head(nodeinfo)

# add the degree as a characteristic of the node to the network to visualize
V(net)$size <- nodeinfo$deg/2
plot(net, 
     vertex.label.cex=0.5,
     vertex.label.dist=0,
     layout=fr, 
     edge.lty=1, edge.width=1,
     margin=c(-0.1,-0.1,0,0))

legend(x=0, y=-0.9, c("Male","Female"), 
       pch=c(0,1), pt.cex=1.5, cex=.8, bty="n", ncol=1,y.intersp = 1.2)
legend(x=-0.9, y=-0.9, c("student", "mentor", "instructor"), 
       pch=16, col=c("red","yellow","blue"), bg=c("red","yellow","blue"), 
       pt.cex=1.5, cex=.8, bty="n", ncol=1, y.intersp= 1.2)

## Network-level characteristics ######################

# Mean Distance: What is the average number of steps between any two pairs of individuals?
mean_distance(net, directed=FALSE)

# Edge Density: What is the proportion of realized links compared to all possible links?
edge_density(net)

# Transitivity: What is the probability that the adjacent vertices of a focal vertex are also connected?
transitivity(net)

# Multiple Regression Quadratic Assignment Procedures (MRQAP) ############

#' The public health authorities have the following question: Are individuals of
#' the same institution or sex more likely to interact with each other than predicted by chance?

## Create dependent variable matrices ####################################

# First, create matrices to describe individual's membership in sex and institution
# Group by Sex
sex.info <- select(nodeinfo, ID = CitizenID, Group = Sex)
gbs <- get_group_by_individual(sex.info, data_format = "individuals")
net.sx <- get_network(gbs, data_format = "GBI")

#inspect how these compare
sex.info[1:6,]
net.sx[1:6,1:6]

# group by institution
inst.info <- select(nodeinfo, ID = CitizenID, Group = Institution)
gbs <- get_group_by_individual(inst.info, data_format = "individuals")
net.inst <- get_network(gbs, data_format = "GBI")

## Statistical Model ###############################

net.model <- mrqap.dsp(network.mat ~ net.sx + net.inst,
                       intercept = TRUE, 
                       directed = "undirected",
                       diagonal = F, test.statistic = "beta", 
                       randomisations = 99)

net.model

# SIR Simulations on Networks ######################

# you can play around with changing beta and gamma
net.sims <- sir(net, beta = 0.2, gamma = 0.2 , no.sim = 1000)

#number of susceptible
plot(net.sims, comp = "NS")
#number of infected
plot(net.sims, comp = "NI")
#number of recovered
plot(net.sims, comp = "NR")
