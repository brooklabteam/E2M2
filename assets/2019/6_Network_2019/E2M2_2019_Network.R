
# Public health officials are concerned about an outbreak of an unknown parasite in a community near the great republic of E2M2
# they want to evaluate the risks of disease spreading in a small community
# and have collected interactions data among a group of 79 scientists from 3 institutions

#each person recorded whether they have spent some time with each other during the past 3 days.

#import dataframe

network<-read.csv("edgelist.csv")
head(network)
summary(network)

library(reshape2)
network<-acast(network,Node1~Node2, fun.aggregate= mean,value.var="Link")
diag(network)<-0
network<-as.matrix(network)

###load package igraph that we installed previously
library(igraph)

#The public health official also recorded other information about each citizen
#ID, Sex, Institution and whether the person had a superpower
#import nodeinfo object
nodeinfo<-read.csv("nodeinfo_all.csv",header = T)
summary(nodeinfo)
str(nodeinfo)

#using the graph.adjacency function we can build a very simple network
net=graph.adjacency(network,mode="undirected",diag=FALSE,weighted=TRUE)
plot(net)

###########################################################################################################

#the next bit of code is to help make a better representation of the information on the network
#for example, we want to assign a color for each institution (student, mentor, instructor), and color for institution 
#match Institution (from the nodeinfo dataframe to each ID from the matrix)
V(net)$Sex=as.character(nodeinfo$Sex[match(V(net)$name,nodeinfo$CitizenID)])
V(net)$shape=V(net)$Sex #assign the "Institution" attribute as the vertex shape
V(net)$shape=gsub("F","circle",V(net)$shape) #Females will be circle
V(net)$shape=gsub("M","square",V(net)$shape) #Males will be square
plot(net, vertex.label=NA, vertex.size=2)

V(net)$Institution=as.character(nodeinfo$Institution[match(V(net)$name,nodeinfo$CitizenID)])
V(net)$color=V(net)$Institution
V(net)$color=gsub("student","red",V(net)$color) #positives will be red
V(net)$color=gsub("mentor","yellow",V(net)$color) #negatives will be black
V(net)$color=gsub("instructor","blue",V(net)$color)

fr <- layout_nicely(net, dim=2) #uses a specific algorithm to draw the network, many posiibilities, *fr*, *kk*)
plot(net, 
     vertex.label=NA,
     vertex.size=3,
     layout=fr, 
     edge.lty=1,edge.width=1,
     margin=c(-0.1,-0.1,0,0))
legend(x=0, y=-0.9, c("Male","Female"), 
       pch=c(0,1), pt.cex=1.5, cex=.8, bty="n", ncol=1,y.intersp = 1.2)
legend(x=-0.9, y=-0.9, c("student", "mentor", "instructor"), 
       pch=16, col=c("red","yellow","blue"), bg=c("red","yellow","blue"), pt.cex=1.5, cex=.8, bty="n", ncol=1, y.intersp= 1.2)
#calculate degree and betweenness for each node
deg<-degree(net)
bet<-betweenness(net)
nodeinfo$deg<-degree(net)
nodeinfo$bet<-betweenness(net, directed=F)
summary(bet)
plot(net, 
     vertex.size=deg/5,
     vertex.label.cex=0.5,
     vertex.label.dist=0,
     layout=fr, 
     edge.lty=1,edge.width=1,
     margin=c(-0.1,-0.1,0,0))

legend(x=0, y=-0.9, c("Male","Female"), 
       pch=c(0,1), pt.cex=1.5, cex=.8, bty="n", ncol=1,y.intersp = 1.2)
legend(x=-0.9, y=-0.9, c("student", "mentor", "instructor"), 
       pch=16, col=c("red","yellow","blue"), bg=c("red","yellow","blue"), pt.cex=1.5, cex=.8, bty="n", ncol=1, y.intersp= 1.2)

######## what is the average number of steps between any two pairs of individuals?
######## what is the proportion of realized links /possible links
mean_distance(net, directed=FALSE)
edge_density(net, loops=FALSE)

########transitivity measures the probability that the adjacent vertices of a focal vertex are connected aka clustering coefficient
transitivity(net) #are my friends also friend with each other?

#We are done with descriptions of the network and now we can answer the questions of our beloved public health officials
#################################################Dyadic, pairwise relationships among individuals
####Are individuals of the same institution group, sex more likely to interact with each other, than predicted by chance?

#######################              MULTIPLE REGRESSION QUADRATIC ASSIGNMENT PROCEDURES-model
library(asnipe)
#multiple regression quadratic assignement procedure
#Calculate the regression coefficient for each input matrix. 
 
#This tests whether y(interaction) is related to
#x1 (sex similarity) while controlling for x2 (group belonging) and vice versa.
#Load the overall network create a Group or sex similarity matrix define group memberships
#####
sx<-data.frame(nodeinfo$CitizenID,nodeinfo$Sex)
gbs<-get_group_by_individual(sx,data_format = "individuals")
net_sx<-get_network(gbs,data_format = "GBI")

instit<-data.frame(nodeinfo$CitizenID,nodeinfo$Institution)
gbinstit<-get_group_by_individual(instit,data_format = "individuals")
net_inst<-get_network(gbinstit,data_format = "GBI")

sp<-data.frame(nodeinfo$CitizenID,nodeinfo$superpower)
gbsp<-get_group_by_individual(sp,data_format = "individuals")
net_sp<-get_network(gbsp,data_format = "GBI")

#mrqap to explain similarity matrix based on distance matrix, similarity in species and sex
m.all<-mrqap.dsp(network~net_sp+net_inst+net_sx,intercept=TRUE,directed="undirected",diagonal=F,test.statistic="beta",randomisations=1000)
m.all

#Here, we will use the built in sir function to make simulations of a SIR model on the network
simulation<-sir(net, 1,1 , no.sim = 1000)
par(mfrow=c(1,1))
plot(simulation, comp="NS")
plot(simulation, comp="NR")
plot(simulation, comp="NI")

