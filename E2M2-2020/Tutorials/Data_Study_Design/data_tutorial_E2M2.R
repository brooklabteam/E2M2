#clear working directory
rm(list = ls())

#set here the directory where your data files are located
setwd("/Users/caraebrook/Documents/R/R_repositories/E2M2/E2M2-2020/Tutorials/Data_Study_Design/")

#load library
library(dplyr)

######## MERGING DATASHEETS ########
#Ekipa Fanihy just returned to Antananarivo from two separate field missions
#Sarah and Angelo sampled 53 bats in Ankarana
#Cara's team sampled 49 bats at two sites near Moramanga
#Sarah and Cara added their data to the compiled project data on the same day
#so now there are two overlapping but incomplete project datasheets that need to be combined
#let's take a look at these two versions of the project datasheet
dat.sarah <- read.csv("Flying_Fox_Capture_Data_12_06_2019_SG.csv", skip=2, header = TRUE)
dat.cara <- read.csv("Flying_Fox_Capture_Data_12_06_2019_CB.csv", skip=2, header = TRUE)

#it's easy to quickly calculate the dimensions of each dataframes
#the number of entries don't match up--their two versions are incomplete and need to be combined
dim(dat.sarah)
dim(dat.cara)

#we can try joining their datasheets with a simple rbind
cap.dat <- rbind(dat.sarah, dat.cara)

#but then the dataframe is suspiciously large...
dim(cap.dat)

#there's a much smaller number of unique entries
unique(cap.dat$Sample.ID)

#we can filter duplicates using the dyplr package
cap.dat <- rbind(dat.sarah, dat.cara) %>% distinct(Sample.ID, .keep_all = T)

#the dimensions match the number of unique entries above!
#now we've combined Sarah and Cara's datasheets into single compiled project datasheet
dim(cap.dat)

######## STOCKING SAMPLES ########
#Ekipa Fanihy stored all of the samples from both missions in the -80C freezer at IPM
#they separate sample types into different boxes and put the samples in numeric order
#they track which sampleIDs are in each box in a large excel spreadsheet
#but it's still easy for samples to get lost
#sometimes samples never make it to the freezer
#or data may have been entered incorrectly
#this time, Sarah thinks that some of her Ankarana samples may have fallen out of their liquid nitrogen canisters 
#and are now floating loosely in the liquid nitrogen tank
#we can identify which samples are still stuck in the liquid nitrogen 
#by checking our freezer stocking datasheet against our field datasheet

#read in freezer stocking data
freeze.dat <- read.csv("data_tutorial/FreezerStock_12_13_2019_SG.csv", header = TRUE)

#remember that we store all of the project data in one large datasheet 
#so let's subset Ankarana data from the most recent mission
Ankarana.dat <- cap.dat[c(cap.dat$Roost.Site=="Ankarana_Andrafiabe" | cap.dat$Roost.Site=="Ankarana_Chauves_Souris"), ]  
Ankarana.dat$Date <- as.Date(Ankarana.dat$Date, format = "%m/%d/%y")
Ankarana <- subset(Ankarana.dat, Date >= "2019-11-27" & Date <= "2019-11-30")

#in order to compare the field and freezer data, we need to find a common variable
names(freeze.dat)  
names(cap.dat) 

#both datasheets have a column for the sample ID, but they don't quite match up
#so let's change the freezer sample ID column name to match the field data
colnames(freeze.dat)[5] <- "Sample.ID"

#we don't always collect every sample type from each bat
#which means even though we caught 53 bats, we likely have less of some sample types
#let's start by checking whether all of our urine samples made it to the freezer
#first, we need to come up with a list of all the bat sampleIDs for which we have urine samples
#how can we subset sample IDs for which we collected a urine sample?
#let's take a look at how we stored info about our urine samples--what values did we record?
unique(Ankarana.dat$Urine)

#looks like for each bat, we recorded the number of urine samples we collected from that individual
#our values range from 0-2
#let's create a all individuals that have at least 1 urine sample 
UR <- Ankarana[!Ankarana$Urine==0,]

#now we can use a filter to find any sample IDs in our list of urine samples 
#that are NOT in the list of freezer samples that have a sample type "Urine"  
UR.missing = UR%>%
  filter(!(Sample.ID %in% freeze.dat$Sample.ID[freeze.dat$Sample.Type=="Urine"]))

#looks like 2 urine samples are missing 
#we can make a note that WAY081 and WAY083 are likely still floating in our liquid nitrogen
UR.missing$Sample.ID

#let's do the same check for our sample types
#but it's annoying to have to keep writing the code above..
#can we try making a function?
fun <- function(samples, type){
  missing = samples%>%
    filter(!(Sample.ID %in% freeze.dat$Sample.ID[freeze.dat$Sample.Type==type]))
  print(missing$Sample.ID)
}

#let's test out our function with the urine samples to make sure it we get the same result
UR.missing <- fun(samples = Ankarana[!Ankarana$Urine==0,], type = "Urine")

#looks like all of our throat samples made it to the freezer
THR.missing <- fun(samples = Ankarana[!Ankarana$Throat==0,], type = "Throat")

#same with our fecal samples
FEC.missing <- fun(Ankarana[!Ankarana$Fecal...UTM==0,], type = "Feces")

#but we're missing a lot of our wing punch samples...
WP.missing <- fun(Ankarana[!Ankarana$Wing.Punches....==0,], type = "Wing Punch")

#I'd be surprised if all of those samples floated into our liquid nitrogen
#let's see if there's something else going on...
#our function pulls sample IDs based on two values: 
#1) the number of wing punches (value in our field data)
#2) the sample type recorded in the freezer data
#first let's check to see whether we have weird values for the number of wing punches
unique(Ankarana.dat$Wing.Punches....) #nope, looks like our only value is "1"

#let's also check the sample type variable in the freezer data
unique(freeze.dat$Sample.Type)

#Looks like Sarah forgot to capitalize consistently! There's 'Wing Punch' and 'Wing punch'
#our function wouldn't recognize the 'Wing punch' samples
#we need to capitalize 'punch'
#we can do using the global substitute function to replace "Wing punch" with "Wing Punch"
freeze.dat$Sample.Type <- gsub("Wing punch", "Wing Punch", freeze.dat$Sample.Type)

#let's check for missing wing punch samples again...
#looks like they're all there now!
WP.missing <- fun(Ankarana[!Ankarana$Wing.Punches....==0,], type = "Wing Punch")

#Typically, we would go through all of our sample types
#but we won't bore you in this tutorial...let's move onto exports

######## EXPORTS ########
#we want to export all of our wing punch samples to UC Berkeley for sequencing
#but in order to prepare the export documents, 
#we need to know exactly how wing punch samples are in our freezer
#and what boxes we need to find
#this should be easy with our freezer stocking datasheet!
#let's first find the boxes
wing.boxes <- freeze.dat$Box.Number[freeze.dat$Sample.Type=="Wing Punch"]
unique(wing.boxes) #we'll need to find these 15 boxes in our freezer...good thing we have a freezer map!

#now let's figure out how many samples we have for each bat species (Eidolon dupreanum, Pteropus rufus, Rousettus madagscariensis)
#we'll have to report these numbers for the export permits
wing.ED <- freeze.dat[freeze.dat$Sample.Type=="Wing Punch" & freeze.dat$Species=="Eidolon dupreanum",] 
length(wing.ED$Sample.ID)

wing.RM <- freeze.dat[freeze.dat$Sample.Type=="Wing Punch" & freeze.dat$Species=="Rousettus madagascariensis",] 
length(wing.RM$Sample.ID)

wing.PR <- freeze.dat[freeze.dat$Sample.Type=="Wing Punch" & freeze.dat$Species=="Pteropus rufus",] 
length(wing.PR$Sample.ID)

