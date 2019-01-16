############################################################################################################
# General worklow:
# 1) Start with pre-cleaned and merged rdt_data_cleaned data file
# 2) Format the data file (include a region_id column and calculated age column)
# 3) Set age bins
# 4) Write a function.age_frequency to determine frequency for each age bin
# 5) Write functions to subset the data:
      # (5.1) function.subsets_rdt_result to pass subsets based on infection status
      # (5.2) function.subsets_site to pass subsets based on site_id
      # (5.3) function.subsets_region to pass subsets based on region
      # (5.4) function.subsets_sex to pass subsets based on sex
      # (5.5) function.nested_subsets to subset by two variables
# 6) Create data tables where the rows are different subsets of data and the columns are age bins
      # (6.1) All regions combined + each region (sex and inf status ignored)
      # (6.2) All regions combined + each region + each site (sex and inf status ignored)
      # (6.3) Table 6.1 + each region by sex (inf status ignored)
      # (6.4) Table 6.1 + each region by inf status (sex ignored)
      # (6.5) Table 6.1 + each region by inf status and by sex 
              # (e.g. rows would be Region 2 overall, Region 2 uninfected men, Region 2 infected men... etc)
# 7) Plot the data as overlayed histograms or as age structure diagrams (women and men mirrored)
############################################################################################################


############################################################################################################
#Set working directory, load necessary packages
############################################################################################################
setwd("/Users/benjaminrice/Dropbox/Lab Projects/Publishing/2018 10 01 Darwin CRS Prevalence Paper/Coding work")
#getwd()
library(tidyverse)
library(reshape2)
############################################################################################################



############################################################################################################
# 1) Start with pre-cleaned and merged rdt_data_cleaned data file
############################################################################################################

#Import cleaned rdt data set; trim out the first column (the counting column that data frames include)
rdt_data_cleaned_untrimmed <- read.csv("rdt_data_cleaned_20181207.csv", stringsAsFactors = FALSE)
rdt_data_cleaned <- rdt_data_cleaned_untrimmed[, 2:7]
############################################################################################################



############################################################################################################
# 2) Format the data file (include a region_id column and calculated age column)
############################################################################################################

#extracting region_id by dividing site_id by 100 and rounding to the nearest integer
age_structure_data <- mutate(rdt_data_cleaned, region_id = as.integer(site_id/100))

#adding an age in years variable by simply substracting the year of birth from 2017 (the year of sampling)
age_structure_data <- mutate(age_structure_data, calc_age = 2017 - dob_yr)

#re-ordering the data file so columns are in a logical order
age_structure_data <- data.frame(
  unique_ind_id = age_structure_data$unique_ind_id,
  region_id = age_structure_data$region_id,
  site_id = age_structure_data$site_id,
  hh_id = age_structure_data$hh_id,
  sex = age_structure_data$sex,
  dob_yr = age_structure_data$dob_yr,
  calc_age = age_structure_data$calc_age,
  rdt_result = age_structure_data$rdt_result
)
############################################################################################################



############################################################################################################
# 3) Set age bins
############################################################################################################

#Create age bin cut-offs. Max age observed is 97 so range will be 0-100.
#2 year increments
cuttoffs_age_bins_2 <- c(seq(0, 100, by=2))

#3 year increments
cuttoffs_age_bins_3 <- c(seq(0, 100, by=3))

#4 year increments
cuttoffs_age_bins_4 <- c(seq(0, 100, by=4))
  
#5 year increments
cuttoffs_age_bins_5 <- c(seq(0, 100, by=5))
  
#10 year increments
cuttoffs_age_bins_10 <- c(seq(0, 100, by=10))
############################################################################################################

  

############################################################################################################
# 4) Write a function.age_frequency to determine frequency for each age bin
############################################################################################################

#Define function.age_frequency with input parameters data and age bin cutoffs

#Psedocode:
#Data input = a list of data frames where each data frame is a subset of data to be analyzed
#Write an outer for loop that loops through the elements of the list (the data subsets being analyzed)
  #Keep the relevant data for the data that needs go in the first column
#Write an inner for loop that loops through the age bins and determines their frequency
#Starting with 0 (cutoff i), count the number of individuals with ages between cuttoff i and cuttoff i + 1
  #divide by n to convert to a frequency
  #Add that bin value as a column to a growing a data frame

#Grouping variables options:
  #all_regions_combined
  #by_region
  #by_site
  #by_sex
  #by_rdt_result
  #by_region_rdt_result
  #by_region_sex
  #by_region_rdt_result_sex

function.age_frequency <- function(data, grouping_variable_input, cuttoffs_age_bins){

  data_list <- data
  grouping_variables <- grouping_variable_input
  age_bin_freq_data_1 <- data.frame()
  column_names <- c(rep(NA, length(cuttoffs_age_bins)))
  row_names1 <- rep(NA, length(data))
  row_names2 <- rep(NA, length(data))
                    
  for(k in 1:length(data_list)){
    
    data_list_k <- data_list[[k]]
    
    if(grouping_variables == "all_regions_combined"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("all_regions_combined")
    }
    if(grouping_variables == "by_region"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("Region", data_list_k[1,2], sep = "_")
    }
    if(grouping_variables == "by_site"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("Site", data_list_k[1,3], sep = "_")
    }
    if(grouping_variables == "by_sex"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("Sex", data_list_k[1,5], sep = "_")
    }
    if(grouping_variables == "by_rdt_result"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("RDT_result", data_list_k[1,8], sep = "_")
    }
    if(grouping_variables == "by_region_rdt_result"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("Region", data_list_k[1,2], "RDT_result", data_list_k[1, 8], sep = "_")
    }
    if(grouping_variables == "by_region_sex"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("Region", data_list_k[1,2], "Sex", data_list_k[1, 5], sep = "_")
    }
    if(grouping_variables == "by_region_rdt_result_sex"){
      row_names1[k] <- paste(grouping_variables)
      row_names2[k] <- paste("Region", data_list_k[1,2], "RDT_result", data_list_k[1, 8],
                             "Sex", data_list_k[1, 5], sep = "_")
    }
    
    for(i in 1:length(cuttoffs_age_bins)){
      column_names[i] <- paste(cuttoffs_age_bins[i],"-", cuttoffs_age_bins[i+1])
      
      bin_i_count <- sum(data_list_k$calc_age > cuttoffs_age_bins[i] & 
                         data_list_k$calc_age <= cuttoffs_age_bins[i+1])
      
      bin_i_freq <- bin_i_count/length(data_list_k$calc_age)
      age_bin_freq_data_1[k, i] <- bin_i_freq
    }
  }
  colnames(age_bin_freq_data_1) <- column_names
  age_bin_freq_data2 <- mutate(age_bin_freq_data_1, grouping = row_names1, subgroup = row_names2)
  age_bin_freq_data3 <- age_bin_freq_data2 %>%
    select(grouping, subgroup, everything())
  age_bin_freq_data <- age_bin_freq_data3[, 1:length(cuttoffs_age_bins)]
  return(age_bin_freq_data)
}
############################################################################################################



############################################################################################################
# 5) Write functions to subset the data:
# (5.1) function.subsets_rdt_result to pass subsets based on infection status
# (5.2) function.subsets_site to pass subsets based on site_id
# (5.3) function.subsets_region to pass subsets based on region
# (5.4) function.subsets_sex to pass subsets based on sex
# (5.5) function.nested_subsets to subset by two variables
############################################################################################################


#Pseudocode:
#Take age_structure_data df and create subset dfs based on variables of interest
#add those subset dfs to a list and output it
#input parameter is a data set

# (5.1) function.subsets_rdt_result to pass subsets based on infection status

function.subsets_rdt_result <- function(data){
  list_subsets_rdt_result = list()
  inf_status_neg <- data[data$rdt_result == 0, ]
  inf_status_pos <- data[data$rdt_result == 1, ]
  list_subsets_rdt_result[[1]] <- inf_status_neg
  list_subsets_rdt_result[[2]] <- inf_status_pos
  return(list_subsets_rdt_result)
}

# (5.2) function.subsets_site to pass subsets based on site_id

function.subsets_site <- function(data){
  list_subsets_site_id = list()
  site_id_v <- unique(data$site_id)
  for(i in 1:length(site_id_v)){
    df_site_i <- data[data$site_id == site_id_v[i], ]
    list_subsets_site_id[[i]] <- df_site_i
  }
  return(list_subsets_site_id)
}

# (5.3) function.subsets_region to pass subsets based on region

function.subsets_region <- function(data){
  list_subsets_region_id = list()
  region_id_v <- unique(data$region_id)
  for(i in 1:length(region_id_v)){
    df_region_i <- data[data$region_id == region_id_v[i], ]
    list_subsets_region_id[[i]] <- df_region_i
  }
  return(list_subsets_region_id)
}

# (5.4) function.subsets_sex to pass subsets based on sex

function.subsets_sex <- function(data){
  list_subsets_sex = list()
  sex_male <- data[data$sex == 'male', ]
  sex_female <- data[data$sex == 'female', ]
  list_subsets_sex[[1]] <- sex_male
  list_subsets_sex[[2]] <- sex_female
  return(list_subsets_sex)
}

# (5.5) function.nested_subsets to subset by two variables

#Inputs:
  #A data frame that contains the data to be analyzed (e.g. age_structure_data)
  #The names of variables to subset by ("rdt_result", "site_id", "region_id", "sex")

#First, loop through the list of variables to extract the ones we are working with ("function choosers")
  #Use those function choosers to run the subsetting/grouping functions written above on the data
  #This makes a list where each element of the list is a data frame that is a group
    #of data of interest ("list_first_grouping")
#Second, take each element of list_first_grouping and pass it to the second grouping/subsetting function
   #(e.g. if first grouped by region, pass the first region to the second grouping function)
  #apply the second grouping function (for example, grouping by RDT result)
  #This gives a second list (e.g. rdt positive, rdt negative) for each element of list_first_grouping
    #The problem is that we create a list (the subgroups) for each element (the group) 
      #of list_first_grouping
    #To solve this problem we can concatenate the subgroups created for each group so that the final output
      #is a properly indexed list containing the subgroups in sequence
        #e.g. region 2 uninfected, region 2 infected, region 3 uninfected, region 3 infected etc
      


function.nested_subsets <- function(data, subsetting_variable_1, subsetting_variable_2){
  
  subsetting_variables_v <- c("rdt_result", "site_id", "region_id", "sex")
  list_first_grouping <- list()
  function1_chooser <- "NA"
  function2_chooser <- "NA"
  holding_list = list()
  nested_subsets_list = list()
  
  for(i in 1:length(subsetting_variables_v)){
    if(subsetting_variable_1 == subsetting_variables_v[i]){
      function1_chooser <- subsetting_variables_v[i]
    }
    if(subsetting_variable_2 == subsetting_variables_v[i]){
      function2_chooser <- subsetting_variables_v[i]
    }
  }
  
  
  if(function1_chooser == 'rdt_result'){
    list_first_grouping <- function.subsets_rdt_result(data)
  }
  if(function1_chooser == 'site_id'){
    list_first_grouping <- function.subsets_site(data)
  }
  if(function1_chooser == 'region_id'){
    list_first_grouping <- function.subsets_region(data)
  }
  if(function1_chooser == 'sex'){
    list_first_grouping <- function.subsets_sex(data)
  }
  
  for(k in 1:length(list_first_grouping)){
    if(function2_chooser == 'rdt_result'){
      list_second_grouping_k <- function.subsets_rdt_result(list_first_grouping[[k]])
    }
    if(function2_chooser == 'sex'){
      list_second_grouping_k <- function.subsets_sex(list_first_grouping[[k]])
    }
    
    holding_list[[k]] <- list_second_grouping_k
  }
  
  if(length(holding_list) == 4){
    nested_subsets_list <- c(holding_list[[1]], holding_list[[2]], holding_list[[3]], holding_list[[4]])
  }
  if(length(holding_list) == 24){
    nested_subsets_list <- c(
      holding_list[[1]],
      holding_list[[2]],
      holding_list[[3]],
      holding_list[[4]],
      holding_list[[5]],
      holding_list[[6]],
      holding_list[[7]],
      holding_list[[8]],
      holding_list[[9]],
      holding_list[[10]],
      holding_list[[11]],
      holding_list[[12]],
      holding_list[[13]],
      holding_list[[14]],
      holding_list[[15]],
      holding_list[[16]],
      holding_list[[17]],
      holding_list[[18]],
      holding_list[[19]],
      holding_list[[20]],
      holding_list[[21]],
      holding_list[[22]],
      holding_list[[23]],
      holding_list[[24]]
      )
  }
  return(nested_subsets_list)
}
############################################################################################################



############################################################################################################
# 6) Create data tables where the rows are different subsets of data and the columns are age bins
############################################################################################################

############################################################################################################
#Creating grouped data lists and calculating age frequencies:
############################################################################################################


#Grouping variables options:
#all_regions_combined
#by_region
#by_site
#by_sex
#by_rdt_result
#by_region_rdt_result
#by_region_sex
#by_region_rdt_result_sex


#All regions combined:
age_structure_data_all_regions_combined <- list(age_structure_data)
age_freqs_all_regions_combined <- function.age_frequency(
  age_structure_data_all_regions_combined, "all_regions_combined", cuttoffs_age_bins_3)

#Grouped by region
age_structure_data_by_region <- function.subsets_region(age_structure_data)
age_freqs_by_region <- function.age_frequency(
  age_structure_data_by_region, "by_region", cuttoffs_age_bins_3)

#Grouped by site
age_structure_data_by_site <- function.subsets_site(age_structure_data)
age_freqs_by_site <- function.age_frequency(
  age_structure_data_by_site, "by_site", cuttoffs_age_bins_3)

#Grouped by sex
age_structure_data_by_sex <- function.subsets_sex(age_structure_data)
age_freqs_by_sex <- function.age_frequency(
  age_structure_data_by_sex, "by_sex", cuttoffs_age_bins_3)

#Grouped by RDT result
age_structure_data_by_rdt_result <- function.subsets_rdt_result(age_structure_data)
age_freqs_by_rdt_result <- function.age_frequency(
  age_structure_data_by_rdt_result, "by_rdt_result", cuttoffs_age_bins_3)

#Grouped by region and RDT result
age_structure_data_by_region_rdt_result <- function.nested_subsets(
  age_structure_data, "region_id", "rdt_result"
)
age_freqs_by_region_rdt_result <- function.age_frequency(
  age_structure_data_by_region_rdt_result, "by_region_rdt_result", cuttoffs_age_bins_3)

#Grouped by region and sex
age_structure_data_by_region_sex <- function.nested_subsets(
  age_structure_data, "region_id", "sex"
  )
age_freqs_by_region_sex <- function.age_frequency(
  age_structure_data_by_region_sex, "by_region_sex", cuttoffs_age_bins_3)

#Grouped by region and RDT_result and sex
list_by_region_rdt_result_sex1 <- function.nested_subsets(
  age_structure_data, "region_id", "rdt_result"
)

function.region_rdt_sex <- function(datalist){
  
  holding_list = list()
  doubly_nested_subsets_list = list()
  
  for(k in 1:length(datalist)){
    list_grouping_k <- function.subsets_sex(datalist[[k]])
    holding_list[[k]] <- list_grouping_k
  }
  
  doubly_nested_subsets_list <- c(
    holding_list[[1]],
    holding_list[[2]],
    holding_list[[3]],
    holding_list[[4]],
    holding_list[[5]],
    holding_list[[6]],
    holding_list[[7]],
    holding_list[[8]]
  )
  
  return(doubly_nested_subsets_list)
}
  
age_structure_data_by_region_rdt_result_sex <- function.region_rdt_sex(list_by_region_rdt_result_sex1)

age_freqs_by_region_rdt_result_sex <- function.age_frequency(
  age_structure_data_by_region_rdt_result_sex, "by_region_rdt_result_sex", cuttoffs_age_bins_3)
############################################################################################################



############################################################################################################
# Creating data tables:
############################################################################################################

# (Table 6.1) All regions combined + each region (sex and inf status ignored)
table_6.1 <- rbind(age_freqs_all_regions_combined, age_freqs_by_region)

# (Table 6.2) All regions combined + each region + each site (sex and inf status ignored)
table_6.2 <- rbind(age_freqs_all_regions_combined, age_freqs_by_region, age_freqs_by_site)

# (Table 6.3) Region by sex (inf status ignored)
table_6.3 <- rbind(age_freqs_all_regions_combined, age_freqs_by_region, age_freqs_by_region_sex)

# (Table 6.4) Regions + each region by inf status (sex ignored)
table_6.4 <- rbind(age_freqs_all_regions_combined, age_freqs_by_region, age_freqs_by_region_rdt_result)

# (Table 6.5) Regions + each region by inf status and by sex 
table_6.5 <- rbind(age_freqs_all_regions_combined, age_freqs_by_region, age_freqs_by_region_rdt_result_sex)

# (Table 6.6) Regions by sex + regions by inf status and by sex
table_6.6 <- rbind(age_freqs_all_regions_combined, age_freqs_by_region_sex, age_freqs_by_region_rdt_result_sex)

#View(table_6.6)
############################################################################################################



############################################################################################################
# 7) Plot the data as overlayed histograms or as age structure diagrams (women and men mirrored)
############################################################################################################

# https://stackoverflow.com/questions/14680075/simpler-population-pyramid-in-ggplot2

#To compare age distribution of infections to age distribution of population for each region

#prep data: 
  #We want a data frame that has:
  #columns:
    #age_bin
    #For each group: sex, freq
      #For example, Region2_gender, Region2_freq

#(multiply frequencies for males by -1 to make them negative to have them going to the left)

#Groups to include:
#Region_2_Male
#Region_2_Female
#Region_2_Male_Positive
#Region_2_Female_Positive
#Region_3_Male
#Region_3_Female
#Region_3_Male_Positive
#Region_3_Female_Positive
#Region_4_Male
#Region_4_Female
#Region_4_Male_Positive
#Region_4_Female_Positive
#Region_5_Male
#Region_5_Female
#Region_5_Male_Positive
#Region_5_Female_Positive


df.plot_age_structure <- data.frame(
  Age_bins = rep(c(seq(3, 96, by=3)), 2),
  Sex = c(rep("Male", length(rep(c(seq(3, 96, by=3))))), rep("Female", length(rep(c(seq(3, 96, by=3)))))),

  Region_2 = c(as.numeric(table_6.6[2, 3:ncol(table_6.6)]), as.numeric(table_6.6[3, 3:ncol(table_6.6)])),
  Region_2_Positive = c(as.numeric(table_6.6[12, 3:ncol(table_6.6)]), as.numeric(table_6.6[13, 3:ncol(table_6.6)])),
  
  Region_3 = c(as.numeric(table_6.6[4, 3:ncol(table_6.6)]), as.numeric(table_6.6[5, 3:ncol(table_6.6)])),
  Region_3_Positive = c(as.numeric(table_6.6[16, 3:ncol(table_6.6)]), as.numeric(table_6.6[17, 3:ncol(table_6.6)])),
  
  Region_4 = c(as.numeric(table_6.6[6, 3:ncol(table_6.6)]), as.numeric(table_6.6[7, 3:ncol(table_6.6)])),
  Region_4_Positive = c(as.numeric(table_6.6[20, 3:ncol(table_6.6)]), as.numeric(table_6.6[21, 3:ncol(table_6.6)])),
  
  Region_5 = c(as.numeric(table_6.6[8, 3:ncol(table_6.6)]), as.numeric(table_6.6[8, 3:ncol(table_6.6)])),
  Region_5_Positive = c(as.numeric(table_6.6[24, 3:ncol(table_6.6)]), as.numeric(table_6.6[25, 3:ncol(table_6.6)]))
)

View(df.plot_age_structure)

#melting the data to prep it for plotting

melted.plot_age_structure <- melt(df.plot_age_structure, id.vars = c("Age_bins", "Sex"))

View(melted.plot_age_structure_region3)

#subsetting by region

melted.plot_age_structure_region2 <- melted.plot_age_structure[
  1:(4 * length(rep(c(seq(3, 96, by=3)))))
, ]

melted.plot_age_structure_region3 <- melted.plot_age_structure[
  (1 + (4 * length(rep(c(seq(3, 96, by=3)))))):(8 * length(rep(c(seq(3, 96, by=3)))))
, ]

melted.plot_age_structure_region4 <- melted.plot_age_structure[
  (1 + (8 * length(rep(c(seq(3, 96, by=3)))))):(12 * length(rep(c(seq(3, 96, by=3)))))
  , ]

melted.plot_age_structure_region5 <- melted.plot_age_structure[
  (1 + (12 * length(rep(c(seq(3, 96, by=3)))))):(16 * length(rep(c(seq(3, 96, by=3)))))
  , ]



#Plotting

p.region2 <- ggplot(data = melted.plot_age_structure_region2,
       aes(x=Age_bins, 
           y=ifelse(test = Sex == "Male", 
                    yes = -value, no = value), 
           fill=variable, 
           color=variable,
           alpha = variable)) +
  geom_bar(stat="identity", position ="identity") +
  coord_flip() +
  scale_colour_manual(values=c("dimgray", "tomato")) +
  scale_fill_manual(values=c("dimgray", "tomato")) +
  scale_alpha_manual(values=c(1, 0.2)) +
  labs(y = "Frequency", x = "Age") +
  scale_y_continuous(labels = abs, limits = 0.25 * c(-1,1))
  
p.region3 <- ggplot(data = melted.plot_age_structure_region3,
                    aes(x=Age_bins, 
                        y=ifelse(test = Sex == "Male", 
                                 yes = -value, no = value), 
                        fill=variable, 
                        color=variable,
                        alpha = variable)) +
  geom_bar(stat="identity", position ="identity") +
  coord_flip() +
  scale_colour_manual(values=c("dimgray", "tomato")) +
  scale_fill_manual(values=c("dimgray", "tomato")) +
  scale_alpha_manual(values=c(1, 0.2)) +
  labs(y = "Frequency", x = "Age") +
  scale_y_continuous(labels = abs, limits = 0.25 * c(-1,1))

p.region4 <- ggplot(data = melted.plot_age_structure_region4,
          aes(x=Age_bins, 
              y=ifelse(test = Sex == "Male", 
                       yes = -value, no = value), 
              fill=variable, 
              color=variable,
              alpha = variable)) +
  geom_bar(stat="identity", position ="identity") +
  coord_flip() +
  scale_colour_manual(values=c("dimgray", "tomato")) +
  scale_fill_manual(values=c("dimgray", "tomato")) +
  scale_alpha_manual(values=c(1, 0.2)) +
  labs(y = "Frequency", x = "Age") +
  scale_y_continuous(labels = abs, limits = 0.25 * c(-1,1))

p.region5 <- ggplot(data = melted.plot_age_structure_region5,
                    aes(x=Age_bins, 
                        y=ifelse(test = Sex == "Male", 
                                 yes = -value, no = value), 
                        fill=variable, 
                        color=variable,
                        alpha = variable)) +
  geom_bar(stat="identity", position ="identity") +
  coord_flip() +
  scale_colour_manual(values=c("dimgray", "tomato")) +
  scale_fill_manual(values=c("dimgray", "tomato")) +
  scale_alpha_manual(values=c(1, 0.2)) +
  labs(y = "Frequency", x = "Age") +
  scale_y_continuous(labels = abs, limits = 0.50 * c(-1,1))
  
#notes: show sample size in legend, left half is males, right half is females



