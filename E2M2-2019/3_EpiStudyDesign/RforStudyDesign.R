
######################################################################################
## Randomly selecting study sites
######################################################################################

# Import a table of study sites
study_site_table <- read.csv("study_sites.csv", stringsAsFactors = FALSE)

# Load the tidyverse package to access some useful functions
library(tidyverse)

# Randommly select 6 sites within each region
## use the group_by() and random() function
selected_sites_table <- group_by(study_site_table, region_id) %>%
  mutate(site_random_number = sample(1:length(region_id), length(region_id))) %>%
  mutate(site_selected_y_n = 
           ifelse(site_random_number <= 6, 1, 0))

selected_sites_table <- selected_sites_table[
  selected_sites_table$site_selected_y_n == 1, ]
           
######################################################################################
## Monitoring data entry
######################################################################################

# Concatenate data fields in a row into one string using a function and a for loop
concatenate.function <- function(dataframe){
  df <- dataframe
  concatenated_df <- rep(0, length(df1[,1]))
  for(i in 1:length(df1[,1])){
    concatenated_df[i] <- paste(
      df[i, 1], df[i, 2], df[i, 3], df[i, 4], df[i, 5], sep = ".")
  }
  return(concatenated_df)
}

df1.concatenated <- concatenate.function(df1)
df2.concatenated <- concatenate.function(df2)

# compare the repeated data entries to see if entries are the same
compare.function <- function(dataframe1, dataframe2){
  rows_with_errors <- rep(NA, length(dataframe1))
  for(i in 1:length(dataframe1)){
    rows_with_errors[i] <- ifelse(dataframe1[i] == dataframe2[i], 0, 1)
  }
  return(rows_with_errors)
}

compare.function(df1.concatenated, df2.concatenated)


######################################################################################
## Merging data sheets
######################################################################################

# Import clinical data and survey data
clinical_data <- read.csv("clinical_data.csv", stringsAsFactors = FALSE)
survey_data <- read.csv("survey_data.csv", stringsAsFactors = FALSE)

# Subset the clinical data to keep the relevant variables
rdt_subset_variables <- c("unique_ind_id", "rdt_result")
clinical_data_trimmed <-clinical_data[rdt_subset_variables]

# Merge data frames using the full_join() function

# First, need to check for duplicate IDs using anyDuplicated()
# If the output of anyDuplicated > 0, then there are duplicate IDs
anyDuplicated(clinical_data_trimmed$unique_ind_id)
anyDuplicated(survey_data$unique_ind_id)

# Merge dataframes using the full_join() function and store as data_joined
data_joined <- full_join(clinical_data_trimmed, survey_data, by = "unique_ind_id")

# Check merged sheet
head(data_joined)
str(data_joined)


######################################################################################
## Identifying errors
######################################################################################

#Data cleaning steps (for trying to clean a data file named "raw_rdt_data_joined")
#(3.1) Use head() and str() to get a look at the data and assess cleaning needs
#(3.2) Trim out the individuals that do not have an RDT result and the individuals with an RDT result 
#but no age or sex data
#(3.3) Create new columns where you decompose unique_ind_id into a site id ("site_id) 
#and a household_id ("hh_id")
#(3.4) Recode rdt_result: 
    #trim out invalid RDT result individuals (2) n = N = negative = 0; pan or panpf = 1
#(3.5) Peek at summary of dob_yr to ensure that no absurd values (less than zero or over 100) are included
#(3.6) If data looks good, store as "rdt_data_cleaned"
#---------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------#
#(3.1) Use head() and str() to get a look at the data and assess cleaning needs
head(raw_rdt_data_joined)
str(raw_rdt_data_joined)
#---------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------#
#(3.2) Trim out the individuals that do not have an RDT result and the individuals with an RDT result 
#but no age or sex data

#Using the is.na function and subsetting to trim individuals who have NA for:
#(1) RDT result (2) sex (3) dob_yr
#Storing the output as raw_rdt_data
raw_rdt_data <- raw_rdt_data_joined[!is.na(raw_rdt_data_joined$rdt_result), ]
raw_rdt_data <- raw_rdt_data[!is.na(raw_rdt_data$sex), ]
raw_rdt_data <- raw_rdt_data[!is.na(raw_rdt_data$dob_yr), ]

#Use str() to check the output
str(raw_rdt_data)

#---------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------#
#(3.3) Create new columns to decompose unique_ind_id into a site id ("site_id) and a household_id ("hh_id")

#use pipes %>% to pass raw_rdt_data to mutate()
#use as.integer() and dividing by 1000 to trim the numeric unique_ind_id to the site_id and hh_id digits

raw_rdt_data <- raw_rdt_data %>% mutate(
  site_id = as.integer(raw_rdt_data$unique_ind_id/10000000),
  hh_id = as.integer(raw_rdt_data$unique_ind_id/1000)
)
#---------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------#
#(3.4) Recode rdt_result: (1) trim out invalid RDT result individuals (2) n = N = negative = 0; pan or panpf = 1

#use subsetting to keep only the individuals with a valid RDT result
raw_rdt_data <- raw_rdt_data[!raw_rdt_data$rdt_result == 'i',]

#Check: use a simple barplot to see if "i"s have been removed
plot_counts1 <- table(raw_rdt_data$rdt_result)
barplot(plot_counts1)

#recode the rdt_results such that they are 1s (positive) or 0s (negatives)
raw_rdt_data$rdt_result <- recode(raw_rdt_data$rdt_result, N = 0, n = 0, pan = 1, panpf = 1)

#Check: use a simple barplot to see if the 1, 0 recoding was correct
plot_counts2 <- table(raw_rdt_data$rdt_result)
barplot(plot_counts2)

#from str() we see that rdt results are stored as num even though just 1s and 0s.Convert to integer:
raw_rdt_data$rdt_result <- as.integer(raw_rdt_data$rdt_result)

str(raw_rdt_data)

#check with a histogram
hist(raw_rdt_data$rdt_result)
#---------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------#
#(3.5) Peek at summary of dob_yr to ensure that no absurd values (less than zero or over 100) are included:
head(raw_rdt_data)
hist(raw_rdt_data$dob_yr)
summary(raw_rdt_data$dob_yr)
#---------------------------------------------------------------------------------------------------------#


#---------------------------------------------------------------------------------------------------------#
#(3.6) Data looks good so store as "rdt_data_cleaned"
rdt_data_cleaned <- raw_rdt_data

#If needed, export raw_rdt_data_cleaned as a CSV
#rite.csv(rdt_data_cleaned, file = "rdt_data_cleaned_20181207")


######################################################################################
## Create summary tables: Prevalence by site
######################################################################################

#(4.1) prev_by_site_table : Prevalence, number RDT positive and n per site

#Data table columns: site_id, n_site, n_rdt_pos_site, prev_site
#site_id: use unique() on site_id
#n_site: group_by site_id and length() on rdt_result
#n_rdt_pos_site: group_by site_id and sum() on rdt_result
#prev_site: mutate() to divide n_site by n_rdt_pos_site and multiply by 100
#n_hh_site: mutate()

prev_by_site_table <- group_by(rdt_data_cleaned, site_id) %>%
  mutate(n_site = length(rdt_result)) %>%
  mutate(n_rdt_pos_site = sum(rdt_result)) %>%
  mutate(prev_site = (n_rdt_pos_site/n_site) * 100)  %>%
  mutate(n_hh_site = length(unique(hh_id)))

prev_by_site_table <- select(
  prev_by_site_table[!duplicated(prev_by_site_table$site_id), ],
  site_id, n_site, n_rdt_pos_site, prev_site, n_hh_site
)

#check the prevalence by site table
View(prev_by_site_table)
plot(prev_by_site_table$site_id, prev_by_site_table$prev_site)
prev_by_site_table

######################################################################################
## Generate figures: Age structure of malaria infections
######################################################################################

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


######################################################################################
## Writing
######################################################################################

# R markdown: https://rmarkdown.rstudio.com/


