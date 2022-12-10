## Tutorial: Study Design in Ecology and Epidemiology
## E2M2: Ecological and Epidemiological Modeling in Madagascar
## December 2022

## Cara Brook, 2022

## Premise: You are writing a grant proposal to the NIH to acquire
## funding for a field study aimed at understanding the mechanistic
## drivers that underpin transmission and persistence of henipaviruses
## in wild Eidolon dupreanum fruit bats in Madagascar. In previous work, 
## you fit models to age-structured seroprevalence data to test hypotheses 
## of the mechanisms underlying patterns witnessed in the data but were 
## unable to distinguish between mechanisms of waning immunity and reinfection
## vs. persistent immunity due to sparse data in older age classes.

## Your goal is to decide how many individual bats need to be sampled in 
## each age class to effectively differentiate between the possible mechanisms.


rm(list=ls())
######################################################################

#load functions from previous work