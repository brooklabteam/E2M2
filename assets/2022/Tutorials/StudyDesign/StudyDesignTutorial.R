## Tutorial: Study Design in Ecology and Epidemiology
## E2M2: Ecological and Epidemiological Modeling in Madagascar
## December 2022

## Cara Brook, 2022

rm(list=ls())

######################################################################
## Part One: Power analysis in Statistical Modeling

## load the pwr package
library(pwr)

######################################################################
## Paired t-test
## Faranky is designing a study to test whether the average mass of several 
## tagged aye-ayes is greater in the wet vs. the dry season. He hypothesizes 
## that aye-ayes might weigh more in the wet season because insects will be 
## more abundant as a result of heavier rains in the wet season. He has already
## captured and found masses of 40 aye-ayes in the dry season, and he is trying 
## to decide whether he needs to capture and weigh more to compare in the 
## wet season. Eventually, he plans to use a paired t-test to assess the
## difference in mass for each aye-aye at different time points. What power 
## would he have to detect a difference at his current sample size?

## effect sizes for t-tests: small ≥ 0.2, medium ≥ 0.5, and large ≥ 0.8.


pwrpt <- pwr.t.test(d=0.2, # small effect size
                    sig.level=0.05,#type I error probability
                    n = 40, # current sample size 
                    type="paired",
                    alternative="greater")
# inspect
pwrpt

## What power does he currently have to detect a difference?

## What would change if his effect size was greater? Or smaller?


######################################################################
## Independent t-test
## Christian is conducting a study to assess whether the average size (length) of
## his favorite Malagasy snake, Langaha madagascariensis, is smaller in Ankarana
## National Park vs. Ranomafana National Park. He measures 175 snakes in Ankarana
## and finds that they appear to be small. He hypothesizes that higher hunting
## pressure in the Ankarana region might result in smaller-sized snakes on average.
## He wonders how many snakes he should measure at Ranomafana to have an 80% 
## probability of demonstrating a significant difference in size if there is one.


# independent t-test
pwr2t <- pwr.t2n.test(d=0.2,
                      n1=175,
                      power = .8,
                      sig.level=0.05,
                      alternative="greater")
# inspect
pwr2t

## How many snakes does he need to measure in Ankarana?

## What would change if the difference in sizes was less dramatic
## than predicted here?

######################################################################
## One Way ANOVA

## Based on his first results, Christian now decides to investigate whether
## Langaha madagascariensis has a different average length across 5 national
## parks in Madagascar. In addition to Ankarana and Ranomafana, he now hopes to
## sample Masoala, Marojejy and Kirindy Mitea national parks. He wonders how many
## snakes he needs to sample in each park to demonstrate these differences

#effect sizes (f) for anovas: small ≥ 0.02, medium ≥ 0.15, and large ≥ 0.35. 


# calculate minimal sample size
pwr.anova.test(k=5,            # number of groups are compared
               f=.25,          # moderate effect size
               sig.level=.05,  # alpha/sig. level = .05
               power=.8)       # confint./power = .8

## How many should he sample per park? Is this number different than above with
## just two populations? If so, why?

######################################################################
## GLM

## ---- Defining Parameters ----

## u is the number of predictors (k) minus 1. The intercepts counts a predictor, too in this case!

## v is defined as the sample size (n) minus the number of predictors (k),
## or n - p. It is the value that is typically associated with the
## degrees of freedom of a model, 

## f2 is the effect size. f2 = r2 / (r - r2)

## power is the probability of finding an effect if it is present, or 
## quantitatively expressed as 1 - P(Type II error). You generally would
## like to have a power value of at least 0.7, which means a 70% change
## of detecting an effect.

## Given the complexity of ecological processes, it is common for linear
## models of ecological systems to have low explanatory power (low r^2), 
## as there are likely many other important factors that are not 
## represented in the the model. It is safe to assume that for any given
## linear model, and r^2 value might be around 0.1 - 0.3. For our examples,
## let's use r2 = 0.1, a conservative estimate of the explanatory power of
## a model. This would make f2 = .1/(1 - .1) = 0.11, or a medium effect size.

## ---- Testing an Example ----

## Tahiry is interested in quantifying predictors of Anopholes occurrence across 
## Madagascar. He is thinking to test a model with the response variable of 
## A. gambiae and the predictors of elevation, temperature, and precipitation.
## He wonders how many sites with disparate elevation/temperature/precipitation
## he will need to sample to find a statistically significant result 70% of the time.

## effect sizes (f2) for GLMs: small ≥ 0.02, medium ≥ 0.15, and large ≥ 0.35

pwrglm <- pwr.f2.test(u = 2,  # number of predictors - 1
                      f2 = .02, # effect size:  f2 = r2 / (r - r2)
                      power=.8, # probability of detecting a significant result if there is one
                      sig.level = 0.05)
# inspect results
pwrglm

## How many sites should he survey?

######################################################################
## Chi-Squared Tests  
## Andres is surveying tourists to Ranomafana National Park, and he wants
## to know if there is variation in self-reported socioeconomic status 
## (low, mid, high income) between Malagasy and vahiny visitors to the park. 
## He hypothesizes that vahiny visitors might tend to be more high income because
## they need to afford an international plane ticket, while Malagasy visitation
## might be more evenly distributed across income classes. He wants to use a
## chi-squared test to address this question, but he is not sure how many 
## tourists in each group (Malagasy vs. vahiny) to survey. He can use the 
## pwr package to help design this study.

## Calculate minimal sample size per group (Malagasy vs. vahiny)


pwrx2 <- pwr.chisq.test(w=0.2, # moderate effect size of socioeconomic status on park visitation
                        df = 2, # degrees of freedom: calculated as: (r-1)(c-1) : (3-1)(2-1)
                        sig.level=0.05, # type I error probability
                        power = .8) # (1 minus type II error probability)

pwrx2
## effect sizes (w) for chi-square: small ≥ 0.1, medium ≥ 0.3, and large ≥ 0.5

## a "moderate" effect size assumes that nationality will interact moderately
## with socioeconomic status to predict park visitation

## What total sample size do you need?
### ???

## What would change if the effect size of nationality on the distribution of 
## visitation by socioeconomic status was smaller? Would you need to sample
## more people per nationality or fewer ?

## What if effect size was even bigger?


######################################################################

## Power analysis for ecological models

## Premise: You are writing a grant proposal to the NIH to acquire
## funding for a field study aimed at understanding the mechanistic
## drivers that underpin transmission and persistence of henipaviruses
## in wild Eidolon dupreanum fruit bats in Madagascar. In previous work, 
## you fit models to age-structured seroprevalence data to test hypotheses 
## of the mechanisms underlying patterns witnessed in the data, but you were 
## unable to distinguish between mechanisms of waning immunity and reinfection
## vs. persistent immunity due to sparse data in older age classes.

## Your goal is to decide how many individual bats need to be sampled in 
## order to effectively differentiate between the possible mechanisms.

## First, load functions from your previous model
source("BatModelFunctions.R")

## Then, load your parameters from your previous model fit
load("bat.mod.par.Rdata")

## Now, simulate the age-seroprevalence in the bat population using your model,
## assuming MSIRN (lifelong immunity).

out.sim <- sim.model(par = par.dat,
                     model = "MSIRN")

# plot output of simulation
p1 <- ggplot(out.sim) + geom_point(aes(x=age, y=seroprevalence, size=N)) +
  geom_line(aes(x=age, y=seroprevalence)) + ylim(c(0,1)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size=18),
                     axis.text = element_text(size=14))                            


print(p1)

## Take a simulated field sample of 1000 bats
sub.df = sub.sample(sim.dat = out.sim, field_sample =  1000)

## plot the data from your field sample
p2 <- ggplot(sub.df) + geom_point(aes(x=age, y=seroprevalence, size=N)) +
  geom_line(aes(x=age, y=seroprevalence)) + ylim(c(0,1)) +
  theme_bw() + theme(panel.grid = element_blank(), axis.title = element_text(size=18),
                     axis.text = element_text(size=14))                            

print(p2)


## How would the data look if you sampled more bats? Fewer bats?
## Try it out.

## Now fit the lifelong immunity hypothesis to your field sampled data.
## You can ignore any warnings returned by R here.
## Be aware that this may take many (5+) minutes to run.

fitted.model.hyp1 <- fit.mod.dat(model="MSIRN",
                                 dat1=sub.df,
                                 par.dat = par.dat,
                                 do.plot=TRUE,
                                 do.save=FALSE, #change to true if you want to save the plot
                                 filename = NA) #if saving the plot, provide a file name

## Based on the plot, how well does the model compare to the data?
## How would this change if you sub-sampled more or less data?

## What AIC is returned? (look at your output)
fitted.model.hyp1

## Now, test the alternative hypothesis of waning immunity
## Beware: This might take even longer! (why?)

fitted.model.hyp2 <- fit.mod.dat(model="MSIRS",
                                 dat1=sub.df,
                                 par.dat = par.dat,
                                 do.plot=TRUE,
                                 do.save=FALSE, #change to true if you want to save the plot
                                 filename = NA) #if saving the plot, provide a file name

## How well does this fit compare to the data? (look at output)
fitted.model.hyp2 

## Is the sample size large enough to differentiate the models?

## What would change if the real dynamics were MSIRS instead of MSIRN? Can you test
## this with your existing tools?
