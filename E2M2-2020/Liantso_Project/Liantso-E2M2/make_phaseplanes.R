rm(list =ls())


library(plyr)
library(dplyr)
library(ggplot2)

#load data
ang.dat <- read.csv(file="life_history_dat_angonoka.csv", header = T,stringsAsFactors = F)
head(ang.dat)

isocline.survival.wrap <- function(df, p_sequence){
  #new version essentially plots isocline wrap based on the age at first rep for each species
  
  zero_growth_line.list <- list()
  for(i in 1:2){
    zero_growth_line.list[[i]] <- lapply(p_sequence, get_isoclines_survival, IBI_fecundity=as.numeric(as.character(df$IBI_fecundity[i])), IBI=as.numeric(as.character(df$IBI[i])), age_1st_rep=as.numeric(as.character(df$age_1st_rep[i])))           
  }
  zero_growth_line <- matrix(unlist(zero_growth_line.list[[1]]), ncol=2, byrow=T)
  zero_growth_line <- data.frame(zero_growth_line)
  names(zero_growth_line) <- c("IBI_infant_survival_rate", "IBI_adult_survival_rate")
  zero_growth_line$scenario <- "best"
  zero_growth_line_worst <- matrix(unlist(zero_growth_line.list[[2]]), ncol=2, byrow=T)
  zero_growth_line_worst <- data.frame(zero_growth_line_worst)  
  names(zero_growth_line_worst) <- c("IBI_infant_survival_rate", "IBI_adult_survival_rate")
  zero_growth_line_worst$scenario <- "worst"
  
  #bind to other
  zero_growth_line <- rbind(zero_growth_line, zero_growth_line_worst)
  
  #now attach some identifier information to this
  zero_growth_line$species <- unique(df$species)
  zero_growth_line$IBI[zero_growth_line$scenario=="best"] <- df$IBI[df$scenario=="best"]
  zero_growth_line$IBI[zero_growth_line$scenario=="worst"] <- df$IBI[df$scenario=="worst"]
  
  zero_growth_line$IBI_fecundity[zero_growth_line$scenario=="best"] <- df$IBI_fecundity[df$scenario=="best"]
  zero_growth_line$IBI_fecundity[zero_growth_line$scenario=="worst"] <- df$IBI_fecundity[df$scenario=="worst"]
  
  zero_growth_line$age_1st_rep[zero_growth_line$scenario=="best"] <- df$age_1st_rep[df$scenario=="best"]
  zero_growth_line$age_1st_rep[zero_growth_line$scenario=="worst"] <- df$age_1st_rep[df$scenario=="worst"]
  
  #zero_growth_line$avg_lifespan[zero_growth_line$scenario=="best"] <- df$avg_lifespan[df$scenario=="best"]
  #zero_growth_line$avg_lifespan[zero_growth_line$scenario=="worst"] <- df$avg_lifespan[df$scenario=="worst"]
  
  zero_growth_line$litter_size[zero_growth_line$scenario=="best"] <- df$litter_size[df$scenario=="best"]
  zero_growth_line$litter_size[zero_growth_line$scenario=="worst"] <- df$litter_size[df$scenario=="worst"]
  
  #and for good measure, include the annual rate too
  #zero_growth_line$lit_infant_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="infant"])/IBI
  #zero_growth_line$lit_juvenile_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="juvenile"])/IBI
  #zero_growth_line$lit_adult_annual_mortality_rate <- (df$stage_mortality_rate[df$stage=="adult"])/IBI
  
  
  return(zero_growth_line)
}
get_isoclines_survival <- function(p, IBI_fecundity, IBI, age_1st_rep){
  
  #first, choose 
  #first record timesteps of adulthood
  #dur_rep <- ((avg_lifespan/IBI)-(age_1st_rep/IBI)) #duration of reproduction in IBI timesteps
  #dur_mat <- ((avg_lifespan/IBI)-1) # durating of juvenile or adult class in IBI timesteps
  
  #and here is where we actually choose our proper IBI fecundity
  
  #and get all our variables in order
  #s was fed in through the function
  #then, adult stage reproductive output (must be multiplied by the probability of escaping this class)
  
  Fa = IBI_fecundity
  
  #our system somehwat breaks down in the case where age at 1st rep = 2 because the first birth cohort must st
  #*(dur_rep)
  #then, age at first rep - this has been carefully chosen to be an even number divisible by the IBI
  a = age_1st_rep/IBI
  
  #L = avg_lifespan/IBI
  
  #then, we solve for i:
  
  #i <- (((-p^(L))+sqrt(p^(2*L)+4*Fa*p^(a+1)))/(2*(Fa*p^(a+1))))
  i <- (1-p)/(Fa*(p^a))
  
  #however, if age at 1st rep=1, then we have no juvenile class and we have to make amends
  #in this case our equation is given by:
  #maybe not. I think our equation still holds. but the leslie matrix will change for sure
  #if(a==1){
  # i <- -(1/Fa)
  #}
  
  
  #adult_surv =  p #IBI adult mortality rate
  #infant_surv = (p*i)
  
  adult_surv =  p
  infant_surv = (p*i)
  #and then return them!
  #return(c(infant_surv,adult_surv))
  return(c(infant_surv,adult_surv))
  #return(c(i,p)) #if you would rather plot annual mortaluty...
}
wrap.isocline.survival <- function(dat1, type, p_sequence, do.plot, do.save, filename){
  
  
  dat1$adult_IBI_mortality <- 1- dat1$adult_IBI_survival 
  dat1$infant_IBI_mortality <- 1- dat1$infant_IBI_survival 
  dat1$infant_IBI_mortality[dat1$infant_IBI_mortality<0 & dat1$scenario=="best"] <- 0
  
  
  #now split to list
  dat.list <- dlply(dat1, .(species))
  
  #now make list of isoclines per species
  #start with
  isocline.list <- lapply(dat.list, isocline.survival.wrap, p_sequence=p_sequence)
  #and bind them all together
  isocline.df <- do.call("rbind", isocline.list)
  rownames(isocline.df) <- c()
  
  #add a grouping ID for each line
  isocline.df$groupID <- paste(isocline.df$species, isocline.df$scenario, sep= "-")
  
  #remove any NAs
  isocline.df <- isocline.df[complete.cases(isocline.df),]
  
  #and set any rates below zero to zero:
  isocline.df$IBI_infant_survival_rate[isocline.df$IBI_infant_surival_rate < 0] <- 0
  isocline.df$IBI_adult_survival_rate[isocline.df$IBI_adult_survival_rate<0] <- 0
  
  #add a column for multiplier
  isocline.df$multiplier <- as.factor(isocline.df$age_1st_rep/isocline.df$IBI)
  isocline.df$IBI <- factor(isocline.df$IBI)
  isocline.df$age_1st_rep <- factor(isocline.df$age_1st_rep)
  
  
  #and we plot
  if(do.plot==TRUE){
    
    isocline.df$species[isocline.df$species == "Astrochelys_yniphora"] <- "Astrochelys yniphora"
    dat1$species[dat1$species=="Astrochelys_yniphora"] <- "Astrochelys yniphora"
    #factor species by body size
    isocline.df$species <- factor(isocline.df$species, levels = c("Astrochelys yniphora"))
    
    #cull those that are over 1
    #isocline.df = subset(isocline.df, IBI_infant_survival_rate<=1)
    
    
    #and we'll shade above line and below, so let's add a mock variable here for each
    #these get reverse from mortality
    plot.lower.df <- subset(isocline.df, scenario=="best") #lowest survival rates to allow for persistence
    
    #and the lower 
    plot.upper.df <- subset(isocline.df, scenario=="worst")
   
    
    #want to fill between each lines
    #redo the levels
    isocline.df$species <- factor(isocline.df$species, levels= c("Astrochelys yniphora"))
    
    plot.upper.df$species <- factor(plot.upper.df$species, levels= c("Astrochelys yniphora"))
    
    
    #for max, make it just above the infant mort
    dat1 = subset(dat1, scenario=="best")
    #then overwrite Pteropus to use Eidolon
    #dat1$infant_annual_survival[dat1$species=="Pteropus rufus"] = dat1$infant_annual_survival[dat1$species=="Eidolon dupreanum"]
    
    quartz()
    
    p <-  ggplot(data=isocline.df) + facet_grid(species~.) +
      geom_ribbon(data =  plot.upper.df, aes(x= IBI_adult_survival_rate, ymin=IBI_infant_survival_rate, ymax=1), colour="black", fill="gray85") + #color above line
      geom_ribbon(data =  plot.lower.df, aes(x= IBI_adult_survival_rate, ymin=0, ymax=IBI_infant_survival_rate), fill="black", colour="black") + #color below line
      geom_line(aes(IBI_adult_survival_rate,IBI_infant_survival_rate, 
                    group=groupID), colour="white", alpha=0) +
      theme_classic()  +coord_cartesian(xlim = c(.25,1), ylim = c(0,1), expand=FALSE) +
      theme(strip.text = element_text(size=20, face = "italic"), 
            axis.title.x=element_text(size=20), axis.title.y=element_text(size=20), 
            axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), panel.spacing = unit(.5, "cm")) + 
      xlab("adult annual survival") + ylab("juvenile annual survival") +
      annotate("text", x=.94, y=.8, label="population\ngrowth", size=3.5, fontface="bold") + 
      annotate("text", x=.4, y=.5, label="population\nextinction",colour="white", size=3.5, fontface="bold") +
      geom_segment(data=dat1, aes(x=adult_IBI_survival, y=0, xend=adult_IBI_survival, yend=1), linetype=2, size=1, color="magenta") +
      geom_point(data=dat1, aes(x=adult_IBI_survival, y=infant_IBI_survival), shape=18, size=6, color="magenta")
    
    print(p)  
    
    if(do.save==TRUE){
      
      ggsave(file = filename,
             units="mm",  
             width=55, 
             height=35, 
             scale=3, 
             dpi=300)
    }
    
  }
  return(list(isocline.df, p))
}

wrap.isocline.survival(dat1=ang.dat, 
                       type="isocline", 
                       p_sequence=seq(0, 1, .01), 
                       do.plot=TRUE, 
                       do.save=TRUE, 
                       filename="Angonoka_phaseplane_plot.pdf")
