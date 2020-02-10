####qPCR_analysis_COPD####

# Analysis of qpcr data #
library(readxl); library(tidyverse); library(knitr); library(emmeans);library(dplyr)
library(nlme)
#load data 

qpcr.dat <- readRDS("./derivedData/qpcr.replicates_copd.Rds")


qdat <- qpcr.dat %>%
  mutate(Ra = eff^-cq) 
         #technical = paste0(subject, cdna), 
         #biological = paste0(subject)) %>%
  #print()


####load sheet with FP and sets info####
library(tidyverse)

info <- read_excel("./data/RNA_stock_matrix.xlsx") %>%
  dplyr:: select(FP, Timepoint, Leg, RM_leg, Ex_nr)%>% 
  print()

####innerjoin on ex_nr####

qpcrdat1 <- qpcrdat %>%
  inner_join(info) %>%
  print() 


#load sheet with FP and sets info

#info <- read.csv("./data/oneThreeSetLeg.csv", sep = ";") %>%
 # pivot_longer(names_to = "sets", 
               #values_to = "leg", multiple:single) %>%
 # print()
#
#qdat <- qdat %>%
 # inner_join(info) %>%
  #print()



# make model

m1 <- lme(log(Ra) ~ 0 + gene + gene:cdna, data = qdat)
          #random = list(subject = ~ 1, technical = ~ 1), 
          #data = qdat) 
          #,control=list(msMaxIter=120,
                      # opt = "nloptwrap", msVerbose=TRUE),
          #method = "REML", na.action = na.exclude)

##### per total-RNA models ######

# Fixed and random effects formulas for the first step of model building
# Compared to Method paper (Hammarstr��t al 2018), the fixed effects are reuced
# to only contain gene-specific time + time:sets



## m1 is the model assuming homoscedastic errors.
m1 <- lme(fixed, random = random, data = qdat.full,
          control=list(msMaxIter=120,
                      opt = "nloptwrap", msVerbose=TRUE),
          method = "REML", na.action = na.exclude) # This allows to check progress


#### Models allowing for heteroscedasticity in the residuals. ####
# The variance is believed to be different from gene to gene due to different expression, primer design etc.
# One can also expect that residual variance to be related to cq-values as higher cq-values 
# implies higher meassurement error due to stochastic variation in the beginning of amplification.

# Variance functions:
varfun1 <- varIdent(form = ~1|target) # Heterogeneity per gene
varfun2 <- varPower(form = ~cq|target) # Power of the variance covariate cq, per gene stratum
varfun3 <- varFixed(~cq) # Variance increases as a function of cq
varfun4 <- varExp(form =~cq|target) # Variance changes as a exponentil of cq-values
varfun5 <- varConstPower(form = ~cq|target)

m2 <- update(m1, weights = varfun1)
# m3 <- update(m1, weights = varfun2)
m4 <- update(m1, weights = varfun3)
m5 <- update(m1, weights = varfun4)
# m6 <- update(m1, weights = varfun5)

### Test models
anova(m1, m2, m5)


m7 <- update(m5, random = list(subject = ~ 1 + target, 
                               biological = ~ 1,
                               technical = ~ 1))

library(emmeans); library(ggplot2)

estimates <- emmeans(m5, specs = ~"target|sets+timepoint")


estimates %>%
  data.frame() %>% 
  mutate(timepoint = factor(timepoint, levels = c("w0", "w2pre", "w2post", "w12"))) %>%
  ggplot(aes(timepoint, emmean, color = sets, group = paste(target, sets))) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  facet_grid(target ~ ., scales = "free" )
