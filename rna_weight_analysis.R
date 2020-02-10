### Total RNA analysis COPD ####

# 2019-09-09: DH 

#install packages
install.packages("emmeans")


# Libraries
library(readxl); library(tidyverse); library(nlme); library(emmeans)


### Read data 



rna_dat <- read_excel("./data/RNA_extraction.xlsx", sheet = 1, na = "") %>%
        dplyr::select(FP, Syk, Placebo, Leg, RM_leg, Timepoint,
                      weight = wet_weight_extracted,
                      stock,
                      conc.stock_ngRNA_ul_1:conc.stock_ngRNA_ul_4) %>%
        gather(replicate, rna, conc.stock_ngRNA_ul_1:conc.stock_ngRNA_ul_4) %>%
        filter(!is.na(rna) & rna > 0) %>% # Removes records with loss of RNA in extraction
        mutate(replicate = paste0("rep", gsub("conc.stock_ngRNA_ul_", "", replicate)),
               elution.volume = if_else(as.numeric(stock) < 20, 20, 30),
               rna = rna * elution.volume) %>%

        print()


# Visualize muscle weight to RNA relationship
rna_dat %>%
        ggplot(aes(log(weight), log(rna))) + geom_point() + geom_smooth(method = "lm")

# There might be sample loss, use the below code to create an outlier variable

# Create linear model 
mod <- with(rna_dat, lm(log(rna) ~ log(weight)))


# Calculates standardized residuals
rna_dat$resid<- resid(mod)/sd(resid(mod))
# store predicted values for diagnostic plotting
rna_dat$pred<- predict(mod)

rna_dat %>%
        ggplot(aes(pred, resid)) + geom_point() +
        scale_y_continuous(breaks = c(-20, -10, -5, -4, 0, 4, 5), 
                           limits = c(-20, 5))

## Adding outlier variable
rna_dat %>%
        mutate(outlier = if_else(resid < -4 | resid > 4 , "out", "in")) %>%# 4 units deviance from zero
        ggplot(aes(pred, resid, color = outlier)) + geom_point() +
        scale_y_continuous(breaks = c(-20, -10, -5, -4, 0, 4, 5), 
                           limits = c(-20, 5))

rna <- rna_dat %>%
        filter(RM_leg %in% c("30", "10")) %>%
        mutate(outlier = if_else(resid < -4 | resid > 4 , "out", "in"), # 4 units deviance from zero
               RM_leg = factor(paste0("RM", RM_leg), levels = c("RM10", "RM30")),  
               Timepoint = factor(Timepoint, levels = c("PreSupp","PreExc", "ThreeW", "PostExc")),
               copd = factor(if_else(Syk == 1, "copd", "healthy"), levels = c("healthy", "copd")),
               vitd = factor(if_else(Placebo == 0, "placebo", "vitd"), levels = c("placebo", "vitd")),
               sample = paste0(FP, Leg, Timepoint), 
               time_copd = paste0(copd, Timepoint),
               rna.weight = rna/weight) %>%
print()

####### Modelling ###########


rna %>%
        filter(RM_leg == "RM10" & copd == "copd") %>%
        group_by(FP, Timepoint, Leg) %>%
        summarise(rna = mean(rna.weight, na.rm = TRUE)) %>%
        ggplot(aes(Timepoint, rna, group = paste(FP, Leg), color = as.factor(FP))) + geom_point() + geom_line()




rna %>%
        group_by(Timepoint, FP, RM_leg, copd, vitd) %>%
        summarise(rna = mean(rna.weight, na.rm = TRUE)) %>%
        ggplot(aes(Timepoint, rna, group = paste(FP, RM_leg), color = RM_leg)) + 
        geom_point() + 
        geom_line() +
        geom_boxplot(aes(group = NULL, color = RM_leg)) +
        facet_grid(copd ~ vitd)


### Is there a supplementation effect?

rna_presup <- rna %>%
        filter(Timepoint %in% c("PreSupp", "PreExc")) %>%
        dplyr::select(FP, Leg, Timepoint, copd, vitd, rna.weight, replicate) %>%
        group_by(FP, Leg, Timepoint, copd, vitd) %>%
        summarise(rna.weight = mean(rna.weight, na.rm = TRUE)) %>%
        spread(Timepoint, rna.weight) %>%
        filter(!is.na(PreSupp)) %>%
        gather(Timepoint, rna.weight, PreSupp:PreExc) %>%
        mutate(Timepoint = factor(Timepoint, levels = c("PreSupp", "PreExc"))) %>%
        print()




m1 <- lme(log(rna.weight) ~ Timepoint * copd * vitd, 
          random = list(FP = ~1),
          data = rna_presup, 
          control=list(msMaxIter=120,
                       opt = "nloptwrap",
                       msVerbose=TRUE), method = "ML",
          na.action = na.omit)


summary(m1) # No evidence for keeping presupp


rna_ex <- rna %>%
        filter(Timepoint != "PreSupp") %>%
        print()
        



m1 <- lme(log(rna.weight) ~ Timepoint * RM_leg * copd * vitd, 
          random = list(FP = pdDiag(~1 + Timepoint),
                        sample = ~1),
          data = rna_ex %>%
                  filter(outlier == "in"), 
          control=list(msMaxIter=120,
                       opt = "nloptwrap",
                       msVerbose=TRUE), method = "ML")





plot(m1, resid(.)~fitted(.)|Timepoint + copd) # Evidence for unequal variance per timepoint?

# Updates model 1 with varIdent function
m2 <- update(m1, weights = varIdent(form = ~1|Timepoint))
m3 <- update(m1, weights = varPower(form = ~fitted(.)|Timepoint))
m4 <- update(m1, weights = varExp(form = ~fitted(.)|Timepoint))


anova(m1, m2, m3) # Evidence for including variance function per timepoint, varPower is best
intervals(m3)

qqnorm(resid(m3), main = "log-values"); qqline(resid(m3))

plot(m3, resid(.)~fitted(.)|Timepoint) # Evidence for unequal variance per timepoint?



summary(m3)




# Removing vitd
m3.1 <- update(m3, fixed = log(rna.weight) ~ (Timepoint * RM_leg) * copd)


# Removing copd
m3.2 <- update(m3, fixed = log(rna.weight) ~ Timepoint * RM_leg * vitd)


anova( m3.1,m3,  m3.2)


# Removing vitd
m3.3 <- update(m3, fixed = log(rna.weight) ~ Timepoint * RM_leg)

anova(m3.1, m3.3)

summary(m3.1)




em_m3 <- emmeans(m3,
                 specs = ~"Timepoint|RM_leg|copd|vitd")

em_m3.1 <- emmeans(m3.1,
                 specs = ~"Timepoint|RM_leg|copd")


em_m3 %>%
        data.frame() %>%
        ggplot(aes(Timepoint, emmean, 
                   group = paste(RM_leg, copd, vitd),
                   color = paste(RM_leg, copd, vitd))) +
        geom_point(position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) + 
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      width = 0.2, 
                      position = position_dodge(width = 0.2))



em_m3.1 %>%
        data.frame() %>%
        ggplot(aes(Timepoint, emmean, 
                   group = paste(RM_leg, copd),
                   color = paste(RM_leg, copd))) +
        geom_point(position = position_dodge(width = 0.2)) +
        geom_line(position = position_dodge(width = 0.2)) + 
        geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                      width = 0.2, 
                      position = position_dodge(width = 0.2))




m2 <- lme(log(rna) ~ Timepoint * vitd * copd, 
          random = list(FP = ~1,
                        sample = ~1),
          data = rna %>%
                  filter(outlier == "in"), 
          control=list(msMaxIter=120,
                       opt = "nloptwrap",
                       msVerbose=TRUE), method = "REML")





qqnorm(resid(m1), main = "total.rna.model log-values"); qqline(resid(m1))






plot(m1)

summary(m1)













