####analysis of qPCR copd testrun####

## qPCR import
library(qpcrpal); library(dplyr); library(qpcR); library(readxl); library(ggplot2); library(stringr); library(tidyr)

### Prepare a batch of imported data 
batch <- prepare_batch("./data/exports/", equipment = "quant", skip = 21) %>%
  model_qpcr()

## Perform model comparisons
model.tests <- test_models(batch)


# Plot best models
data.frame(table(model.tests$best.model, model.tests$target)) %>%
  ggplot(aes(Var2, Freq, fill = Var1)) + geom_bar(stat = "identity") + coord_flip()


# Best model per target are stored for modeling with best model
# Plot best model per target
data.frame(table(model.tests$best.model, model.tests$target)) %>%
  ggplot(aes(Var2, Freq, fill = Var1)) + geom_bar(stat="identity") + coord_flip()

# Best model per target are stored for modeling with best model
best.models <- data.frame(target = names(apply(table(model.tests$best.model, model.tests$target), 2, which.max)),
                          model = as.character(row.names(table(model.tests$best.model, model.tests$target))[apply(table(model.tests$best.model, model.tests$target), 2, which.max)]))

## load data with best model
qpcrbatch <- prepare_batch("./data/exports/", equipment = "quant", skip = 21) 

results <- list()

# Loop through all targets in best.models data frame
for(i in 1:nrow(best.models)){
  
  results[[i]] <- qpcrbatch %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>% # use the best model in each model_qpcr
    analyze_models() # analyze models for cpD2
  
}


# combine all results and str split id variables
qpcrdat <- rbind_all(results) 
id.var <- str_split_fixed(qpcrdat$ID, "_", 5) 
colnames(id.var) <- c("subject", "cdna", "timepoint", "gene","leg")  
qpcrdat <- cbind(id.var, qpcrdat[,-1])



## estimate efficiencies ##
efficiencies <- list()

# use the same loop to analyze efficiencies
for(i in  1:nrow(best.models)){
  
  efficiencies[[i]] <- qpcrbatch %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>%
    analyze_efficiency(method = "cpD2", model = "linexp", offset = -3)
  
}


# combine results and use str split to extract id variables

efficiencies <- rbind_all(efficiencies) 
id.var <- str_split_fixed(efficiencies$ID, "_", 5) 
colnames(id.var) <- c("subject", "cdna", "timepoint", "gene","leg")  
efficiencies <- cbind(id.var, efficiencies[,-1])

efficiencies %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  #separate(target, into = c("target", "cDNA"), sep = "_") %>%
  group_by(gene)%>%
  summarise(efficiency = mean(eff, na.rm = TRUE),
            max.eff = max(eff, na.rm = TRUE),
            min.eff = min(eff, na.rm = TRUE),
            sd.eff = sd(eff, na.rm = TRUE))%>%
  ggplot(aes(gene, efficiency)) + geom_point() + 
  coord_flip() 

#MYHC2a is removed in this step

effs <- efficiencies %>%
  filter(eff > 1.5 & eff < 2.5)%>% # remove outliers from efficiency estimation
  #separate(target, into = c("target", "cdna"), sep = "_") %>%
  group_by(gene)%>%
  summarise(eff = mean(eff, na.rm = TRUE))



## Extract information on replicates

replicates <- data.frame(str_split(qpcrdat$gene, "_", simplify = TRUE))

colnames(replicates) <- c("gene")


qpcrdat <- cbind(replicates, qpcrdat[, -4])

head(qpcrdat)
qpcrdat

## Combine all qPCR parameters in the study to a data.frame containing all replicates
qpcrdat.replicates <- qpcrdat %>%
  inner_join(effs, by = "gene") %>%
  dplyr::select(subject, gene, cdna, cpD2, eff.y) %>%
  mutate(cq = cpD2,
         eff = eff.y) %>%
  ungroup() %>%
  data.frame()

qdat <- qpcrdat.replicates %>%
  mutate(Ra = eff^-cq) 


#make model#

m1 <- lm(Ra ~ gene + cdna + subject, data = qdat)

library(emmeans); library(ggplot2)

estimates <- emmeans(m1, specs = ~ "gene|cdna+subject")


estimates %>%
  data.frame() %>% 
  ggplot(aes(cdna, emmean, color = subject, group = gene)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), 
                width = 0.2, position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  facet_grid(gene ~ ., scales = "free" )



####expression of cq values####

qdat_analysis <- qdat %>% 
  dplyr::select(subject, gene, cdna, cq) %>%
  group_by(cdna, gene) %>%
  summarise(m = mean(cq, na.rm = TRUE), 
            s = sd(cq, na.rm = TRUE)) 
  

qdat_analysis %>% 
  ggplot(aes(cdna, m, color = gene)) +
  geom_errorbar(aes(ymin = m-s, ymax = m+s),
                width = 0.2,
                position = position_dodge(width = 0.2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_line(position = position_dodge(width = 0.2)) +
  facet_grid(gene ~ ., scales = "free" )




