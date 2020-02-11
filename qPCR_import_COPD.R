####qPCR import COPD####
This code works now. You only have to upload new exports to the exports folder under data.


## qPCR import
library(qpcrpal); library(dplyr); library(qpcR); library(readxl); library(ggplot2); library(stringr); library(tidyr); 

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

# Remove targets that are not to be used (bad primers)
#best.models <- best.models %>%
  #filter(!(target %in% c("NEAT1 F2R2", "MTs16 F3R3", 
                         #"MTs12 F2R2", "MALAT1 F2R2", 
                         #"Lnc31 F3R3", "LincRAM F2R2", 
                         #"LincPINT F2R2", "LincP21 F3R3", 
                         #"H19 F1R1", "DRRRNA F3R3"))) %>%
  #print()

## load data with best model
qpcrbatch <- prepare_batch("./data/exports/", equipment = "quant", skip = 21) 

results <- list()

# Loop through all targets in best.models data frame
for(i in 1:nrow(best.models)){
  
  results[[i]] <- qpcrbatch %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>% # use the best model in each model_qpcr
    analyze_models() # analyze models for cpD2
  
}




 #combine all results and str split id variables
qpcrdat <- rbind_all(results)
id.var <- str_split_fixed(qpcrdat$ID, "_", 5) 
colnames(id.var) <- c("Ex_nr", "cdna", "timepoint", "gene", "leg")  
qpcrdat <- cbind(id.var, qpcrdat[,-1])




## estimate efficiencies ##
efficiencies <- list()

# use the same loop to analyze efficiencies
for(i in  1:nrow(best.models)){
  
  efficiencies[[i]] <- qpcrbatch %>%
    filter(target == best.models[i,1]) %>%
    model_qpcr(model = eval(parse(text = as.character(best.models[i,2]))), replicate = FALSE) %>%
    analyze_efficiency(method = "cpD2", model = "linexp", offset = -3)
  
}


# combine results and use str split to extract id variables

efficiencies <- rbind_all(efficiencies) 
id.var <- str_split_fixed(efficiencies$ID, "_", 5) 
colnames(id.var) <- c("Ex_nr", "cdna", "timepoint", "gene", "leg")  
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
  dplyr::select(Ex_nr, gene, cdna, cpD2, eff.y) %>%
  mutate(cq = cpD2,
         eff = eff.y) %>%
  ungroup() %>%
  data.frame()


### removed this code and then it worked:
#separate(target, into = c("target", "cdna"), sep = "_") %>%

####load sheet with FP and sets info####
library(tidyverse)

info <- read_excel("./data/RNA_stock_matrix.xlsx") %>%
  dplyr:: select(FP, Timepoint, Leg, RM_leg, Ex_nr)%>% 
  print()

####innerjoin on ex_nr####

qpcrdat.replicates <- qpcrdat.replicates %>%
  inner_join(info) %>%
  print() 

### Save data ####

saveRDS(qpcrdat.replicates, "./derivedData/qpcr.replicates_copd.Rds")


####write excel spreadsheet####
library(writexl)
write_xlsx(qpcrdat.replicates, "./derivedData/qpcr_r_output.xlsx")



# Remove objects from environment
rm(list = ls())
