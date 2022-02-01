### Forest plot comparing the agnostic and mechanistic KTSP in the CV part using all metrics

rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)

## Agnostic XGB Average performance

# Load the agnostic models (IndGenes)
AgnosticEMTAB_Out_IndGenes <- load("./Objs/XGB/EMTAB_Out_XGB_IndvGenes_AgnosticPerformance.rda")
AgnosticGSE13507_Out_IndGenes <- load("./Objs/XGB/GSE13507_Out_XGB_IndvGenes_AgnosticPerformance.rda")
AgnosticGSE32894_Out_IndGenes <- load("./Objs/XGB/GSE32894_Out_XGB_IndvGenes_AgnosticPerformance.rda")
AgnosticGSE57813_Out_IndGenes <- load("./Objs/XGB/GSE57813_Out_XGB_IndvGenes_AgnosticPerformance.rda")
AgnosticPMID_Out_IndGenes <- load("./Objs/XGB/PMID_Out_XGB_IndvGenes_AgnosticPerformance.rda")

# Load the agnostic models (Pairs)
AgnosticEMTAB_Out_Pairs <- load("./Objs/XGB/EMTAB_Out_XGB_AgnosticPerformance.rda")
AgnosticGSE13507_Out_Pairs <- load("./Objs/XGB/GSE13507_Out_XGB_AgnosticPerformance.rda")
AgnosticGSE32894_Out_Pairs <- load("./Objs/XGB/GSE32894_Out_XGB_AgnosticPerformance.rda")
AgnosticGSE57813_Out_Pairs <- load("./Objs/XGB/GSE57813_Out_XGB_AgnosticPerformance.rda")
AgnosticPMID_Out_Pairs <- load("./Objs/XGB/PMID_Out_XGB_AgnosticPerformance.rda")

# Load the Mech models
MechEMTAB_Out <- load("./Objs/XGB/EMTAB_Out_XGB_MechPerformance.rda")
MechGSE13507_Out <- load("./Objs/XGB/GSE13507_Out_XGB_MechPerformance.rda")
MechGSE32894_Out <- load("./Objs/XGB/GSE32894_Out_XGB_MechPerformance.rda")
MechGSE57813_Out <- load("./Objs/XGB/GSE57813_Out_XGB_MechPerformance.rda")
MechPMID_Out <- load("./Objs/XGB/PMID_Out_XGB_MechPerformance.rda")




# # Mean and SD of the Agnostic AUC in the training data
# MeanAUCTrainingAgnostic <- mean(c(EMTAB_Out_RF_AgnosticPerformance$Training[1], 
#                         GSE13507_Out_RF_AgnosticPerformance$Training[1],
#                         GSE32894_Out_RF_AgnosticPerformance$Training[1],
#                         GSE57813_Out_RF_AgnosticPerformance$Training[1],
#                         PMID_Out_RF_AgnosticPerformance$Training[1]
#                         ))
# 
# SDAUCTrainingAgnostic <- sd(c(EMTAB_Out_RF_AgnosticPerformance$Training[1], 
#                           GSE13507_Out_RF_AgnosticPerformance$Training[1],
#                           GSE32894_Out_RF_AgnosticPerformance$Training[1],
#                           GSE57813_Out_RF_AgnosticPerformance$Training[1],
#                           PMID_Out_RF_AgnosticPerformance$Training[1]
#                          ))
# 
# # Mean and SD of the Agnostic AUC in the testing data
# MeanAUCTestingAgnostic <- mean(c(EMTAB_Out_RF_AgnosticPerformance$Testing[1], 
#                          GSE13507_Out_RF_AgnosticPerformance$Testing[1],
#                          GSE32894_Out_RF_AgnosticPerformance$Testing[1],
#                          GSE57813_Out_RF_AgnosticPerformance$Testing[1],
#                          PMID_Out_RF_AgnosticPerformance$Testing[1]
#                         ))
# 
# 
# SDAUCTestingAgnostic <- sd(c(EMTAB_Out_RF_AgnosticPerformance$Testing[1], 
#                                  GSE13507_Out_RF_AgnosticPerformance$Testing[1],
#                                  GSE32894_Out_RF_AgnosticPerformance$Testing[1],
#                                  GSE57813_Out_RF_AgnosticPerformance$Testing[1],
#                                  PMID_Out_RF_AgnosticPerformance$Testing[1]
#                         ))

##############################################################
## Mech RF Average performance


# Mean and SD of the Mechanistic AUC in the training data
# MeanAUCTrainingMech <- mean(c(EMTAB_Out_RF_MechPerformance$Training[1], 
#                                   GSE13507_Out_RF_MechPerformance$Training[1],
#                                   GSE32894_Out_RF_MechPerformance$Training[1],
#                                   GSE57813_Out_RF_MechPerformance$Training[1],
#                                   PMID_Out_RF_MechPerformance$Training[1]
#                         ))
# 
# SDAUCTrainingMech <- sd(c(EMTAB_Out_RF_MechPerformance$Training[1], 
#                               GSE13507_Out_RF_MechPerformance$Training[1],
#                               GSE32894_Out_RF_MechPerformance$Training[1],
#                               GSE57813_Out_RF_MechPerformance$Training[1],
#                               PMID_Out_RF_MechPerformance$Training[1]
#                         ))
# 
# # Mean and SD of the Mechanistic AUC in the training data
# MeanAUCTestingMech <- mean(c(EMTAB_Out_RF_MechPerformance$Testing[1], 
#                                  GSE13507_Out_RF_MechPerformance$Testing[1],
#                                  GSE32894_Out_RF_MechPerformance$Testing[1],
#                                  GSE57813_Out_RF_MechPerformance$Testing[1],
#                                  PMID_Out_RF_MechPerformance$Testing[1]
#                         ))
# 
# 
# SDAUCTestingMech <- sd(c(EMTAB_Out_RF_MechPerformance$Testing[1], 
#                              GSE13507_Out_RF_MechPerformance$Testing[1],
#                              GSE32894_Out_RF_MechPerformance$Testing[1],
#                              GSE57813_Out_RF_MechPerformance$Testing[1],
#                              PMID_Out_RF_MechPerformance$Testing[1]
#                        ))



###################

# EMTAB_out
EMTAB_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("EMTAB", nrow(EMTAB_Out_XGB_IndvGenes_AgnosticPerformance))
EMTAB_Out_XGB_AgnosticPerformance$out_study <- rep("EMTAB", nrow(EMTAB_Out_XGB_AgnosticPerformance))
EMTAB_Out_XGB_MechPerformance$out_study <- rep("EMTAB", nrow(EMTAB_Out_XGB_MechPerformance))

EMTAB_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(EMTAB_Out_XGB_IndvGenes_AgnosticPerformance))
EMTAB_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(EMTAB_Out_XGB_AgnosticPerformance))
EMTAB_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(EMTAB_Out_XGB_MechPerformance))

# GSE13507_out
GSE13507_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("GSE13507", nrow(GSE13507_Out_XGB_IndvGenes_AgnosticPerformance))
GSE13507_Out_XGB_AgnosticPerformance$out_study <- rep("GSE13507", nrow(GSE13507_Out_XGB_AgnosticPerformance))
GSE13507_Out_XGB_MechPerformance$out_study <- rep("GSE13507", nrow(GSE13507_Out_XGB_MechPerformance))

GSE13507_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(GSE13507_Out_XGB_IndvGenes_AgnosticPerformance))
GSE13507_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(GSE13507_Out_XGB_AgnosticPerformance))
GSE13507_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE13507_Out_XGB_MechPerformance))

# GSE32894_out
GSE32894_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("GSE32894", nrow(GSE32894_Out_XGB_IndvGenes_AgnosticPerformance))
GSE32894_Out_XGB_AgnosticPerformance$out_study <- rep("GSE32894", nrow(GSE32894_Out_XGB_AgnosticPerformance))
GSE32894_Out_XGB_MechPerformance$out_study <- rep("GSE32894", nrow(GSE32894_Out_XGB_MechPerformance))

GSE32894_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(GSE32894_Out_XGB_IndvGenes_AgnosticPerformance))
GSE32894_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(GSE32894_Out_XGB_AgnosticPerformance))
GSE32894_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE32894_Out_XGB_MechPerformance))

# GSE57813_out
GSE57813_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("GSE57813", nrow(GSE57813_Out_XGB_IndvGenes_AgnosticPerformance))
GSE57813_Out_XGB_AgnosticPerformance$out_study <- rep("GSE57813", nrow(GSE57813_Out_XGB_AgnosticPerformance))
GSE57813_Out_XGB_MechPerformance$out_study <- rep("GSE57813", nrow(GSE57813_Out_XGB_MechPerformance))

GSE57813_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(GSE57813_Out_XGB_IndvGenes_AgnosticPerformance))
GSE57813_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(GSE57813_Out_XGB_AgnosticPerformance))
GSE57813_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE57813_Out_XGB_MechPerformance))

# PMID_out
PMID_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("PMID", nrow(PMID_Out_XGB_IndvGenes_AgnosticPerformance))
PMID_Out_XGB_AgnosticPerformance$out_study <- rep("PMID", nrow(PMID_Out_XGB_AgnosticPerformance))
PMID_Out_XGB_MechPerformance$out_study <- rep("PMID", nrow(PMID_Out_XGB_MechPerformance))

PMID_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(PMID_Out_XGB_IndvGenes_AgnosticPerformance))
PMID_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(PMID_Out_XGB_AgnosticPerformance))
PMID_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(PMID_Out_XGB_MechPerformance))


All_XGB_CV <- rbind(EMTAB_Out_XGB_IndvGenes_AgnosticPerformance, EMTAB_Out_XGB_AgnosticPerformance, EMTAB_Out_XGB_MechPerformance,
                    GSE13507_Out_XGB_IndvGenes_AgnosticPerformance, GSE13507_Out_XGB_AgnosticPerformance, GSE13507_Out_XGB_MechPerformance,
                    GSE32894_Out_XGB_IndvGenes_AgnosticPerformance, GSE32894_Out_XGB_AgnosticPerformance, GSE32894_Out_XGB_MechPerformance, 
                    GSE57813_Out_XGB_IndvGenes_AgnosticPerformance, GSE57813_Out_XGB_AgnosticPerformance, GSE57813_Out_XGB_MechPerformance, 
                    PMID_Out_XGB_IndvGenes_AgnosticPerformance, PMID_Out_XGB_AgnosticPerformance, PMID_Out_XGB_MechPerformance)

All_XGB_CV$metric <- rownames(All_XGB_CV)
All_XGB_CV$metric <- gsub('[[:digit:]]+', '', All_XGB_CV$metric )


##########################################################
##########################################################
## Calculate the average statistics (Across the 5 repeats)

#####
# average stats in the agnostic models (Ind genes)

## AUC
# Mean AUC in both the training and testing data
MeanAUC_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                             filter(model_type == "Agnostic_IndGenes") %>%
                             filter(metric == "AUC") %>%
                             select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                            filter(model_type == "Agnostic_IndGenes") %>%
                            filter(metric == "AUC") %>%
                            select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                  filter(model_type == "Agnostic_IndGenes") %>%
                                  filter(metric == "Accuracy") %>%
                                  select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                 filter(model_type == "Agnostic_IndGenes") %>%
                                 filter(metric == "Accuracy") %>%
                                 select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                     filter(model_type == "Agnostic_IndGenes") %>%
                                     filter(metric == "Bal.Accuracy") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                    filter(model_type == "Agnostic_IndGenes") %>%
                                    filter(metric == "Bal.Accuracy") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                     filter(model_type == "Agnostic_IndGenes") %>%
                                     filter(metric == "Sensitivity") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                    filter(model_type == "Agnostic_IndGenes") %>%
                                    filter(metric == "Sensitivity") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                     filter(model_type == "Agnostic_IndGenes") %>%
                                     filter(metric == "Specificity") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                                    filter(model_type == "Agnostic_IndGenes") %>%
                                    filter(metric == "Specificity") %>%
                                    select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                             filter(model_type == "Agnostic_IndGenes") %>%
                             filter(metric == "MCC") %>%
                             select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_IndGenes <-  apply(All_XGB_CV %>%
                            filter(model_type == "Agnostic_IndGenes") %>%
                            filter(metric == "MCC") %>%
                            select(c("Training", "Testing")), 2, sd)



## Group them together
MeansAgnostic_IndGenes <- rbind(MeanAUC_Agnostic_IndGenes, SD_AUC_Agnostic_IndGenes, 
                       MeanAccuracy_Agnostic_IndGenes, SD_Accuracy_Agnostic_IndGenes, 
                       MeanBalAccuracy_Agnostic_IndGenes, SD_BalAccuracy_Agnostic_IndGenes, 
                       MeanSensitivity_Agnostic_IndGenes, SD_Sensitivity_Agnostic_IndGenes, 
                       MeanSpecificity_Agnostic_IndGenes, SD_Specificity_Agnostic_IndGenes,
                       MeanMCC_Agnostic_IndGenes, SD_MCC_Agnostic_IndGenes)
# Round the Numbers
MeansAgnostic_IndGenes <- apply(MeansAgnostic_IndGenes, 2, round, digits = 2) 

################################################################
## Agnostic pairs

MeanAUC_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                      filter(model_type == "Agnostic_Pairs") %>%
                                      filter(metric == "AUC") %>%
                                      select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                     filter(model_type == "Agnostic_Pairs") %>%
                                     filter(metric == "AUC") %>%
                                     select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                           filter(model_type == "Agnostic_Pairs") %>%
                                           filter(metric == "Accuracy") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                          filter(model_type == "Agnostic_Pairs") %>%
                                          filter(metric == "Accuracy") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                              filter(model_type == "Agnostic_Pairs") %>%
                                              filter(metric == "Bal.Accuracy") %>%
                                              select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                             filter(model_type == "Agnostic_Pairs") %>%
                                             filter(metric == "Bal.Accuracy") %>%
                                             select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                              filter(model_type == "Agnostic_Pairs") %>%
                                              filter(metric == "Sensitivity") %>%
                                              select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                             filter(model_type == "Agnostic_Pairs") %>%
                                             filter(metric == "Sensitivity") %>%
                                             select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                              filter(model_type == "Agnostic_Pairs") %>%
                                              filter(metric == "Specificity") %>%
                                              select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                             filter(model_type == "Agnostic_Pairs") %>%
                                             filter(metric == "Specificity") %>%
                                             select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                      filter(model_type == "Agnostic_Pairs") %>%
                                      filter(metric == "MCC") %>%
                                      select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_Pairs <-  apply(All_XGB_CV %>%
                                     filter(model_type == "Agnostic_Pairs") %>%
                                     filter(metric == "MCC") %>%
                                     select(c("Training", "Testing")), 2, sd)



## Group them together
MeansAgnostic_Pairs <- rbind(MeanAUC_Agnostic_Pairs, SD_AUC_Agnostic_Pairs, 
                                MeanAccuracy_Agnostic_Pairs, SD_Accuracy_Agnostic_Pairs, 
                                MeanBalAccuracy_Agnostic_Pairs, SD_BalAccuracy_Agnostic_Pairs, 
                                MeanSensitivity_Agnostic_Pairs, SD_Sensitivity_Agnostic_Pairs, 
                                MeanSpecificity_Agnostic_Pairs, SD_Specificity_Agnostic_Pairs,
                                MeanMCC_Agnostic_Pairs, SD_MCC_Agnostic_Pairs)
# Round the Numbers
MeansAgnostic_Pairs <- apply(MeansAgnostic_Pairs, 2, round, digits = 2) 

############################################
##
#####
# average stats in the Mechanistic models

## AUC
# Mean AUC in both the training and testing data
MeanAUC_Mechanistic <-  apply(All_XGB_CV %>%
                                filter(model_type == "Mechanistic") %>%
                                filter(metric == "AUC") %>%
                                select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Mechanistic <-  apply(All_XGB_CV %>%
                               filter(model_type == "Mechanistic") %>%
                               filter(metric == "AUC") %>%
                               select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Mechanistic <-  apply(All_XGB_CV %>%
                                     filter(model_type == "Mechanistic") %>%
                                     filter(metric == "Accuracy") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Mechanistic <-  apply(All_XGB_CV %>%
                                    filter(model_type == "Mechanistic") %>%
                                    filter(metric == "Accuracy") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Mechanistic <-  apply(All_XGB_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Bal.Accuracy") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Mechanistic <-  apply(All_XGB_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Bal.Accuracy") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Mechanistic <-  apply(All_XGB_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Sensitivity") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Mechanistic <-  apply(All_XGB_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Sensitivity") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Mechanistic <-  apply(All_XGB_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Specificity") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Mechanistic <-  apply(All_XGB_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Specificity") %>%
                                       select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Mechanistic <-  apply(All_XGB_CV %>%
                                filter(model_type == "Mechanistic") %>%
                                filter(metric == "MCC") %>%
                                select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Mechanistic <-  apply(All_XGB_CV %>%
                               filter(model_type == "Mechanistic") %>%
                               filter(metric == "MCC") %>%
                               select(c("Training", "Testing")), 2, sd)



## Group them together
MeansMechanistic <- rbind(MeanAUC_Mechanistic, SD_AUC_Mechanistic, 
                          MeanAccuracy_Mechanistic, SD_Accuracy_Mechanistic, 
                          MeanBalAccuracy_Mechanistic, SD_BalAccuracy_Mechanistic, 
                          MeanSensitivity_Mechanistic, SD_Sensitivity_Mechanistic, 
                          MeanSpecificity_Mechanistic, SD_Specificity_Mechanistic,
                          MeanMCC_Mechanistic, SD_MCC_Mechanistic)
# Round the Numbers
MeansMechanistic <- apply(MeansMechanistic, 2, round, digits = 2) 

######################
## Bind mean agnostic and mean mechanistic together
colnames(MeansAgnostic_IndGenes) <- paste("Agnostic_IndGenes", colnames(MeansAgnostic_IndGenes), sep = "_")
colnames(MeansAgnostic_Pairs) <- paste("Agnostic_Pairs", colnames(MeansAgnostic_Pairs), sep = "_")
colnames(MeansMechanistic) <- paste("Mechanistic", colnames(MeansMechanistic), sep = "_")

# Modify the rownames
rownames(MeansAgnostic_IndGenes) <- c("Mean_AUC", "SD_AUC", "Mean_Accuracy", "SD_Accuracy", "Mean_BalAccuracy", "SD_BalAccuracy", "Mean_Sensitivity", "SD_Sensitivity", "Mean_Specificity", "SD_Specificity", "Mean_MCC", "SD_MCC")
rownames(MeansAgnostic_Pairs) <- c("Mean_AUC", "SD_AUC", "Mean_Accuracy", "SD_Accuracy", "Mean_BalAccuracy", "SD_BalAccuracy", "Mean_Sensitivity", "SD_Sensitivity", "Mean_Specificity", "SD_Specificity", "Mean_MCC", "SD_MCC")
rownames(MeansMechanistic) <- c("Mean_AUC", "SD_AUC", "Mean_Accuracy", "SD_Accuracy", "Mean_BalAccuracy", "SD_BalAccuracy", "Mean_Sensitivity", "SD_Sensitivity", "Mean_Specificity", "SD_Specificity", "Mean_MCC", "SD_MCC")

XGB_CrossStudyValidation <- cbind(MeansAgnostic_IndGenes, MeansAgnostic_Pairs, MeansMechanistic)
save(XGB_CrossStudyValidation, file = "./Objs/XGB/XGB_CrossStudyValidation.rda")

load("./Objs/XGB/XGB_CrossStudyValidation.rda")
############################
## calculate the 95% CI 

### Agnostic models (Ind Genes)

## Training data

## CIs for AUC
AUCs_Train_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "AUC", model_type == "Agnostic_IndGenes") %>%
  select(Training)

CI_AUCs_Train_Agnostic_IndGenes <- quantile(AUCs_Train_Agnostic_IndGenes$Training, c(0.025, 0.975))
CI_AUCs_Train_Agnostic_IndGenes

# CIs for Accuracy
Accuracy_Train_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Accuracy", model_type == "Agnostic_IndGenes") %>%
  select(Training)

CI_Accuracy_Train_Agnostic_IndGenes <- quantile(Accuracy_Train_Agnostic_IndGenes$Training, c(0.025, 0.975))
CI_Accuracy_Train_Agnostic_IndGenes

# CIs for Bal. Accuracy
Bal.Accuracy_Train_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Agnostic_IndGenes") %>%
  select(Training)


CI_Bal.Accuracy_Train_Agnostic_IndGenes <- quantile(Bal.Accuracy_Train_Agnostic_IndGenes$Training, c(0.025, 0.975))
CI_Bal.Accuracy_Train_Agnostic_IndGenes

# CIs for sensitivity
Sensitivity_Train_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Sensitivity", model_type == "Agnostic_IndGenes") %>%
  select(Training)

CI_Sensitivity_Train_Agnostic_IndGenes <- quantile(Sensitivity_Train_Agnostic_IndGenes$Training, c(0.025, 0.975))
CI_Sensitivity_Train_Agnostic_IndGenes

# CIs for Specificity
Specificity_Train_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Specificity", model_type == "Agnostic_IndGenes") %>%
  select(Training)

CI_Specificity_Train_Agnostic_IndGenes <- quantile(Specificity_Train_Agnostic_IndGenes$Training, c(0.025, 0.975))
CI_Specificity_Train_Agnostic_IndGenes

# CIs for MCC
MCC_Train_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "MCC", model_type == "Agnostic_IndGenes") %>%
  select(Training)

CI_MCC_Train_Agnostic_IndGenes <- quantile(MCC_Train_Agnostic_IndGenes$Training, c(0.025, 0.975))
CI_MCC_Train_Agnostic_IndGenes

####################
## Testing data

## CIs for AUC
AUCs_Test_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "AUC", model_type == "Agnostic_IndGenes") %>%
  select(Testing)

CI_AUCs_Test_Agnostic_IndGenes <- quantile(AUCs_Test_Agnostic_IndGenes$Testing, c(0.025, 0.975))
CI_AUCs_Test_Agnostic_IndGenes

# CIs for Accuracy
Accuracy_Test_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Accuracy", model_type == "Agnostic_IndGenes") %>%
  select(Testing)

CI_Accuracy_Test_Agnostic_IndGenes <- quantile(Accuracy_Test_Agnostic_IndGenes$Testing, c(0.025, 0.975))
CI_Accuracy_Test_Agnostic_IndGenes

# CIs for Bal. Accuracy
Bal.Accuracy_Test_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Agnostic_IndGenes") %>%
  select(Testing)


CI_Bal.Accuracy_Test_Agnostic_IndGenes <- quantile(Bal.Accuracy_Test_Agnostic_IndGenes$Testing, c(0.025, 0.975))
CI_Bal.Accuracy_Test_Agnostic_IndGenes

# CIs for sensitivity
Sensitivity_Test_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Sensitivity", model_type == "Agnostic_IndGenes") %>%
  select(Testing)

CI_Sensitivity_Test_Agnostic_IndGenes <- quantile(Sensitivity_Test_Agnostic_IndGenes$Testing, c(0.025, 0.975))
CI_Sensitivity_Test_Agnostic_IndGenes

# CIs for Specificity
Specificity_Test_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "Specificity", model_type == "Agnostic_IndGenes") %>%
  select(Testing)

CI_Specificity_Test_Agnostic_IndGenes <- quantile(Specificity_Test_Agnostic_IndGenes$Testing, c(0.025, 0.975))
CI_Specificity_Test_Agnostic_IndGenes

# CIs for MCC
MCC_Test_Agnostic_IndGenes <- All_XGB_CV %>%
  filter(metric == "MCC", model_type == "Agnostic_IndGenes") %>%
  select(Testing)

CI_MCC_Test_Agnostic_IndGenes <- quantile(MCC_Test_Agnostic_IndGenes$Testing, c(0.025, 0.975))
CI_MCC_Test_Agnostic_IndGenes

#########################################
### Agnostic models (Pairs)

## Training data

## CIs for AUC
AUCs_Train_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "AUC", model_type == "Agnostic_Pairs") %>%
  select(Training)

CI_AUCs_Train_Agnostic_Pairs <- quantile(AUCs_Train_Agnostic_Pairs$Training, c(0.025, 0.975))
CI_AUCs_Train_Agnostic_Pairs

# CIs for Accuracy
Accuracy_Train_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Accuracy", model_type == "Agnostic_Pairs") %>%
  select(Training)

CI_Accuracy_Train_Agnostic_Pairs <- quantile(Accuracy_Train_Agnostic_Pairs$Training, c(0.025, 0.975))
CI_Accuracy_Train_Agnostic_Pairs

# CIs for Bal. Accuracy
Bal.Accuracy_Train_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Agnostic_Pairs") %>%
  select(Training)


CI_Bal.Accuracy_Train_Agnostic_Pairs <- quantile(Bal.Accuracy_Train_Agnostic_Pairs$Training, c(0.025, 0.975))
CI_Bal.Accuracy_Train_Agnostic_Pairs

# CIs for sensitivity
Sensitivity_Train_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Sensitivity", model_type == "Agnostic_Pairs") %>%
  select(Training)

CI_Sensitivity_Train_Agnostic_Pairs <- quantile(Sensitivity_Train_Agnostic_Pairs$Training, c(0.025, 0.975))
CI_Sensitivity_Train_Agnostic_Pairs

# CIs for Specificity
Specificity_Train_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Specificity", model_type == "Agnostic_Pairs") %>%
  select(Training)

CI_Specificity_Train_Agnostic_Pairs <- quantile(Specificity_Train_Agnostic_Pairs$Training, c(0.025, 0.975))
CI_Specificity_Train_Agnostic_Pairs

# CIs for MCC
MCC_Train_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "MCC", model_type == "Agnostic_Pairs") %>%
  select(Training)

CI_MCC_Train_Agnostic_Pairs <- quantile(MCC_Train_Agnostic_Pairs$Training, c(0.025, 0.975))
CI_MCC_Train_Agnostic_Pairs

####################
## Testing data

## CIs for AUC
AUCs_Test_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "AUC", model_type == "Agnostic_Pairs") %>%
  select(Testing)

CI_AUCs_Test_Agnostic_Pairs <- quantile(AUCs_Test_Agnostic_Pairs$Testing, c(0.025, 0.975))
CI_AUCs_Test_Agnostic_Pairs

# CIs for Accuracy
Accuracy_Test_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Accuracy", model_type == "Agnostic_Pairs") %>%
  select(Testing)

CI_Accuracy_Test_Agnostic_Pairs <- quantile(Accuracy_Test_Agnostic_Pairs$Testing, c(0.025, 0.975))
CI_Accuracy_Test_Agnostic_Pairs

# CIs for Bal. Accuracy
Bal.Accuracy_Test_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Agnostic_Pairs") %>%
  select(Testing)


CI_Bal.Accuracy_Test_Agnostic_Pairs <- quantile(Bal.Accuracy_Test_Agnostic_Pairs$Testing, c(0.025, 0.975))
CI_Bal.Accuracy_Test_Agnostic_Pairs

# CIs for sensitivity
Sensitivity_Test_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Sensitivity", model_type == "Agnostic_Pairs") %>%
  select(Testing)

CI_Sensitivity_Test_Agnostic_Pairs <- quantile(Sensitivity_Test_Agnostic_Pairs$Testing, c(0.025, 0.975))
CI_Sensitivity_Test_Agnostic_Pairs

# CIs for Specificity
Specificity_Test_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "Specificity", model_type == "Agnostic_Pairs") %>%
  select(Testing)

CI_Specificity_Test_Agnostic_Pairs <- quantile(Specificity_Test_Agnostic_Pairs$Testing, c(0.025, 0.975))
CI_Specificity_Test_Agnostic_Pairs

# CIs for MCC
MCC_Test_Agnostic_Pairs <- All_XGB_CV %>%
  filter(metric == "MCC", model_type == "Agnostic_Pairs") %>%
  select(Testing)

CI_MCC_Test_Agnostic_Pairs <- quantile(MCC_Test_Agnostic_Pairs$Testing, c(0.025, 0.975))
CI_MCC_Test_Agnostic_Pairs

##################
##############################
### Mechanistic models

## Training data

## CIs for AUC
AUCs_Train_Mechanistic <- All_XGB_CV %>%
  filter(metric == "AUC", model_type == "Mechanistic") %>%
  select(Training)

CI_AUCs_Train_Mechanistic <- quantile(AUCs_Train_Mechanistic$Training, c(0.025, 0.975))
CI_AUCs_Train_Mechanistic

# CIs for Accuracy
Accuracy_Train_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Accuracy", model_type == "Mechanistic") %>%
  select(Training)

CI_Accuracy_Train_Mechanistic <- quantile(Accuracy_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Accuracy_Train_Mechanistic

# CIs for Bal. Accuracy
Bal.Accuracy_Train_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Mechanistic") %>%
  select(Training)


CI_Bal.Accuracy_Train_Mechanistic <- quantile(Bal.Accuracy_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Bal.Accuracy_Train_Mechanistic

# CIs for sensitivity
Sensitivity_Train_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Sensitivity", model_type == "Mechanistic") %>%
  select(Training)

CI_Sensitivity_Train_Mechanistic <- quantile(Sensitivity_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Sensitivity_Train_Mechanistic

# CIs for Specificity
Specificity_Train_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Specificity", model_type == "Mechanistic") %>%
  select(Training)

CI_Specificity_Train_Mechanistic <- quantile(Specificity_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Specificity_Train_Mechanistic

# CIs for MCC
MCC_Train_Mechanistic <- All_XGB_CV %>%
  filter(metric == "MCC", model_type == "Mechanistic") %>%
  select(Training)

CI_MCC_Train_Mechanistic <- quantile(MCC_Train_Mechanistic$Training, c(0.025, 0.975))
CI_MCC_Train_Mechanistic

####################
## Testing data

## CIs for AUC
AUCs_Test_Mechanistic <- All_XGB_CV %>%
  filter(metric == "AUC", model_type == "Mechanistic") %>%
  select(Testing)

CI_AUCs_Test_Mechanistic <- quantile(AUCs_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_AUCs_Test_Mechanistic

# CIs for Accuracy
Accuracy_Test_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Accuracy", model_type == "Mechanistic") %>%
  select(Testing)

CI_Accuracy_Test_Mechanistic <- quantile(Accuracy_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Accuracy_Test_Mechanistic

# CIs for Bal. Accuracy
Bal.Accuracy_Test_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Mechanistic") %>%
  select(Testing)


CI_Bal.Accuracy_Test_Mechanistic <- quantile(Bal.Accuracy_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Bal.Accuracy_Test_Mechanistic

# CIs for sensitivity
Sensitivity_Test_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Sensitivity", model_type == "Mechanistic") %>%
  select(Testing)

CI_Sensitivity_Test_Mechanistic <- quantile(Sensitivity_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Sensitivity_Test_Mechanistic

# CIs for Specificity
Specificity_Test_Mechanistic <- All_XGB_CV %>%
  filter(metric == "Specificity", model_type == "Mechanistic") %>%
  select(Testing)

CI_Specificity_Test_Mechanistic <- quantile(Specificity_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Specificity_Test_Mechanistic

# CIs for MCC
MCC_Test_Mechanistic <- All_XGB_CV %>%
  filter(metric == "MCC", model_type == "Mechanistic") %>%
  select(Testing)

CI_MCC_Test_Mechanistic <- quantile(MCC_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_MCC_Test_Mechanistic


######################################
MeansAgnostic_IndGenes <- MeansAgnostic_IndGenes[-grep("SD", rownames(MeansAgnostic_IndGenes)), ]
MeansAgnostic_Pairs <- MeansAgnostic_Pairs[-grep("SD", rownames(MeansAgnostic_Pairs)), ]
MeansMechanistic <- MeansMechanistic[-grep("SD", rownames(MeansMechanistic)), ]

MeansAgnostic_IndGenes <- data.frame(MeansAgnostic_IndGenes)
MeansAgnostic_IndGenes$metric <- rownames(MeansAgnostic_IndGenes)
colnames(MeansAgnostic_IndGenes) <- c("Training", "Testing", "metric")

MeansAgnostic_Pairs <- data.frame(MeansAgnostic_Pairs)
MeansAgnostic_Pairs$metric <- rownames(MeansAgnostic_Pairs)
colnames(MeansAgnostic_Pairs) <- c("Training", "Testing", "metric")

MeansMechanistic <- data.frame(MeansMechanistic)
MeansMechanistic$metric <- rownames(MeansMechanistic)
colnames(MeansMechanistic) <- c("Training", "Testing", "metric")




MeansAgnostic_IndGenes <- MeansAgnostic_IndGenes %>%
  add_row(Training = CI_AUCs_Train_Agnostic_IndGenes[1], Testing = CI_AUCs_Test_Agnostic_IndGenes[1], metric = "AUC_CILow", .before = 2) %>%
  add_row(Training = CI_AUCs_Train_Agnostic_IndGenes[2], Testing = CI_AUCs_Test_Agnostic_IndGenes[2], metric = "AUC_CIHigh", .before = 3) %>%
  add_row(Training = CI_Accuracy_Train_Agnostic_IndGenes[1], Testing = CI_Accuracy_Test_Agnostic_IndGenes[1], metric = "Accuracy_CILow", .before = 5) %>%
  add_row(Training = CI_Accuracy_Train_Agnostic_IndGenes[2], Testing = CI_Accuracy_Test_Agnostic_IndGenes[2], metric = "Accuracy_CIHigh", .before = 6) %>%
  add_row(Training = CI_Bal.Accuracy_Train_Agnostic_IndGenes[1], Testing = CI_Bal.Accuracy_Test_Agnostic_IndGenes[1], metric = "BalAccuracy_CILow", .before = 8) %>%
  add_row(Training = CI_Bal.Accuracy_Train_Agnostic_IndGenes[2], Testing = CI_Bal.Accuracy_Test_Agnostic_IndGenes[2], metric = "BalAccuracy_CIHigh", .before = 9) %>%
  add_row(Training = CI_Sensitivity_Train_Agnostic_IndGenes[1], Testing = CI_Sensitivity_Test_Agnostic_IndGenes[1], metric = "Sensitivity_CILow", .before = 11) %>%
  add_row(Training = CI_Sensitivity_Train_Agnostic_IndGenes[2], Testing = CI_Sensitivity_Test_Agnostic_IndGenes[2], metric = "Sensitvity_CIHigh", .before = 12) %>%
  add_row(Training = CI_Specificity_Train_Agnostic_IndGenes[1], Testing = CI_Specificity_Test_Agnostic_IndGenes[1], metric = "Specificity_CILow", .before = 14) %>%
  add_row(Training = CI_Specificity_Train_Agnostic_IndGenes[2], Testing = CI_Specificity_Test_Agnostic_IndGenes[2], metric = "Specificity_CIHigh", .before = 15) %>%
  add_row(Training = CI_MCC_Train_Agnostic_IndGenes[1], Testing = CI_MCC_Test_Agnostic_IndGenes[1], metric = "MCC_CILow", .before = 17) %>%
  add_row(Training = CI_MCC_Train_Agnostic_IndGenes[2], Testing = CI_MCC_Test_Agnostic_IndGenes[2], metric = "MCC_CIHigh", .before = 18)

rownames(MeansAgnostic_IndGenes) <- MeansAgnostic_IndGenes$metric  
MeansAgnostic_IndGenes$model_type <- rep("Agnostic_IndGenes", 18)
MeansAgnostic_IndGenes$Algorithm <- rep("XGB", 18)
 
MeansAgnostic_Pairs <- MeansAgnostic_Pairs %>%
  add_row(Training = CI_AUCs_Train_Agnostic_Pairs[1], Testing = CI_AUCs_Test_Agnostic_Pairs[1], metric = "AUC_CILow", .before = 2) %>%
  add_row(Training = CI_AUCs_Train_Agnostic_Pairs[2], Testing = CI_AUCs_Test_Agnostic_Pairs[2], metric = "AUC_CIHigh", .before = 3) %>%
  add_row(Training = CI_Accuracy_Train_Agnostic_Pairs[1], Testing = CI_Accuracy_Test_Agnostic_Pairs[1], metric = "Accuracy_CILow", .before = 5) %>%
  add_row(Training = CI_Accuracy_Train_Agnostic_Pairs[2], Testing = CI_Accuracy_Test_Agnostic_Pairs[2], metric = "Accuracy_CIHigh", .before = 6) %>%
  add_row(Training = CI_Bal.Accuracy_Train_Agnostic_Pairs[1], Testing = CI_Bal.Accuracy_Test_Agnostic_Pairs[1], metric = "BalAccuracy_CILow", .before = 8) %>%
  add_row(Training = CI_Bal.Accuracy_Train_Agnostic_Pairs[2], Testing = CI_Bal.Accuracy_Test_Agnostic_Pairs[2], metric = "BalAccuracy_CIHigh", .before = 9) %>%
  add_row(Training = CI_Sensitivity_Train_Agnostic_Pairs[1], Testing = CI_Sensitivity_Test_Agnostic_Pairs[1], metric = "Sensitivity_CILow", .before = 11) %>%
  add_row(Training = CI_Sensitivity_Train_Agnostic_Pairs[2], Testing = CI_Sensitivity_Test_Agnostic_Pairs[2], metric = "Sensitvity_CIHigh", .before = 12) %>%
  add_row(Training = CI_Specificity_Train_Agnostic_Pairs[1], Testing = CI_Specificity_Test_Agnostic_Pairs[1], metric = "Specificity_CILow", .before = 14) %>%
  add_row(Training = CI_Specificity_Train_Agnostic_Pairs[2], Testing = CI_Specificity_Test_Agnostic_Pairs[2], metric = "Specificity_CIHigh", .before = 15) %>%
  add_row(Training = CI_MCC_Train_Agnostic_Pairs[1], Testing = CI_MCC_Test_Agnostic_Pairs[1], metric = "MCC_CILow", .before = 17) %>%
  add_row(Training = CI_MCC_Train_Agnostic_Pairs[2], Testing = CI_MCC_Test_Agnostic_Pairs[2], metric = "MCC_CIHigh", .before = 18)

rownames(MeansAgnostic_Pairs) <- MeansAgnostic_Pairs$metric  
MeansAgnostic_Pairs$model_type <- rep("Agnostic_Pairs", 18)
MeansAgnostic_Pairs$Algorithm <- rep("XGB", 18)


MeansMechanistic <- MeansMechanistic %>%
  add_row(Training = CI_AUCs_Train_Mechanistic[1], Testing = CI_AUCs_Test_Mechanistic[1], metric = "AUC_CILow", .before = 2) %>%
  add_row(Training = CI_AUCs_Train_Mechanistic[2], Testing = CI_AUCs_Test_Mechanistic[2], metric = "AUC_CIHigh", .before = 3) %>%
  add_row(Training = CI_Accuracy_Train_Mechanistic[1], Testing = CI_Accuracy_Test_Mechanistic[1], metric = "Accuracy_CILow", .before = 5) %>%
  add_row(Training = CI_Accuracy_Train_Mechanistic[2], Testing = CI_Accuracy_Test_Mechanistic[2], metric = "Accuracy_CIHigh", .before = 6) %>%
  add_row(Training = CI_Bal.Accuracy_Train_Mechanistic[1], Testing = CI_Bal.Accuracy_Test_Mechanistic[1], metric = "BalAccuracy_CILow", .before = 8) %>%
  add_row(Training = CI_Bal.Accuracy_Train_Mechanistic[2], Testing = CI_Bal.Accuracy_Test_Mechanistic[2], metric = "BalAccuracy_CIHigh", .before = 9) %>%
  add_row(Training = CI_Sensitivity_Train_Mechanistic[1], Testing = CI_Sensitivity_Test_Mechanistic[1], metric = "Sensitivity_CILow", .before = 11) %>%
  add_row(Training = CI_Sensitivity_Train_Mechanistic[2], Testing = CI_Sensitivity_Test_Mechanistic[2], metric = "Sensitvity_CIHigh", .before = 12) %>%
  add_row(Training = CI_Specificity_Train_Mechanistic[1], Testing = CI_Specificity_Test_Mechanistic[1], metric = "Specificity_CILow", .before = 14) %>%
  add_row(Training = CI_Specificity_Train_Mechanistic[2], Testing = CI_Specificity_Test_Mechanistic[2], metric = "Specificity_CIHigh", .before = 15) %>%
  add_row(Training = CI_MCC_Train_Mechanistic[1], Testing = CI_MCC_Test_Mechanistic[1], metric = "MCC_CILow", .before = 17) %>%
  add_row(Training = CI_MCC_Train_Mechanistic[2], Testing = CI_MCC_Test_Mechanistic[2], metric = "MCC_CIHigh", .before = 18)

rownames(MeansMechanistic) <- MeansMechanistic$metric  
MeansMechanistic$model_type <- rep("Mechanistic", 18)
MeansMechanistic$Algorithm <- rep("XGB", 18)
 
 
XGB_CV_AllMetrics <- rbind(MeansAgnostic_IndGenes, MeansAgnostic_Pairs, MeansMechanistic)

XGB_CV_AllMetrics <- XGB_CV_AllMetrics %>%
  pivot_longer(cols = c(1,2), names_to = c("data_type"))

XGB_CV_AllMetrics$type <- XGB_CV_AllMetrics$metric
XGB_CV_AllMetrics$metric <- gsub("Mean_", "", XGB_CV_AllMetrics$metric)
#XGB_CV_AllMetrics$metric <- gsub("_Agnostic", "", XGB_CV_AllMetrics$metric)
#XGB_CV_AllMetrics$metric <- gsub("_Mechanistic", "", XGB_CV_AllMetrics$metric)
XGB_CV_AllMetrics$metric <- gsub("_CILow", "", XGB_CV_AllMetrics$metric)
XGB_CV_AllMetrics$metric <- gsub("_CIHigh", "", XGB_CV_AllMetrics$metric)

XGB_CV_AllMetrics$type[grep("Mean", XGB_CV_AllMetrics$type)] <- "Mean"
XGB_CV_AllMetrics$type[grep("CILow", XGB_CV_AllMetrics$type)] <- "CILow"
XGB_CV_AllMetrics$type[grep("CIHigh", XGB_CV_AllMetrics$type)] <- "CIHigh"

XGB_CV_AllMetrics$metric[XGB_CV_AllMetrics$metric == "Sensitvity"] <- "Sensitivity"

XGB_CV_AllMetrics <- XGB_CV_AllMetrics %>% pivot_wider(names_from = "type", values_from = c("value"))
 
XGB_CV_AllMetrics$data_type <- factor(XGB_CV_AllMetrics$data_type, levels = c("Training", "Testing"))  
XGB_CV_AllMetrics$model_type <- factor(XGB_CV_AllMetrics$model_type, levels = c("Mechanistic", "Agnostic_Pairs", "Agnostic_IndGenes"))
XGB_CV_AllMetrics$metric <- factor(XGB_CV_AllMetrics$metric, levels = c("AUC", "Accuracy", "BalAccuracy", "Sensitivity", "Specificity", "MCC"))
# 
# ## Plot the figure
png(filename = "./Figs/XGB/XGBCV_ALlMetrics.png", width = 3000, height = 2000, res = 300)
ggplot(data = XGB_CV_AllMetrics, aes(x = Mean, y = model_type, shape =  model_type, fill = model_type)) +
  geom_linerange(aes(xmin = CILow, xmax = CIHigh)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = 3) +
  facet_grid(metric ~ data_type, scales = "free") +
  scale_x_continuous(lim = c(-0.5, 1)) +
  scale_shape_manual(values = c(21, 22, 23)) +
  scale_fill_manual(values = c("white", "grey20", "grey50")) +
  labs(y = NULL,
       x = "Mean (95% Confidence Interval)") +
  theme(legend.position = "none")

dev.off()
 
## Save the Dataframe 
save(XGB_CV_AllMetrics, file = "./Objs/XGB/XGB_CV_AllMetrics.rda")
# 
