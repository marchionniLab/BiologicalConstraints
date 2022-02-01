### Forest plot comparing the agnostic and mechanistic KTSP in the CV part using all metrics

rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)

## Agnostic SVM Average performance

# Load the agnostic models (IndGenes)
AgnosticEMTAB_Out_IndGenes <- load("./Objs/SVM/EMTAB_Out_SVM_IndvGenes_AgnosticPerformance.rda")
AgnosticGSE13507_Out_IndGenes <- load("./Objs/SVM/GSE13507_Out_SVM_IndvGenes_AgnosticPerformance.rda")
AgnosticGSE32894_Out_IndGenes <- load("./Objs/SVM/GSE32894_Out_SVM_IndvGenes_AgnosticPerformance.rda")
AgnosticGSE57813_Out_IndGenes <- load("./Objs/SVM/GSE57813_Out_SVM_IndvGenes_AgnosticPerformance.rda")
AgnosticPMID_Out_IndGenes <- load("./Objs/SVM/PMID_Out_SVM_IndvGenes_AgnosticPerformance.rda")

# Load the agnostic models (Pairs)
AgnosticEMTAB_Out_Pairs <- load("./Objs/SVM/EMTAB_Out_SVM_AgnosticPerformance.rda")
AgnosticGSE13507_Out_Pairs <- load("./Objs/SVM/GSE13507_Out_SVM_AgnosticPerformance.rda")
AgnosticGSE32894_Out_Pairs <- load("./Objs/SVM/GSE32894_Out_SVM_AgnosticPerformance.rda")
AgnosticGSE57813_Out_Pairs <- load("./Objs/SVM/GSE57813_Out_SVM_AgnosticPerformance.rda")
AgnosticPMID_Out_Pairs <- load("./Objs/SVM/PMID_Out_SVM_AgnosticPerformance.rda")

# Load the Mech models
MechEMTAB_Out <- load("./Objs/SVM/EMTAB_Out_SVM_MechPerformance.rda")
MechGSE13507_Out <- load("./Objs/SVM/GSE13507_Out_SVM_MechPerformance.rda")
MechGSE32894_Out <- load("./Objs/SVM/GSE32894_Out_SVM_MechPerformance.rda")
MechGSE57813_Out <- load("./Objs/SVM/GSE57813_Out_SVM_MechPerformance.rda")
MechPMID_Out <- load("./Objs/SVM/PMID_Out_SVM_MechPerformance.rda")




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
EMTAB_Out_SVM_IndvGenes_AgnosticPerformance$out_study <- rep("EMTAB", nrow(EMTAB_Out_SVM_IndvGenes_AgnosticPerformance))
EMTAB_Out_SVM_AgnosticPerformance$out_study <- rep("EMTAB", nrow(EMTAB_Out_SVM_AgnosticPerformance))
EMTAB_Out_SVM_MechPerformance$out_study <- rep("EMTAB", nrow(EMTAB_Out_SVM_MechPerformance))

EMTAB_Out_SVM_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(EMTAB_Out_SVM_IndvGenes_AgnosticPerformance))
EMTAB_Out_SVM_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(EMTAB_Out_SVM_IndvGenes_AgnosticPerformance))
EMTAB_Out_SVM_MechPerformance$model_type <- rep("Mechanistic", nrow(EMTAB_Out_SVM_MechPerformance))

# GSE13507_out
GSE13507_Out_SVM_IndvGenes_AgnosticPerformance$out_study <- rep("GSE13507", nrow(GSE13507_Out_SVM_IndvGenes_AgnosticPerformance))
GSE13507_Out_SVM_AgnosticPerformance$out_study <- rep("GSE13507", nrow(GSE13507_Out_SVM_AgnosticPerformance))
GSE13507_Out_SVM_MechPerformance$out_study <- rep("GSE13507", nrow(GSE13507_Out_SVM_MechPerformance))

GSE13507_Out_SVM_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(GSE13507_Out_SVM_IndvGenes_AgnosticPerformance))
GSE13507_Out_SVM_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(GSE13507_Out_SVM_AgnosticPerformance))
GSE13507_Out_SVM_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE13507_Out_SVM_MechPerformance))

# GSE32894_out
GSE32894_Out_SVM_IndvGenes_AgnosticPerformance$out_study <- rep("GSE32894", nrow(GSE32894_Out_SVM_IndvGenes_AgnosticPerformance))
GSE32894_Out_SVM_AgnosticPerformance$out_study <- rep("GSE32894", nrow(GSE32894_Out_SVM_AgnosticPerformance))
GSE32894_Out_SVM_MechPerformance$out_study <- rep("GSE32894", nrow(GSE32894_Out_SVM_MechPerformance))

GSE32894_Out_SVM_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(GSE32894_Out_SVM_IndvGenes_AgnosticPerformance))
GSE32894_Out_SVM_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(GSE32894_Out_SVM_AgnosticPerformance))
GSE32894_Out_SVM_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE32894_Out_SVM_MechPerformance))

# GSE57813_out
GSE57813_Out_SVM_IndvGenes_AgnosticPerformance$out_study <- rep("GSE57813", nrow(GSE57813_Out_SVM_IndvGenes_AgnosticPerformance))
GSE57813_Out_SVM_AgnosticPerformance$out_study <- rep("GSE57813", nrow(GSE57813_Out_SVM_AgnosticPerformance))
GSE57813_Out_SVM_MechPerformance$out_study <- rep("GSE57813", nrow(GSE57813_Out_SVM_MechPerformance))

GSE57813_Out_SVM_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(GSE57813_Out_SVM_IndvGenes_AgnosticPerformance))
GSE57813_Out_SVM_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(GSE57813_Out_SVM_AgnosticPerformance))
GSE57813_Out_SVM_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE57813_Out_SVM_MechPerformance))

# PMID_out
PMID_Out_SVM_IndvGenes_AgnosticPerformance$out_study <- rep("PMID", nrow(PMID_Out_SVM_IndvGenes_AgnosticPerformance))
PMID_Out_SVM_AgnosticPerformance$out_study <- rep("PMID", nrow(PMID_Out_SVM_AgnosticPerformance))
PMID_Out_SVM_MechPerformance$out_study <- rep("PMID", nrow(PMID_Out_SVM_MechPerformance))

PMID_Out_SVM_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_IndGenes", nrow(PMID_Out_SVM_IndvGenes_AgnosticPerformance))
PMID_Out_SVM_AgnosticPerformance$model_type <- rep("Agnostic_Pairs", nrow(PMID_Out_SVM_AgnosticPerformance))
PMID_Out_SVM_MechPerformance$model_type <- rep("Mechanistic", nrow(PMID_Out_SVM_MechPerformance))


All_SVM_CV <- rbind(EMTAB_Out_SVM_IndvGenes_AgnosticPerformance, EMTAB_Out_SVM_AgnosticPerformance, EMTAB_Out_SVM_MechPerformance,
                    GSE13507_Out_SVM_IndvGenes_AgnosticPerformance, GSE13507_Out_SVM_AgnosticPerformance, GSE13507_Out_SVM_MechPerformance,
                    GSE32894_Out_SVM_IndvGenes_AgnosticPerformance, GSE32894_Out_SVM_AgnosticPerformance, GSE32894_Out_SVM_MechPerformance, 
                    GSE57813_Out_SVM_IndvGenes_AgnosticPerformance, GSE57813_Out_SVM_AgnosticPerformance, GSE57813_Out_SVM_MechPerformance, 
                    PMID_Out_SVM_IndvGenes_AgnosticPerformance, PMID_Out_SVM_AgnosticPerformance, PMID_Out_SVM_MechPerformance)

All_SVM_CV$metric <- rownames(All_SVM_CV)
All_SVM_CV$metric <- gsub('[[:digit:]]+', '', All_SVM_CV$metric )


##########################################################
##########################################################
## Calculate the average statistics (Across the 5 repeats)

#####
# average stats in the agnostic models (Ind Genes)

## AUC
# Mean AUC in both the training and testing data
MeanAUC_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                             filter(model_type == "Agnostic_IndGenes") %>%
                             filter(metric == "AUC") %>%
                             select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                            filter(model_type == "Agnostic_IndGenes") %>%
                            filter(metric == "AUC") %>%
                            select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                  filter(model_type == "Agnostic_IndGenes") %>%
                                  filter(metric == "Accuracy") %>%
                                  select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                 filter(model_type == "Agnostic_IndGenes") %>%
                                 filter(metric == "Accuracy") %>%
                                 select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                     filter(model_type == "Agnostic_IndGenes") %>%
                                     filter(metric == "Bal.Accuracy") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                    filter(model_type == "Agnostic_IndGenes") %>%
                                    filter(metric == "Bal.Accuracy") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                     filter(model_type == "Agnostic_IndGenes") %>%
                                     filter(metric == "Sensitivity") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                    filter(model_type == "Agnostic_IndGenes") %>%
                                    filter(metric == "Sensitivity") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                     filter(model_type == "Agnostic_IndGenes") %>%
                                     filter(metric == "Specificity") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                                    filter(model_type == "Agnostic_IndGenes") %>%
                                    filter(metric == "Specificity") %>%
                                    select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
                             filter(model_type == "Agnostic_IndGenes") %>%
                             filter(metric == "MCC") %>%
                             select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_IndGenes <-  apply(All_SVM_CV %>%
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

################################################
#####
# average stats in the agnostic models (Pairs)
## AUC
# Mean AUC in both the training and testing data
MeanAUC_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                      filter(model_type == "Agnostic_Pairs") %>%
                                      filter(metric == "AUC") %>%
                                      select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                     filter(model_type == "Agnostic_Pairs") %>%
                                     filter(metric == "AUC") %>%
                                     select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                           filter(model_type == "Agnostic_Pairs") %>%
                                           filter(metric == "Accuracy") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                          filter(model_type == "Agnostic_Pairs") %>%
                                          filter(metric == "Accuracy") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                              filter(model_type == "Agnostic_Pairs") %>%
                                              filter(metric == "Bal.Accuracy") %>%
                                              select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                             filter(model_type == "Agnostic_Pairs") %>%
                                             filter(metric == "Bal.Accuracy") %>%
                                             select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                              filter(model_type == "Agnostic_Pairs") %>%
                                              filter(metric == "Sensitivity") %>%
                                              select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                             filter(model_type == "Agnostic_Pairs") %>%
                                             filter(metric == "Sensitivity") %>%
                                             select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                              filter(model_type == "Agnostic_Pairs") %>%
                                              filter(metric == "Specificity") %>%
                                              select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                             filter(model_type == "Agnostic_Pairs") %>%
                                             filter(metric == "Specificity") %>%
                                             select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_Pairs <-  apply(All_SVM_CV %>%
                                      filter(model_type == "Agnostic_Pairs") %>%
                                      filter(metric == "MCC") %>%
                                      select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_Pairs <-  apply(All_SVM_CV %>%
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
MeanAUC_Mechanistic <-  apply(All_SVM_CV %>%
                                filter(model_type == "Mechanistic") %>%
                                filter(metric == "AUC") %>%
                                select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Mechanistic <-  apply(All_SVM_CV %>%
                               filter(model_type == "Mechanistic") %>%
                               filter(metric == "AUC") %>%
                               select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Mechanistic <-  apply(All_SVM_CV %>%
                                     filter(model_type == "Mechanistic") %>%
                                     filter(metric == "Accuracy") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Mechanistic <-  apply(All_SVM_CV %>%
                                    filter(model_type == "Mechanistic") %>%
                                    filter(metric == "Accuracy") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Mechanistic <-  apply(All_SVM_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Bal.Accuracy") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Mechanistic <-  apply(All_SVM_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Bal.Accuracy") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Mechanistic <-  apply(All_SVM_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Sensitivity") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Mechanistic <-  apply(All_SVM_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Sensitivity") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Mechanistic <-  apply(All_SVM_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Specificity") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Mechanistic <-  apply(All_SVM_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Specificity") %>%
                                       select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Mechanistic <-  apply(All_SVM_CV %>%
                                filter(model_type == "Mechanistic") %>%
                                filter(metric == "MCC") %>%
                                select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Mechanistic <-  apply(All_SVM_CV %>%
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

########
## Bind mean agnostic and mean mechanistic together
colnames(MeansAgnostic_IndGenes) <- paste("Agnostic_IndGenes", colnames(MeansAgnostic_IndGenes), sep = "_")
colnames(MeansAgnostic_Pairs) <- paste("Agnostic_Pairs", colnames(MeansAgnostic_Pairs), sep = "_")
colnames(MeansMechanistic) <- paste("Mechanistic", colnames(MeansMechanistic), sep = "_")

# Modify the rownames
rownames(MeansAgnostic_IndGenes) <- c("Mean_AUC", "SD_AUC", "Mean_Accuracy", "SD_Accuracy", "Mean_BalAccuracy", "SD_BalAccuracy", "Mean_Sensitivity", "SD_Sensitivity", "Mean_Specificity", "SD_Specificity", "Mean_MCC", "SD_MCC")
rownames(MeansAgnostic_Pairs) <- c("Mean_AUC", "SD_AUC", "Mean_Accuracy", "SD_Accuracy", "Mean_BalAccuracy", "SD_BalAccuracy", "Mean_Sensitivity", "SD_Sensitivity", "Mean_Specificity", "SD_Specificity", "Mean_MCC", "SD_MCC")
rownames(MeansMechanistic) <- c("Mean_AUC", "SD_AUC", "Mean_Accuracy", "SD_Accuracy", "Mean_BalAccuracy", "SD_BalAccuracy", "Mean_Sensitivity", "SD_Sensitivity", "Mean_Specificity", "SD_Specificity", "Mean_MCC", "SD_MCC")

SVM_CrossStudyValidation <- cbind(MeansAgnostic_IndGenes, MeansAgnostic_Pairs, MeansMechanistic)
save(SVM_CrossStudyValidation, file = "./Objs/SVM/SVM_CrossStudyValidation.rda")

load("./Objs/SVM/SVM_CrossStudyValidation.rda")

################
# Functions to calculate the 95% CI from the mean and SD

# CalculateCIHigh <- function(x, y){Z <- x + 1.96 * (y/sqrt(5))
# return(Z)
# }
# 
# 
# CalculateCILow <- function(x, y){Z <- x - 1.96 * (y/sqrt(5))
# return(Z)
# }
# 
# 
# MeansAgnostic <- data.frame(MeansAgnostic)
# MeansAgnostic$metric <- rownames(MeansAgnostic)
# 
# MeansMechanistic <- data.frame(MeansMechanistic)
# MeansMechanistic$metric <- rownames(MeansMechanistic)

############################3
## Use the mean and SD to calculate the 95% CI 

### Agnostic models

## Training data

## CIs for AUC
AUCs_Train_Agnostic <- All_SVM_CV %>%
  filter(metric == "AUC", model_type == "Agnostic") %>%
  select(Training)

CI_AUCs_Train_Agnostic <- quantile(AUCs_Train_Agnostic$Training, c(0.025, 0.975))
CI_AUCs_Train_Agnostic

# CIs for Accuracy
Accuracy_Train_Agnostic <- All_SVM_CV %>%
  filter(metric == "Accuracy", model_type == "Agnostic") %>%
  select(Training)

CI_Accuracy_Train_Agnostic <- quantile(Accuracy_Train_Agnostic$Training, c(0.025, 0.975))
CI_Accuracy_Train_Agnostic

# CIs for Bal. Accuracy
Bal.Accuracy_Train_Agnostic <- All_SVM_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Agnostic") %>%
  select(Training)


CI_Bal.Accuracy_Train_Agnostic <- quantile(Bal.Accuracy_Train_Agnostic$Training, c(0.025, 0.975))
CI_Bal.Accuracy_Train_Agnostic

# CIs for sensitivity
Sensitivity_Train_Agnostic <- All_SVM_CV %>%
  filter(metric == "Sensitivity", model_type == "Agnostic") %>%
  select(Training)

CI_Sensitivity_Train_Agnostic <- quantile(Sensitivity_Train_Agnostic$Training, c(0.025, 0.975))
CI_Sensitivity_Train_Agnostic

# CIs for Specificity
Specificity_Train_Agnostic <- All_SVM_CV %>%
  filter(metric == "Specificity", model_type == "Agnostic") %>%
  select(Training)

CI_Specificity_Train_Agnostic <- quantile(Specificity_Train_Agnostic$Training, c(0.025, 0.975))
CI_Specificity_Train_Agnostic

# CIs for MCC
MCC_Train_Agnostic <- All_SVM_CV %>%
  filter(metric == "MCC", model_type == "Agnostic") %>%
  select(Training)

CI_MCC_Train_Agnostic <- quantile(MCC_Train_Agnostic$Training, c(0.025, 0.975))
CI_MCC_Train_Agnostic

####################
## Testing data

## CIs for AUC
AUCs_Test_Agnostic <- All_SVM_CV %>%
  filter(metric == "AUC", model_type == "Agnostic") %>%
  select(Testing)

CI_AUCs_Test_Agnostic <- quantile(AUCs_Test_Agnostic$Testing, c(0.025, 0.975))
CI_AUCs_Test_Agnostic

# CIs for Accuracy
Accuracy_Test_Agnostic <- All_SVM_CV %>%
  filter(metric == "Accuracy", model_type == "Agnostic") %>%
  select(Testing)

CI_Accuracy_Test_Agnostic <- quantile(Accuracy_Test_Agnostic$Testing, c(0.025, 0.975))
CI_Accuracy_Test_Agnostic

# CIs for Bal. Accuracy
Bal.Accuracy_Test_Agnostic <- All_SVM_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Agnostic") %>%
  select(Testing)


CI_Bal.Accuracy_Test_Agnostic <- quantile(Bal.Accuracy_Test_Agnostic$Testing, c(0.025, 0.975))
CI_Bal.Accuracy_Test_Agnostic

# CIs for sensitivity
Sensitivity_Test_Agnostic <- All_SVM_CV %>%
  filter(metric == "Sensitivity", model_type == "Agnostic") %>%
  select(Testing)

CI_Sensitivity_Test_Agnostic <- quantile(Sensitivity_Test_Agnostic$Testing, c(0.025, 0.975))
CI_Sensitivity_Test_Agnostic

# CIs for Specificity
Specificity_Test_Agnostic <- All_SVM_CV %>%
  filter(metric == "Specificity", model_type == "Agnostic") %>%
  select(Testing)

CI_Specificity_Test_Agnostic <- quantile(Specificity_Test_Agnostic$Testing, c(0.025, 0.975))
CI_Specificity_Test_Agnostic

# CIs for MCC
MCC_Test_Agnostic <- All_SVM_CV %>%
  filter(metric == "MCC", model_type == "Agnostic") %>%
  select(Testing)

CI_MCC_Test_Agnostic <- quantile(MCC_Test_Agnostic$Testing, c(0.025, 0.975))
CI_MCC_Test_Agnostic

##############################
### Mechanistic models

## Training data

## CIs for AUC
AUCs_Train_Mechanistic <- All_SVM_CV %>%
  filter(metric == "AUC", model_type == "Mechanistic") %>%
  select(Training)

CI_AUCs_Train_Mechanistic <- quantile(AUCs_Train_Mechanistic$Training, c(0.025, 0.975))
CI_AUCs_Train_Mechanistic

# CIs for Accuracy
Accuracy_Train_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Accuracy", model_type == "Mechanistic") %>%
  select(Training)

CI_Accuracy_Train_Mechanistic <- quantile(Accuracy_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Accuracy_Train_Mechanistic

# CIs for Bal. Accuracy
Bal.Accuracy_Train_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Mechanistic") %>%
  select(Training)


CI_Bal.Accuracy_Train_Mechanistic <- quantile(Bal.Accuracy_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Bal.Accuracy_Train_Mechanistic

# CIs for sensitivity
Sensitivity_Train_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Sensitivity", model_type == "Mechanistic") %>%
  select(Training)

CI_Sensitivity_Train_Mechanistic <- quantile(Sensitivity_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Sensitivity_Train_Mechanistic

# CIs for Specificity
Specificity_Train_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Specificity", model_type == "Mechanistic") %>%
  select(Training)

CI_Specificity_Train_Mechanistic <- quantile(Specificity_Train_Mechanistic$Training, c(0.025, 0.975))
CI_Specificity_Train_Mechanistic

# CIs for MCC
MCC_Train_Mechanistic <- All_SVM_CV %>%
  filter(metric == "MCC", model_type == "Mechanistic") %>%
  select(Training)

CI_MCC_Train_Mechanistic <- quantile(MCC_Train_Mechanistic$Training, c(0.025, 0.975))
CI_MCC_Train_Mechanistic

####################
## Testing data

## CIs for AUC
AUCs_Test_Mechanistic <- All_SVM_CV %>%
  filter(metric == "AUC", model_type == "Mechanistic") %>%
  select(Testing)

CI_AUCs_Test_Mechanistic <- quantile(AUCs_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_AUCs_Test_Mechanistic

# CIs for Accuracy
Accuracy_Test_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Accuracy", model_type == "Mechanistic") %>%
  select(Testing)

CI_Accuracy_Test_Mechanistic <- quantile(Accuracy_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Accuracy_Test_Mechanistic

# CIs for Bal. Accuracy
Bal.Accuracy_Test_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Bal.Accuracy", model_type == "Mechanistic") %>%
  select(Testing)


CI_Bal.Accuracy_Test_Mechanistic <- quantile(Bal.Accuracy_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Bal.Accuracy_Test_Mechanistic

# CIs for sensitivity
Sensitivity_Test_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Sensitivity", model_type == "Mechanistic") %>%
  select(Testing)

CI_Sensitivity_Test_Mechanistic <- quantile(Sensitivity_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Sensitivity_Test_Mechanistic

# CIs for Specificity
Specificity_Test_Mechanistic <- All_SVM_CV %>%
  filter(metric == "Specificity", model_type == "Mechanistic") %>%
  select(Testing)

CI_Specificity_Test_Mechanistic <- quantile(Specificity_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_Specificity_Test_Mechanistic

# CIs for MCC
MCC_Test_Mechanistic <- All_SVM_CV %>%
  filter(metric == "MCC", model_type == "Mechanistic") %>%
  select(Testing)

CI_MCC_Test_Mechanistic <- quantile(MCC_Test_Mechanistic$Testing, c(0.025, 0.975))
CI_MCC_Test_Mechanistic


######################################
# MeansAgnostic <- MeansAgnostic[-grep("SD", rownames(MeansAgnostic)), ]
# MeansMechanistic <- MeansMechanistic[-grep("SD", rownames(MeansMechanistic)), ]
# 
# MeansAgnostic <- MeansAgnostic %>%
#   add_row(Training = Agnostic_Train_AUC_CI_Low, Testing = Agnostic_Test_AUC_CI_Low, metric = "AUC_CILow", .before = 2) %>%
#   add_row(Training = Agnostic_Train_AUC_CI_High, Testing = Agnostic_Test_AUC_CI_High, metric = "AUC_CIHigh", .before = 3) %>%
#   add_row(Training = Agnostic_Train_Accuracy_CI_Low, Testing = Agnostic_Test_Accuracy_CI_Low, metric = "Accuracy_CILow", .before = 5) %>%
#   add_row(Training = Agnostic_Train_Accuracy_CI_High, Testing = Agnostic_Test_Accuracy_CI_High, metric = "Accuracy_CIHigh", .before = 6) %>%
#   add_row(Training = Agnostic_Train_BalAccuracy_CI_Low, Testing = Agnostic_Test_BalAccuracy_CI_Low, metric = "BalAccuracy_CILow", .before = 8) %>%
#   add_row(Training = Agnostic_Train_BalAccuracy_CI_High, Testing = Agnostic_Test_BalAccuracy_CI_High, metric = "BalAccuracy_CIHigh", .before = 9) %>%
#   add_row(Training = Agnostic_Train_Sensitivity_CI_Low, Testing = Agnostic_Test_Sensitivity_CI_Low, metric = "Sensitivity_CILow", .before = 11) %>%
#   add_row(Training = Agnostic_Train_Sensitivity_CI_High, Testing = Agnostic_Test_Sensitivity_CI_High, metric = "Sensitvity_CIHigh", .before = 12) %>%
#   add_row(Training = Agnostic_Train_Specificity_CI_Low, Testing = Agnostic_Test_Specificity_CI_Low, metric = "Specificity_CILow", .before = 14) %>%
#   add_row(Training = Agnostic_Train_Specificity_CI_High, Testing = Agnostic_Test_Specificity_CI_High, metric = "Specificity_CIHigh", .before = 15) %>%
#   add_row(Training = Agnostic_Train_MCC_CI_Low, Testing = Agnostic_Test_MCC_CI_Low, metric = "MCC_CILow", .before = 17) %>%
#   add_row(Training = Agnostic_Train_MCC_CI_High, Testing = Agnostic_Test_MCC_CI_High, metric = "MCC_CIHigh", .before = 18)
# 
# rownames(MeansAgnostic) <- MeansAgnostic$metric  
# MeansAgnostic$model_type <- rep("Agnostic", 18)
# MeansAgnostic$Algorithm <- rep("SVM", 18)
# 
# MeansMechanistic <- MeansMechanistic %>%
#   add_row(Training = Mechanistic_Train_AUC_CI_Low, Testing = Mechanistic_Test_AUC_CI_Low, metric = "AUC_CILow", .before = 2) %>%
#   add_row(Training = Mechanistic_Train_AUC_CI_High, Testing = Mechanistic_Test_AUC_CI_High, metric = "AUC_CIHigh", .before = 3) %>%
#   add_row(Training = Mechanistic_Train_Accuracy_CI_Low, Testing = Mechanistic_Test_Accuracy_CI_Low, metric = "Accuracy_CILow", .before = 5) %>%
#   add_row(Training = Mechanistic_Train_Accuracy_CI_High, Testing = Mechanistic_Test_Accuracy_CI_High, metric = "Accuracy_CIHigh", .before = 6) %>%
#   add_row(Training = Mechanistic_Train_BalAccuracy_CI_Low, Testing = Mechanistic_Test_BalAccuracy_CI_Low, metric = "BalAccuracy_CILow", .before = 8) %>%
#   add_row(Training = Mechanistic_Train_BalAccuracy_CI_High, Testing = Mechanistic_Test_BalAccuracy_CI_High, metric = "BalAccuracy_CIHigh", .before = 9) %>%
#   add_row(Training = Mechanistic_Train_Sensitivity_CI_Low, Testing = Mechanistic_Test_Sensitivity_CI_Low, metric = "Sensitivity_CILow", .before = 11) %>%
#   add_row(Training = Mechanistic_Train_Sensitivity_CI_High, Testing = Mechanistic_Test_Sensitivity_CI_High, metric = "Sensitvity_CIHigh", .before = 12) %>%
#   add_row(Training = Mechanistic_Train_Specificity_CI_Low, Testing = Mechanistic_Test_Specificity_CI_Low, metric = "Specificity_CILow", .before = 14) %>%
#   add_row(Training = Mechanistic_Train_Specificity_CI_High, Testing = Mechanistic_Test_Specificity_CI_High, metric = "Specificity_CIHigh", .before = 15) %>%
#   add_row(Training = Mechanistic_Train_MCC_CI_Low, Testing = Mechanistic_Test_MCC_CI_Low, metric = "MCC_CILow", .before = 17) %>%
#   add_row(Training = Mechanistic_Train_MCC_CI_High, Testing = Mechanistic_Test_MCC_CI_High, metric = "MCC_CIHigh", .before = 18)
# 
# rownames(MeansMechanistic) <- MeansMechanistic$metric  
# MeansMechanistic$model_type <- rep("Mechanistic", 18)
# MeansMechanistic$Algorithm <- rep("SVM", 18)
# 
# 
# SVM_CV_AllMetrics <- rbind(MeansAgnostic, MeansMechanistic)
# 
# SVM_CV_AllMetrics <- SVM_CV_AllMetrics %>%  
#   pivot_longer(cols = c(1,2), names_to = c("data_type"))
# 
# SVM_CV_AllMetrics$type <- SVM_CV_AllMetrics$metric
# SVM_CV_AllMetrics$metric <- gsub("Mean", "", SVM_CV_AllMetrics$metric)
# SVM_CV_AllMetrics$metric <- gsub("_Agnostic", "", SVM_CV_AllMetrics$metric)
# SVM_CV_AllMetrics$metric <- gsub("_Mechanistic", "", SVM_CV_AllMetrics$metric)
# SVM_CV_AllMetrics$metric <- gsub("_CILow", "", SVM_CV_AllMetrics$metric)
# SVM_CV_AllMetrics$metric <- gsub("_CIHigh", "", SVM_CV_AllMetrics$metric) 
# 
# SVM_CV_AllMetrics$type[grep("Mean", SVM_CV_AllMetrics$type)] <- "Mean" 
# SVM_CV_AllMetrics$type[grep("CILow", SVM_CV_AllMetrics$type)] <- "CILow" 
# SVM_CV_AllMetrics$type[grep("CIHigh", SVM_CV_AllMetrics$type)] <- "CIHigh" 
# 
# SVM_CV_AllMetrics$metric[SVM_CV_AllMetrics$metric == "Sensitvity"] <- "Sensitivity"
# 
# SVM_CV_AllMetrics <- SVM_CV_AllMetrics %>% pivot_wider(names_from = "type", values_from = c("value"))
# 
# SVM_CV_AllMetrics$data_type <- factor(SVM_CV_AllMetrics$data_type, levels = c("Training", "Testing"))  
# SVM_CV_AllMetrics$model_type <- factor(SVM_CV_AllMetrics$model_type, levels = c("Mechanistic", "Agnostic"))
# SVM_CV_AllMetrics$metric <- factor(SVM_CV_AllMetrics$metric, levels = c("AUC", "Accuracy", "BalAccuracy", "Sensitivity", "Specificity", "MCC"))
# 
# ## Plot the figure
# png(filename = "./Figs/SVM/SVMCV_ALlMetrics.png", width = 3000, height = 2000, res = 300)
# ggplot(data = SVM_CV_AllMetrics, aes(x = Mean, y = model_type, shape =  model_type, fill = model_type)) +
#   geom_linerange(aes(xmin = CILow, xmax = CIHigh)) +
#   geom_point(size = 2.5) +
#   geom_vline(xintercept = 0.5, linetype = 2) + 
#   facet_grid(metric ~ data_type, scales = "free") +
#   scale_x_continuous(lim = c(0, 1)) +
#   scale_shape_manual(values = c(21, 22)) +
#   scale_fill_manual(values = c("white", "grey50")) +
#   labs(y = NULL,
#        x = "Mean (95% Confidence Interval)") +
#   theme(legend.position = "none")
# 
# dev.off()
# 
# ## Save the Dataframe 
# save(SVM_CV_AllMetrics, file = "./Objs/SVM/SVM_CV_AllMetrics.rda")
# 
