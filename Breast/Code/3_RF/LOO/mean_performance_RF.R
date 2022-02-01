### Forest plot comparing the agnostic and mechanistic RF in the CV part using all metrics

rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)

## Agnostic RF Average performance

# Load the agnostic pairs models
AgnosticGSE25055_Out_pairs <- load("./Objs/RF/GSE25055_Out_RF_AgnosticPerformance.rda")
AgnosticGSE25065_Out_pairs <- load("./Objs/RF/GSE25065_Out_RF_AgnosticPerformance.rda")
AgnosticGSE140494_Out_pairs <- load("./Objs/RF/GSE140494_Out_RF_AgnosticPerformance.rda")
AgnosticGSE103668_Out_pairs <- load("./Objs/RF/GSE103668_Out_RF_AgnosticPerformance.rda")
AgnosticGSE20194_Out_pairs <- load("./Objs/RF/GSE20194_Out_RF_AgnosticPerformance.rda")
AgnosticGSE20271_Out_pairs <- load("./Objs/RF/GSE20271_Out_RF_AgnosticPerformance.rda")
AgnosticGSE32646_Out_pairs <- load("./Objs/RF/GSE32646_Out_RF_AgnosticPerformance.rda")

# Load the agnostic genes models
AgnosticGSE25055_Out_indvGenes <- load("./Objs/RF/GSE25055_Out_RF_indvGenes_AgnosticPerformance.rda")
AgnosticGSE25065_Out_indvGenes <- load("./Objs/RF/GSE25065_Out_RF_indvGenes_AgnosticPerformance.rda")
AgnosticGSE140494_Out_indvGenes <- load("./Objs/RF/GSE140494_Out_RF_indvGenes_AgnosticPerformance.rda")
AgnosticGSE103668_Out_indvGenes <- load("./Objs/RF/GSE103668_Out_RF_indvGenes_AgnosticPerformance.rda")
AgnosticGSE20194_Out_indvGenes <- load("./Objs/RF/GSE20194_Out_RF_indvGenes_AgnosticPerformance.rda")
AgnosticGSE20271_Out_indvGenes <- load("./Objs/RF/GSE20271_Out_RF_indvGenes_AgnosticPerformance.rda")
AgnosticGSE32646_Out_indvGenes <- load("./Objs/RF/GSE32646_Out_RF_indvGenes_AgnosticPerformance.rda")

# Load the Mech models
MechGSE25055_Out <- load("./Objs/RF/GSE25055_Out_RF_MechPerformance.rda")
MechGSE25065_Out <- load("./Objs/RF/GSE25065_Out_RF_MechPerformance.rda")
MechGSEGSE140494_Out <- load("./Objs/RF/GSE140494_Out_RF_MechPerformance.rda")
MechGSE103668_Out <- load("./Objs/RF/GSE103668_Out_RF_MechPerformance.rda")
MechGSE20194_Out <- load("./Objs/RF/GSE20194_Out_RF_MechPerformance.rda")
MechGSE20271_Out <- load("./Objs/RF/GSE20271_Out_RF_MechPerformance.rda")
MechGSE32646_Out <- load("./Objs/RF/GSE32646_Out_RF_MechPerformance.rda")




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

# GSE25055_out
GSE25055_Out_RF_AgnosticPerformance$out_study <- rep("GSE25055", nrow(GSE25055_Out_RF_AgnosticPerformance))
GSE25055_Out_RF_IndvGenes_AgnosticPerformance$out_study <- rep("GSE25055", nrow(GSE25055_Out_RF_IndvGenes_AgnosticPerformance))
GSE25055_Out_RF_MechPerformance$out_study <- rep("GSE25055", nrow(GSE25055_Out_RF_MechPerformance))

GSE25055_Out_RF_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE25055_Out_RF_AgnosticPerformance))
GSE25055_Out_RF_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE25055_Out_RF_IndvGenes_AgnosticPerformance))
GSE25055_Out_RF_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE25055_Out_RF_MechPerformance))

##############
# GSE25065_out
GSE25065_Out_RF_AgnosticPerformance$out_study <- rep("GSE25065", nrow(GSE25065_Out_RF_AgnosticPerformance))
GSE25065_Out_RF_IndvGenes_AgnosticPerformance$out_study <- rep("GSE25065", nrow(GSE25065_Out_RF_IndvGenes_AgnosticPerformance))
GSE25065_Out_RF_MechPerformance$out_study <- rep("GSE25065", nrow(GSE25065_Out_RF_MechPerformance))

GSE25065_Out_RF_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE25065_Out_RF_AgnosticPerformance))
GSE25065_Out_RF_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE25065_Out_RF_IndvGenes_AgnosticPerformance))
GSE25065_Out_RF_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE25065_Out_RF_MechPerformance))

##############
# GSE140494_out
GSE140494_Out_RF_AgnosticPerformance$out_study <- rep("GSE140494", nrow(GSE140494_Out_RF_AgnosticPerformance))
GSE140494_Out_RF_IndvGenes_AgnosticPerformance$out_study <- rep("GSE140494", nrow(GSE140494_Out_RF_IndvGenes_AgnosticPerformance))
GSE140494_Out_RF_MechPerformance$out_study <- rep("GSE140494", nrow(GSE140494_Out_RF_MechPerformance))

GSE140494_Out_RF_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE140494_Out_RF_AgnosticPerformance))
GSE140494_Out_RF_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE140494_Out_RF_IndvGenes_AgnosticPerformance))
GSE140494_Out_RF_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE140494_Out_RF_MechPerformance))

###############
# GSE103668_out
GSE103668_Out_RF_AgnosticPerformance$out_study <- rep("GSE103668", nrow(GSE103668_Out_RF_AgnosticPerformance))
GSE103668_Out_RF_IndvGenes_AgnosticPerformance$out_study <- rep("GSE103668", nrow(GSE103668_Out_RF_IndvGenes_AgnosticPerformance))
GSE103668_Out_RF_MechPerformance$out_study <- rep("GSE103668", nrow(GSE103668_Out_RF_MechPerformance))

GSE103668_Out_RF_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE103668_Out_RF_AgnosticPerformance))
GSE103668_Out_RF_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE103668_Out_RF_IndvGenes_AgnosticPerformance))
GSE103668_Out_RF_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE103668_Out_RF_MechPerformance))

##############
# GSE20194_out
GSE20194_Out_RF_AgnosticPerformance$out_study <- rep("GSE20194", nrow(GSE20194_Out_RF_AgnosticPerformance))
GSE20194_Out_RF_IndvGenes_AgnosticPerformance$out_study <- rep("GSE20194", nrow(GSE20194_Out_RF_IndvGenes_AgnosticPerformance))
GSE20194_Out_RF_MechPerformance$out_study <- rep("GSE20194", nrow(GSE20194_Out_RF_MechPerformance))

GSE20194_Out_RF_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE20194_Out_RF_AgnosticPerformance))
GSE20194_Out_RF_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE20194_Out_RF_IndvGenes_AgnosticPerformance))
GSE20194_Out_RF_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE20194_Out_RF_MechPerformance))

##############
# GSE20271_out
GSE20271_Out_RF_AgnosticPerformance$out_study <- rep("GSE20271", nrow(GSE20271_Out_RF_AgnosticPerformance))
GSE20271_Out_RF_IndvGenes_AgnosticPerformance$out_study <- rep("GSE20271", nrow(GSE20271_Out_RF_IndvGenes_AgnosticPerformance))
GSE20271_Out_RF_MechPerformance$out_study <- rep("GSE20271", nrow(GSE20271_Out_RF_MechPerformance))

GSE20271_Out_RF_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE20271_Out_RF_AgnosticPerformance))
GSE20271_Out_RF_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE20271_Out_RF_IndvGenes_AgnosticPerformance))
GSE20271_Out_RF_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE20271_Out_RF_MechPerformance))

#############
# GSE32646_out
GSE32646_Out_RF_AgnosticPerformance$out_study <- rep("GSE32646", nrow(GSE32646_Out_RF_AgnosticPerformance))
GSE32646_Out_RF_IndvGenes_AgnosticPerformance$out_study <- rep("GSE32646", nrow(GSE32646_Out_RF_IndvGenes_AgnosticPerformance))
GSE32646_Out_RF_MechPerformance$out_study <- rep("GSE32646", nrow(GSE32646_Out_RF_MechPerformance))

GSE32646_Out_RF_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE32646_Out_RF_AgnosticPerformance))
GSE32646_Out_RF_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE32646_Out_RF_IndvGenes_AgnosticPerformance))
GSE32646_Out_RF_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE32646_Out_RF_MechPerformance))


###############
All_RF_CV <- rbind(GSE25055_Out_RF_AgnosticPerformance, GSE25055_Out_RF_IndvGenes_AgnosticPerformance, GSE25055_Out_RF_MechPerformance,
                     GSE25065_Out_RF_AgnosticPerformance, GSE25065_Out_RF_IndvGenes_AgnosticPerformance, GSE25065_Out_RF_MechPerformance,
                     GSE140494_Out_RF_AgnosticPerformance, GSE140494_Out_RF_IndvGenes_AgnosticPerformance, GSE140494_Out_RF_MechPerformance, 
                     GSE103668_Out_RF_AgnosticPerformance, GSE103668_Out_RF_IndvGenes_AgnosticPerformance, GSE103668_Out_RF_MechPerformance, 
                     GSE20194_Out_RF_AgnosticPerformance, GSE20194_Out_RF_IndvGenes_AgnosticPerformance, GSE20194_Out_RF_MechPerformance,
                     GSE20271_Out_RF_AgnosticPerformance, GSE20271_Out_RF_IndvGenes_AgnosticPerformance, GSE20271_Out_RF_MechPerformance,
                     GSE32646_Out_RF_AgnosticPerformance, GSE32646_Out_RF_IndvGenes_AgnosticPerformance, GSE32646_Out_RF_MechPerformance)

All_RF_CV$metric <- rownames(All_RF_CV)
All_RF_CV$metric <- gsub('[[:digit:]]+', '', All_RF_CV$metric )


##########################################################
##########################################################
## Calculate the average statistics (Across the 7 repeats)

#####
# average stats in the agnostic pairs models

## AUC
# Mean AUC in both the training and testing data
MeanAUC_Agnostic_pairs <-  apply(All_RF_CV %>%
                             filter(model_type == "Agnostic_pairs") %>%
                             filter(metric == "AUC") %>%
                             select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_pairs <-  apply(All_RF_CV %>%
                            filter(model_type == "Agnostic_pairs") %>%
                            filter(metric == "AUC") %>%
                            select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_pairs <-  apply(All_RF_CV %>%
                                  filter(model_type == "Agnostic_pairs") %>%
                                  filter(metric == "Accuracy") %>%
                                  select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_pairs <-  apply(All_RF_CV %>%
                                 filter(model_type == "Agnostic_pairs") %>%
                                 filter(metric == "Accuracy") %>%
                                 select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_pairs <-  apply(All_RF_CV %>%
                                     filter(model_type == "Agnostic_pairs") %>%
                                     filter(metric == "Bal.Accuracy") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_pairs <-  apply(All_RF_CV %>%
                                    filter(model_type == "Agnostic_pairs") %>%
                                    filter(metric == "Bal.Accuracy") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_pairs <-  apply(All_RF_CV %>%
                                     filter(model_type == "Agnostic_pairs") %>%
                                     filter(metric == "Sensitivity") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_pairs <-  apply(All_RF_CV %>%
                                    filter(model_type == "Agnostic_pairs") %>%
                                    filter(metric == "Sensitivity") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_pairs <-  apply(All_RF_CV %>%
                                     filter(model_type == "Agnostic_pairs") %>%
                                     filter(metric == "Specificity") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_pairs <-  apply(All_RF_CV %>%
                                    filter(model_type == "Agnostic_pairs") %>%
                                    filter(metric == "Specificity") %>%
                                    select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_pairs <-  apply(All_RF_CV %>%
                             filter(model_type == "Agnostic_pairs") %>%
                             filter(metric == "MCC") %>%
                             select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_pairs <-  apply(All_RF_CV %>%
                            filter(model_type == "Agnostic_pairs") %>%
                            filter(metric == "MCC") %>%
                            select(c("Training", "Testing")), 2, sd)



## Group them together
MeansAgnostic_pairs <- rbind(MeanAUC_Agnostic_pairs, SD_AUC_Agnostic_pairs, 
                       MeanAccuracy_Agnostic_pairs, SD_Accuracy_Agnostic_pairs, 
                       MeanBalAccuracy_Agnostic_pairs, SD_BalAccuracy_Agnostic_pairs, 
                       MeanSensitivity_Agnostic_pairs, SD_Sensitivity_Agnostic_pairs, 
                       MeanSpecificity_Agnostic_pairs, SD_Specificity_Agnostic_pairs,
                       MeanMCC_Agnostic_pairs, SD_MCC_Agnostic_pairs)
# Round the Numbers
MeansAgnostic_pairs <- apply(MeansAgnostic_pairs, 2, round, digits = 2) 

##############################################
# average stats in the agnostic genes models

## AUC
# Mean AUC in both the training and testing data
MeanAUC_Agnostic_genes <-  apply(All_RF_CV %>%
                                   filter(model_type == "Agnostic_genes") %>%
                                   filter(metric == "AUC") %>%
                                   select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_genes <-  apply(All_RF_CV %>%
                                  filter(model_type == "Agnostic_genes") %>%
                                  filter(metric == "AUC") %>%
                                  select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_genes <-  apply(All_RF_CV %>%
                                        filter(model_type == "Agnostic_genes") %>%
                                        filter(metric == "Accuracy") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_genes <-  apply(All_RF_CV %>%
                                       filter(model_type == "Agnostic_genes") %>%
                                       filter(metric == "Accuracy") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_genes <-  apply(All_RF_CV %>%
                                           filter(model_type == "Agnostic_genes") %>%
                                           filter(metric == "Bal.Accuracy") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_genes <-  apply(All_RF_CV %>%
                                          filter(model_type == "Agnostic_genes") %>%
                                          filter(metric == "Bal.Accuracy") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_genes <-  apply(All_RF_CV %>%
                                           filter(model_type == "Agnostic_genes") %>%
                                           filter(metric == "Sensitivity") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_genes <-  apply(All_RF_CV %>%
                                          filter(model_type == "Agnostic_genes") %>%
                                          filter(metric == "Sensitivity") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_genes <-  apply(All_RF_CV %>%
                                           filter(model_type == "Agnostic_genes") %>%
                                           filter(metric == "Specificity") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_genes <-  apply(All_RF_CV %>%
                                          filter(model_type == "Agnostic_genes") %>%
                                          filter(metric == "Specificity") %>%
                                          select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_genes <-  apply(All_RF_CV %>%
                                   filter(model_type == "Agnostic_genes") %>%
                                   filter(metric == "MCC") %>%
                                   select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_genes <-  apply(All_RF_CV %>%
                                  filter(model_type == "Agnostic_genes") %>%
                                  filter(metric == "MCC") %>%
                                  select(c("Training", "Testing")), 2, sd)



## Group them together
MeansAgnostic_genes <- rbind(MeanAUC_Agnostic_genes, SD_AUC_Agnostic_genes, 
                             MeanAccuracy_Agnostic_genes, SD_Accuracy_Agnostic_genes, 
                             MeanBalAccuracy_Agnostic_genes, SD_BalAccuracy_Agnostic_genes, 
                             MeanSensitivity_Agnostic_genes, SD_Sensitivity_Agnostic_genes, 
                             MeanSpecificity_Agnostic_genes, SD_Specificity_Agnostic_genes,
                             MeanMCC_Agnostic_genes, SD_MCC_Agnostic_genes)
# Round the Numbers
MeansAgnostic_genes <- apply(MeansAgnostic_genes, 2, round, digits = 2) 


################################################################
############################################
##
#####
# average stats in the Mechanistic models

## AUC
# Mean AUC in both the training and testing data
MeanAUC_Mechanistic <-  apply(All_RF_CV %>%
                                filter(model_type == "Mechanistic") %>%
                                filter(metric == "AUC") %>%
                                select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Mechanistic <-  apply(All_RF_CV %>%
                               filter(model_type == "Mechanistic") %>%
                               filter(metric == "AUC") %>%
                               select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Mechanistic <-  apply(All_RF_CV %>%
                                     filter(model_type == "Mechanistic") %>%
                                     filter(metric == "Accuracy") %>%
                                     select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Mechanistic <-  apply(All_RF_CV %>%
                                    filter(model_type == "Mechanistic") %>%
                                    filter(metric == "Accuracy") %>%
                                    select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Mechanistic <-  apply(All_RF_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Bal.Accuracy") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Mechanistic <-  apply(All_RF_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Bal.Accuracy") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Mechanistic <-  apply(All_RF_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Sensitivity") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Mechanistic <-  apply(All_RF_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Sensitivity") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Mechanistic <-  apply(All_RF_CV %>%
                                        filter(model_type == "Mechanistic") %>%
                                        filter(metric == "Specificity") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Mechanistic <-  apply(All_RF_CV %>%
                                       filter(model_type == "Mechanistic") %>%
                                       filter(metric == "Specificity") %>%
                                       select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Mechanistic <-  apply(All_RF_CV %>%
                                filter(model_type == "Mechanistic") %>%
                                filter(metric == "MCC") %>%
                                select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Mechanistic <-  apply(All_RF_CV %>%
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

################
# Functions to calculate the 95% CI from the mean and SD

#CalculateCIHigh <- function(x, y){Z <- x + 1.96 * (y/sqrt(5))
#                         return(Z)
#                         }


#CalculateCILow <- function(x, y){Z <- x - 1.96 * (y/sqrt(5))
#return(Z)
#}


MeansAgnostic_pairs <- data.frame(MeansAgnostic_pairs)
MeansAgnostic_pairs$metric <- rownames(MeansAgnostic_pairs)

MeansAgnostic_genes <- data.frame(MeansAgnostic_genes)
MeansAgnostic_genes$metric <- rownames(MeansAgnostic_genes)

MeansMechanistic <- data.frame(MeansMechanistic)
MeansMechanistic$metric <- rownames(MeansMechanistic)

############################3
## Use the mean and SD to calculate the 95% CI 

### Agnostic models

## Training data

## CIs for AUC
# Agnostic_Train_AUC_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanAUC_Agnostic","Training"], y = MeansAgnostic["SD_AUC_Agnostic", "Training"])
# Agnostic_Train_AUC_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanAUC_Agnostic","Training"], y = MeansAgnostic["SD_AUC_Agnostic", "Training"])
# 
# # CIs for Accuracy
# Agnostic_Train_Accuracy_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanAccuracy_Agnostic","Training"], y = MeansAgnostic["SD_Accuracy_Agnostic", "Training"])
# Agnostic_Train_Accuracy_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanAccuracy_Agnostic","Training"], y = MeansAgnostic["SD_Accuracy_Agnostic", "Training"])
# 
# # CIs for Bal. Accuracy
# Agnostic_Train_BalAccuracy_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanBalAccuracy_Agnostic","Training"], y = MeansAgnostic["SD_BalAccuracy_Agnostic", "Training"])
# Agnostic_Train_BalAccuracy_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanBalAccuracy_Agnostic","Training"], y = MeansAgnostic["SD_BalAccuracy_Agnostic", "Training"])
# 
# # CIs for sensitivity
# Agnostic_Train_Sensitivity_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanSensitivity_Agnostic","Training"], y = MeansAgnostic["SD_Sensitivity_Agnostic", "Training"])
# Agnostic_Train_Sensitivity_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanSensitivity_Agnostic","Training"], y = MeansAgnostic["SD_Sensitivity_Agnostic", "Training"])
# 
# # CIs for Specificity
# Agnostic_Train_Specificity_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanSpecificity_Agnostic","Training"], y = MeansAgnostic["SD_Specificity_Agnostic", "Training"])
# Agnostic_Train_Specificity_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanSpecificity_Agnostic","Training"], y = MeansAgnostic["SD_Specificity_Agnostic", "Training"])
# 
# # CIs for MCC
# Agnostic_Train_MCC_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanMCC_Agnostic","Training"], y = MeansAgnostic["SD_MCC_Agnostic", "Training"])
# Agnostic_Train_MCC_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanMCC_Agnostic","Training"], y = MeansAgnostic["SD_MCC_Agnostic", "Training"])
# 
# 
# ## Testing data
# 
# ## CIs for AUC
# Agnostic_Test_AUC_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanAUC_Agnostic","Testing"], y = MeansAgnostic["SD_AUC_Agnostic", "Testing"])
# Agnostic_Test_AUC_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanAUC_Agnostic","Testing"], y = MeansAgnostic["SD_AUC_Agnostic", "Testing"])
# 
# # CIs for Accuracy
# Agnostic_Test_Accuracy_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanAccuracy_Agnostic","Testing"], y = MeansAgnostic["SD_Accuracy_Agnostic", "Testing"])
# Agnostic_Test_Accuracy_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanAccuracy_Agnostic","Testing"], y = MeansAgnostic["SD_Accuracy_Agnostic", "Testing"])
# 
# # CIs for Bal. Accuracy
# Agnostic_Test_BalAccuracy_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanBalAccuracy_Agnostic","Testing"], y = MeansAgnostic["SD_BalAccuracy_Agnostic", "Testing"])
# Agnostic_Test_BalAccuracy_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanBalAccuracy_Agnostic","Testing"], y = MeansAgnostic["SD_BalAccuracy_Agnostic", "Testing"])
# 
# # CIs for sensitivity
# Agnostic_Test_Sensitivity_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanSensitivity_Agnostic","Testing"], y = MeansAgnostic["SD_Sensitivity_Agnostic", "Testing"])
# Agnostic_Test_Sensitivity_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanSensitivity_Agnostic","Testing"], y = MeansAgnostic["SD_Sensitivity_Agnostic", "Testing"])
# 
# # CIs for Specificity
# Agnostic_Test_Specificity_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanSpecificity_Agnostic","Testing"], y = MeansAgnostic["SD_Specificity_Agnostic", "Testing"])
# Agnostic_Test_Specificity_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanSpecificity_Agnostic","Testing"], y = MeansAgnostic["SD_Specificity_Agnostic", "Testing"])
# 
# # CIs for MCC
# Agnostic_Test_MCC_CI_Low <- CalculateCILow(x = MeansAgnostic["MeanMCC_Agnostic","Testing"], y = MeansAgnostic["SD_MCC_Agnostic", "Testing"])
# Agnostic_Test_MCC_CI_High <- CalculateCIHigh(x = MeansAgnostic["MeanMCC_Agnostic","Testing"], y = MeansAgnostic["SD_MCC_Agnostic", "Testing"])
# 
# ##############################
# ### Mechanistic models
# 
# ## Training data
# 
# ## CIs for AUC
# Mechanistic_Train_AUC_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanAUC_Mechanistic","Training"], y = MeansMechanistic["SD_AUC_Mechanistic", "Training"])
# Mechanistic_Train_AUC_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanAUC_Mechanistic","Training"], y = MeansMechanistic["SD_AUC_Mechanistic", "Training"])
# 
# # CIs for Accuracy
# Mechanistic_Train_Accuracy_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanAccuracy_Mechanistic","Training"], y = MeansMechanistic["SD_Accuracy_Mechanistic", "Training"])
# Mechanistic_Train_Accuracy_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanAccuracy_Mechanistic","Training"], y = MeansMechanistic["SD_Accuracy_Mechanistic", "Training"])
# 
# # CIs for Bal. Accuracy
# Mechanistic_Train_BalAccuracy_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanBalAccuracy_Mechanistic","Training"], y = MeansMechanistic["SD_BalAccuracy_Mechanistic", "Training"])
# Mechanistic_Train_BalAccuracy_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanBalAccuracy_Mechanistic","Training"], y = MeansMechanistic["SD_BalAccuracy_Mechanistic", "Training"])
# 
# # CIs for sensitivity
# Mechanistic_Train_Sensitivity_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanSensitivity_Mechanistic","Training"], y = MeansMechanistic["SD_Sensitivity_Mechanistic", "Training"])
# Mechanistic_Train_Sensitivity_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanSensitivity_Mechanistic","Training"], y = MeansMechanistic["SD_Sensitivity_Mechanistic", "Training"])
# 
# # CIs for Specificity
# Mechanistic_Train_Specificity_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanSpecificity_Mechanistic","Training"], y = MeansMechanistic["SD_Specificity_Mechanistic", "Training"])
# Mechanistic_Train_Specificity_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanSpecificity_Mechanistic","Training"], y = MeansMechanistic["SD_Specificity_Mechanistic", "Training"])
# 
# # CIs for MCC
# Mechanistic_Train_MCC_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanMCC_Mechanistic","Training"], y = MeansMechanistic["SD_MCC_Mechanistic", "Training"])
# Mechanistic_Train_MCC_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanMCC_Mechanistic","Training"], y = MeansMechanistic["SD_MCC_Mechanistic", "Training"])
# 
# 
# ## Testing data
# 
# ## CIs for AUC
# Mechanistic_Test_AUC_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanAUC_Mechanistic","Testing"], y = MeansMechanistic["SD_AUC_Mechanistic", "Testing"])
# Mechanistic_Test_AUC_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanAUC_Mechanistic","Testing"], y = MeansMechanistic["SD_AUC_Mechanistic", "Testing"])
# 
# # CIs for Accuracy
# Mechanistic_Test_Accuracy_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanAccuracy_Mechanistic","Testing"], y = MeansMechanistic["SD_Accuracy_Mechanistic", "Testing"])
# Mechanistic_Test_Accuracy_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanAccuracy_Mechanistic","Testing"], y = MeansMechanistic["SD_Accuracy_Mechanistic", "Testing"])
# 
# # CIs for Bal. Accuracy
# Mechanistic_Test_BalAccuracy_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanBalAccuracy_Mechanistic","Testing"], y = MeansMechanistic["SD_BalAccuracy_Mechanistic", "Testing"])
# Mechanistic_Test_BalAccuracy_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanBalAccuracy_Mechanistic","Testing"], y = MeansMechanistic["SD_BalAccuracy_Mechanistic", "Testing"])
# 
# # CIs for sensitivity
# Mechanistic_Test_Sensitivity_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanSensitivity_Mechanistic","Testing"], y = MeansMechanistic["SD_Sensitivity_Mechanistic", "Testing"])
# Mechanistic_Test_Sensitivity_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanSensitivity_Mechanistic","Testing"], y = MeansMechanistic["SD_Sensitivity_Mechanistic", "Testing"])
# 
# # CIs for Specificity
# Mechanistic_Test_Specificity_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanSpecificity_Mechanistic","Testing"], y = MeansMechanistic["SD_Specificity_Mechanistic", "Testing"])
# Mechanistic_Test_Specificity_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanSpecificity_Mechanistic","Testing"], y = MeansMechanistic["SD_Specificity_Mechanistic", "Testing"])
# 
# # CIs for MCC
# Mechanistic_Test_MCC_CI_Low <- CalculateCILow(x = MeansMechanistic["MeanMCC_Mechanistic","Testing"], y = MeansMechanistic["SD_MCC_Mechanistic", "Testing"])
# Mechanistic_Test_MCC_CI_High <- CalculateCIHigh(x = MeansMechanistic["MeanMCC_Mechanistic","Testing"], y = MeansMechanistic["SD_MCC_Mechanistic", "Testing"])
# 


######################################
MeansAgnostic <- MeansAgnostic[-grep("SD", rownames(MeansAgnostic)), ]
MeansMechanistic <- MeansMechanistic[-grep("SD", rownames(MeansMechanistic)), ]

MeansAgnostic <- MeansAgnostic %>%
  add_row(Training = Agnostic_Train_AUC_CI_Low, Testing = Agnostic_Test_AUC_CI_Low, metric = "AUC_CILow", .before = 2) %>%
  add_row(Training = Agnostic_Train_AUC_CI_High, Testing = Agnostic_Test_AUC_CI_High, metric = "AUC_CIHigh", .before = 3) %>%
  add_row(Training = Agnostic_Train_Accuracy_CI_Low, Testing = Agnostic_Test_Accuracy_CI_Low, metric = "Accuracy_CILow", .before = 5) %>%
  add_row(Training = Agnostic_Train_Accuracy_CI_High, Testing = Agnostic_Test_Accuracy_CI_High, metric = "Accuracy_CIHigh", .before = 6) %>%
  add_row(Training = Agnostic_Train_BalAccuracy_CI_Low, Testing = Agnostic_Test_BalAccuracy_CI_Low, metric = "BalAccuracy_CILow", .before = 8) %>%
  add_row(Training = Agnostic_Train_BalAccuracy_CI_High, Testing = Agnostic_Test_BalAccuracy_CI_High, metric = "BalAccuracy_CIHigh", .before = 9) %>%
  add_row(Training = Agnostic_Train_Sensitivity_CI_Low, Testing = Agnostic_Test_Sensitivity_CI_Low, metric = "Sensitivity_CILow", .before = 11) %>%
  add_row(Training = Agnostic_Train_Sensitivity_CI_High, Testing = Agnostic_Test_Sensitivity_CI_High, metric = "Sensitvity_CIHigh", .before = 12) %>%
  add_row(Training = Agnostic_Train_Specificity_CI_Low, Testing = Agnostic_Test_Specificity_CI_Low, metric = "Specificity_CILow", .before = 14) %>%
  add_row(Training = Agnostic_Train_Specificity_CI_High, Testing = Agnostic_Test_Specificity_CI_High, metric = "Specificity_CIHigh", .before = 15) %>%
  add_row(Training = Agnostic_Train_MCC_CI_Low, Testing = Agnostic_Test_MCC_CI_Low, metric = "MCC_CILow", .before = 17) %>%
  add_row(Training = Agnostic_Train_MCC_CI_High, Testing = Agnostic_Test_MCC_CI_High, metric = "MCC_CIHigh", .before = 18)

rownames(MeansAgnostic) <- MeansAgnostic$metric  
MeansAgnostic$model_type <- rep("Agnostic", 18)
MeansAgnostic$Algorithm <- rep("RF", 18)

MeansMechanistic <- MeansMechanistic %>%
  add_row(Training = Mechanistic_Train_AUC_CI_Low, Testing = Mechanistic_Test_AUC_CI_Low, metric = "AUC_CILow", .before = 2) %>%
  add_row(Training = Mechanistic_Train_AUC_CI_High, Testing = Mechanistic_Test_AUC_CI_High, metric = "AUC_CIHigh", .before = 3) %>%
  add_row(Training = Mechanistic_Train_Accuracy_CI_Low, Testing = Mechanistic_Test_Accuracy_CI_Low, metric = "Accuracy_CILow", .before = 5) %>%
  add_row(Training = Mechanistic_Train_Accuracy_CI_High, Testing = Mechanistic_Test_Accuracy_CI_High, metric = "Accuracy_CIHigh", .before = 6) %>%
  add_row(Training = Mechanistic_Train_BalAccuracy_CI_Low, Testing = Mechanistic_Test_BalAccuracy_CI_Low, metric = "BalAccuracy_CILow", .before = 8) %>%
  add_row(Training = Mechanistic_Train_BalAccuracy_CI_High, Testing = Mechanistic_Test_BalAccuracy_CI_High, metric = "BalAccuracy_CIHigh", .before = 9) %>%
  add_row(Training = Mechanistic_Train_Sensitivity_CI_Low, Testing = Mechanistic_Test_Sensitivity_CI_Low, metric = "Sensitivity_CILow", .before = 11) %>%
  add_row(Training = Mechanistic_Train_Sensitivity_CI_High, Testing = Mechanistic_Test_Sensitivity_CI_High, metric = "Sensitvity_CIHigh", .before = 12) %>%
  add_row(Training = Mechanistic_Train_Specificity_CI_Low, Testing = Mechanistic_Test_Specificity_CI_Low, metric = "Specificity_CILow", .before = 14) %>%
  add_row(Training = Mechanistic_Train_Specificity_CI_High, Testing = Mechanistic_Test_Specificity_CI_High, metric = "Specificity_CIHigh", .before = 15) %>%
  add_row(Training = Mechanistic_Train_MCC_CI_Low, Testing = Mechanistic_Test_MCC_CI_Low, metric = "MCC_CILow", .before = 17) %>%
  add_row(Training = Mechanistic_Train_MCC_CI_High, Testing = Mechanistic_Test_MCC_CI_High, metric = "MCC_CIHigh", .before = 18)

rownames(MeansMechanistic) <- MeansMechanistic$metric  
MeansMechanistic$model_type <- rep("Mechanistic", 18)
MeansMechanistic$Algorithm <- rep("RF", 18)


RF_CV_AllMetrics <- rbind(MeansAgnostic, MeansMechanistic)

RF_CV_AllMetrics <- RF_CV_AllMetrics %>%  
  pivot_longer(cols = c(1,2), names_to = c("data_type"))

RF_CV_AllMetrics$type <- RF_CV_AllMetrics$metric
RF_CV_AllMetrics$metric <- gsub("Mean", "", RF_CV_AllMetrics$metric)
RF_CV_AllMetrics$metric <- gsub("_Agnostic", "", RF_CV_AllMetrics$metric)
RF_CV_AllMetrics$metric <- gsub("_Mechanistic", "", RF_CV_AllMetrics$metric)
RF_CV_AllMetrics$metric <- gsub("_CILow", "", RF_CV_AllMetrics$metric)
RF_CV_AllMetrics$metric <- gsub("_CIHigh", "", RF_CV_AllMetrics$metric) 

RF_CV_AllMetrics$type[grep("Mean", RF_CV_AllMetrics$type)] <- "Mean" 
RF_CV_AllMetrics$type[grep("CILow", RF_CV_AllMetrics$type)] <- "CILow" 
RF_CV_AllMetrics$type[grep("CIHigh", RF_CV_AllMetrics$type)] <- "CIHigh" 

RF_CV_AllMetrics$metric[RF_CV_AllMetrics$metric == "Sensitvity"] <- "Sensitivity"

RF_CV_AllMetrics <- RF_CV_AllMetrics %>% pivot_wider(names_from = "type", values_from = c("value"))

RF_CV_AllMetrics$data_type <- factor(RF_CV_AllMetrics$data_type, levels = c("Training", "Testing"))  
RF_CV_AllMetrics$model_type <- factor(RF_CV_AllMetrics$model_type, levels = c("Mechanistic", "Agnostic"))
RF_CV_AllMetrics$metric <- factor(RF_CV_AllMetrics$metric, levels = c("AUC", "Accuracy", "BalAccuracy", "Sensitivity", "Specificity", "MCC"))

## Plot the figure
png(filename = "./Figs/RF/RFCV_ALlMetrics.png", width = 3000, height = 2000, res = 300)
ggplot(data = RF_CV_AllMetrics, aes(x = Mean, y = model_type, shape =  model_type, fill = model_type)) +
  geom_linerange(aes(xmin = CILow, xmax = CIHigh)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = 2) + 
  facet_grid(metric ~ data_type, scales = "free") +
  scale_x_continuous(lim = c(0, 1)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = c("white", "grey50")) +
  labs(y = NULL,
       x = "Mean (95% Confidence Interval)") +
  theme(legend.position = "none")

dev.off()

## Save the Dataframe 
save(RF_CV_AllMetrics, file = "./Objs/RF_CV_AllMetrics.rda")

