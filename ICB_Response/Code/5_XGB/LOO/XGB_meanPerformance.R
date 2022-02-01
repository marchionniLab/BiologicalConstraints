### Forest plot comparing the agnostic and mechanistic XGB in the CV part using all metrics

rm(list = ls())

library(dplyr)
library(tidyverse)
library(ggplot2)

## Agnostic XGB Average performance

# Load the agnostic pairs models
AgnosticGSE78220_Out <- load("./Objs/XGB/GSE78220_Out_XGB_AgnosticPerformance.rda")
AgnosticGSE91061_Out <- load("./Objs/XGB/GSE91061_Out_XGB_AgnosticPerformance.rda")
AgnosticGSE115821_Out <- load("./Objs/XGB/GSE115821_Out_XGB_AgnosticPerformance.rda")
AgnosticTCGA_Out <- load("./Objs/XGB/TCGA_Out_XGB_AgnosticPerformance.rda")
AgnosticVanAllen_Out <- load("./Objs/XGB/VanAllen_Out_XGB_AgnosticPerformance.rda")

# Load the agnostic genes models
AgnosticGSE78220_Out_indvGenes <- load("./Objs/XGB/GSE78220_Out_XGB_indvGenes_AgnosticPerformance.rda")
AgnosticGSE91061_Out_indvGenes <- load("./Objs/XGB/GSE91061_Out_XGB_indvGenes_AgnosticPerformance.rda")
AgnosticGSE115821_Out_indvGenes <- load("./Objs/XGB/GSE115821_Out_XGB_indvGenes_AgnosticPerformance.rda")
AgnosticTCGA_Out_indvGenes <- load("./Objs/XGB/TCGA_Out_XGB_indvGenes_AgnosticPerformance.rda")
AgnosticVanAllen_Out_indvGenes <- load("./Objs/XGB/VanAllen_Out_XGB_indvGenes_AgnosticPerformance.rda")

# Load the Mech models
MechGSE78220_Out <- load("./Objs/XGB/GSE78220_Out_XGB_MechPerformance.rda")
MechGSE91061_Out <- load("./Objs/XGB/GSE91061_Out_XGB_MechPerformance.rda")
MechGSEGSE115821_Out <- load("./Objs/XGB/GSE115821_Out_XGB_MechPerformance.rda")
MechTCGA_Out <- load("./Objs/XGB/TCGA_Out_XGB_MechPerformance.rda")
MechVanAllen_Out <- load("./Objs/XGB/VanAllen_Out_XGB_MechPerformance.rda")

###################

# GSE78220_out
GSE78220_Out_XGB_AgnosticPerformance$out_study <- rep("GSE78220", nrow(GSE78220_Out_XGB_AgnosticPerformance))
GSE78220_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("GSE78220", nrow(GSE78220_Out_XGB_IndvGenes_AgnosticPerformance))
GSE78220_Out_XGB_MechPerformance$out_study <- rep("GSE78220", nrow(GSE78220_Out_XGB_MechPerformance))

GSE78220_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE78220_Out_XGB_AgnosticPerformance))
GSE78220_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE78220_Out_XGB_IndvGenes_AgnosticPerformance))
GSE78220_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE78220_Out_XGB_MechPerformance))

# GSE91061_out
GSE91061_Out_XGB_AgnosticPerformance$out_study <- rep("GSE91061", nrow(GSE91061_Out_XGB_AgnosticPerformance))
GSE91061_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("GSE91061", nrow(GSE91061_Out_XGB_IndvGenes_AgnosticPerformance))
GSE91061_Out_XGB_MechPerformance$out_study <- rep("GSE91061", nrow(GSE91061_Out_XGB_MechPerformance))

GSE91061_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE91061_Out_XGB_AgnosticPerformance))
GSE91061_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE91061_Out_XGB_IndvGenes_AgnosticPerformance))
GSE91061_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE91061_Out_XGB_MechPerformance))

# GSE115821_out
GSE115821_Out_XGB_AgnosticPerformance$out_study <- rep("GSE115821", nrow(GSE115821_Out_XGB_AgnosticPerformance))
GSE115821_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("GSE115821", nrow(GSE115821_Out_XGB_IndvGenes_AgnosticPerformance))
GSE115821_Out_XGB_MechPerformance$out_study <- rep("GSE115821", nrow(GSE115821_Out_XGB_MechPerformance))

GSE115821_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(GSE115821_Out_XGB_AgnosticPerformance))
GSE115821_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(GSE115821_Out_XGB_IndvGenes_AgnosticPerformance))
GSE115821_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(GSE115821_Out_XGB_MechPerformance))

# TCGA_out
TCGA_Out_XGB_AgnosticPerformance$out_study <- rep("TCGA", nrow(TCGA_Out_XGB_AgnosticPerformance))
TCGA_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("TCGA", nrow(TCGA_Out_XGB_IndvGenes_AgnosticPerformance))
TCGA_Out_XGB_MechPerformance$out_study <- rep("TCGA", nrow(TCGA_Out_XGB_MechPerformance))

TCGA_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(TCGA_Out_XGB_AgnosticPerformance))
TCGA_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(TCGA_Out_XGB_IndvGenes_AgnosticPerformance))
TCGA_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(TCGA_Out_XGB_MechPerformance))

# VanAllen_out
VanAllen_Out_XGB_AgnosticPerformance$out_study <- rep("VanAllen", nrow(VanAllen_Out_XGB_AgnosticPerformance))
VanAllen_Out_XGB_IndvGenes_AgnosticPerformance$out_study <- rep("VanAllen", nrow(VanAllen_Out_XGB_IndvGenes_AgnosticPerformance))
VanAllen_Out_XGB_MechPerformance$out_study <- rep("VanAllen", nrow(VanAllen_Out_XGB_MechPerformance))

VanAllen_Out_XGB_AgnosticPerformance$model_type <- rep("Agnostic_pairs", nrow(VanAllen_Out_XGB_AgnosticPerformance))
VanAllen_Out_XGB_IndvGenes_AgnosticPerformance$model_type <- rep("Agnostic_genes", nrow(VanAllen_Out_XGB_IndvGenes_AgnosticPerformance))
VanAllen_Out_XGB_MechPerformance$model_type <- rep("Mechanistic", nrow(VanAllen_Out_XGB_MechPerformance))


All_XGB_CV <- rbind(GSE78220_Out_XGB_AgnosticPerformance, GSE78220_Out_XGB_IndvGenes_AgnosticPerformance, GSE78220_Out_XGB_MechPerformance,
                   GSE91061_Out_XGB_AgnosticPerformance, GSE91061_Out_XGB_IndvGenes_AgnosticPerformance, GSE91061_Out_XGB_MechPerformance,
                   GSE115821_Out_XGB_AgnosticPerformance, GSE115821_Out_XGB_IndvGenes_AgnosticPerformance, GSE115821_Out_XGB_MechPerformance, 
                   TCGA_Out_XGB_AgnosticPerformance, TCGA_Out_XGB_IndvGenes_AgnosticPerformance, TCGA_Out_XGB_MechPerformance, 
                   VanAllen_Out_XGB_AgnosticPerformance, VanAllen_Out_XGB_IndvGenes_AgnosticPerformance, VanAllen_Out_XGB_MechPerformance
)

All_XGB_CV$metric <- rownames(All_XGB_CV)
All_XGB_CV$metric <- gsub('[[:digit:]]+', '', All_XGB_CV$metric )


##########################################################
##########################################################
## Calculate the average statistics (Across the 5 repeats)


#####
# average stats in the agnostic pairs models

## AUC
# Mean AUC in both the training and testing data
MeanAUC_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                   filter(model_type == "Agnostic_pairs") %>%
                                   filter(metric == "AUC") %>%
                                   select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                  filter(model_type == "Agnostic_pairs") %>%
                                  filter(metric == "AUC") %>%
                                  select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                        filter(model_type == "Agnostic_pairs") %>%
                                        filter(metric == "Accuracy") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                       filter(model_type == "Agnostic_pairs") %>%
                                       filter(metric == "Accuracy") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                           filter(model_type == "Agnostic_pairs") %>%
                                           filter(metric == "Bal.Accuracy") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                          filter(model_type == "Agnostic_pairs") %>%
                                          filter(metric == "Bal.Accuracy") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                           filter(model_type == "Agnostic_pairs") %>%
                                           filter(metric == "Sensitivity") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                          filter(model_type == "Agnostic_pairs") %>%
                                          filter(metric == "Sensitivity") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                           filter(model_type == "Agnostic_pairs") %>%
                                           filter(metric == "Specificity") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                          filter(model_type == "Agnostic_pairs") %>%
                                          filter(metric == "Specificity") %>%
                                          select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_pairs <-  apply(All_XGB_CV %>%
                                   filter(model_type == "Agnostic_pairs") %>%
                                   filter(metric == "MCC") %>%
                                   select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_pairs <-  apply(All_XGB_CV %>%
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
MeanAUC_Agnostic_genes <-  apply(All_XGB_CV %>%
                                   filter(model_type == "Agnostic_genes") %>%
                                   filter(metric == "AUC") %>%
                                   select(c("Training", "Testing")), 2, mean)


# AUC SD in both the training and testing data
SD_AUC_Agnostic_genes <-  apply(All_XGB_CV %>%
                                  filter(model_type == "Agnostic_genes") %>%
                                  filter(metric == "AUC") %>%
                                  select(c("Training", "Testing")), 2, sd)

### Accuracy
# Mean Accuracy in both the training and testing data
MeanAccuracy_Agnostic_genes <-  apply(All_XGB_CV %>%
                                        filter(model_type == "Agnostic_genes") %>%
                                        filter(metric == "Accuracy") %>%
                                        select(c("Training", "Testing")), 2, mean)

# Accuracy SD in both the training and testing data
SD_Accuracy_Agnostic_genes <-  apply(All_XGB_CV %>%
                                       filter(model_type == "Agnostic_genes") %>%
                                       filter(metric == "Accuracy") %>%
                                       select(c("Training", "Testing")), 2, sd)

### Balanced Accuracy
# Mean Balanced Accuracy in both the training and testing data
MeanBalAccuracy_Agnostic_genes <-  apply(All_XGB_CV %>%
                                           filter(model_type == "Agnostic_genes") %>%
                                           filter(metric == "Bal.Accuracy") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Balanced Accuracy SD in both the training and testing data
SD_BalAccuracy_Agnostic_genes <-  apply(All_XGB_CV %>%
                                          filter(model_type == "Agnostic_genes") %>%
                                          filter(metric == "Bal.Accuracy") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Sensitivity
# Mean Sensitivity in both the training and testing data
MeanSensitivity_Agnostic_genes <-  apply(All_XGB_CV %>%
                                           filter(model_type == "Agnostic_genes") %>%
                                           filter(metric == "Sensitivity") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Sensitivity SD in both the training and testing data
SD_Sensitivity_Agnostic_genes <-  apply(All_XGB_CV %>%
                                          filter(model_type == "Agnostic_genes") %>%
                                          filter(metric == "Sensitivity") %>%
                                          select(c("Training", "Testing")), 2, sd)

### Specificity
# Mean Specificity in both the training and testing data
MeanSpecificity_Agnostic_genes <-  apply(All_XGB_CV %>%
                                           filter(model_type == "Agnostic_genes") %>%
                                           filter(metric == "Specificity") %>%
                                           select(c("Training", "Testing")), 2, mean)

# Specificity SD in both the training and testing data
SD_Specificity_Agnostic_genes <-  apply(All_XGB_CV %>%
                                          filter(model_type == "Agnostic_genes") %>%
                                          filter(metric == "Specificity") %>%
                                          select(c("Training", "Testing")), 2, sd)


### MCC
# Mean MCC in both the training and testing data
MeanMCC_Agnostic_genes <-  apply(All_XGB_CV %>%
                                   filter(model_type == "Agnostic_genes") %>%
                                   filter(metric == "MCC") %>%
                                   select(c("Training", "Testing")), 2, mean)

# MCC SD in both the training and testing data
SD_MCC_Agnostic_genes <-  apply(All_XGB_CV %>%
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

################

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
MeansAgnostic$Algorithm <- rep("XGB", 18)

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
MeansMechanistic$Algorithm <- rep("XGB", 18)


XGB_CV_AllMetrics <- rbind(MeansAgnostic, MeansMechanistic)

XGB_CV_AllMetrics <- XGB_CV_AllMetrics %>%  
  pivot_longer(cols = c(1,2), names_to = c("data_type"))

XGB_CV_AllMetrics$type <- XGB_CV_AllMetrics$metric
XGB_CV_AllMetrics$metric <- gsub("Mean", "", XGB_CV_AllMetrics$metric)
XGB_CV_AllMetrics$metric <- gsub("_Agnostic", "", XGB_CV_AllMetrics$metric)
XGB_CV_AllMetrics$metric <- gsub("_Mechanistic", "", XGB_CV_AllMetrics$metric)
XGB_CV_AllMetrics$metric <- gsub("_CILow", "", XGB_CV_AllMetrics$metric)
XGB_CV_AllMetrics$metric <- gsub("_CIHigh", "", XGB_CV_AllMetrics$metric) 

XGB_CV_AllMetrics$type[grep("Mean", XGB_CV_AllMetrics$type)] <- "Mean" 
XGB_CV_AllMetrics$type[grep("CILow", XGB_CV_AllMetrics$type)] <- "CILow" 
XGB_CV_AllMetrics$type[grep("CIHigh", XGB_CV_AllMetrics$type)] <- "CIHigh" 

XGB_CV_AllMetrics$metric[XGB_CV_AllMetrics$metric == "Sensitvity"] <- "Sensitivity"

XGB_CV_AllMetrics <- XGB_CV_AllMetrics %>% pivot_wider(names_from = "type", values_from = c("value"))

XGB_CV_AllMetrics$data_type <- factor(XGB_CV_AllMetrics$data_type, levels = c("Training", "Testing"))  
XGB_CV_AllMetrics$model_type <- factor(XGB_CV_AllMetrics$model_type, levels = c("Mechanistic", "Agnostic"))
XGB_CV_AllMetrics$metric <- factor(XGB_CV_AllMetrics$metric, levels = c("AUC", "Accuracy", "BalAccuracy", "Sensitivity", "Specificity", "MCC"))

## Plot the figure
png(filename = "./Figs/XGB/XGBCV_ALlMetrics.png", width = 3000, height = 2000, res = 300)
ggplot(data = XGB_CV_AllMetrics, aes(x = Mean, y = model_type, shape =  model_type, fill = model_type)) +
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
save(XGB_CV_AllMetrics, file = "./Objs/XGB_CV_AllMetrics.rda")

