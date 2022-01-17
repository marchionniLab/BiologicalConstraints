################################################################################
### Mohamed Omar
### April 2, 2020 
### Goal: Creating unrestricted Pairs
###############################################################################

rm(list = ls()) 

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

####################### 
##Load required packages
library(switchBox)
library(Biobase)
library(limma)
library(pROC)
library(caret)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(plotROC)
library(xtable)
library(mltools)

#####################################################################
### Load data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")
##########################################################
## Quantile normalize the datasets
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

###########################################################################
### TRAINING using all expressed genes

## Set Feature number and max K
featN <- nrow(usedTrainMat) # the same as in the mechanistic classifier 

###########
### Train a classifier using the default filter function
set.seed(333)

ktspPredictorUnRes_500 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = 500, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= featN)
ktspPredictorUnRes_500

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_500 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_500, CombineFunc = sum)
summary(ktspStatsTrainUnRes_500$statistics)


KTSP_STATs_Train_Agnostic_500 <- t(ktspStatsTrainUnRes_500$comparisons)
KTSP_STATs_Train_Agnostic_500[KTSP_STATs_Train_Agnostic_500 == FALSE] <- 0

################################################################################
###############################################################################
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_500 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_500, CombineFunc = sum)
summary(ktspStatsTestUnRes_500$statistics)


KTSP_STATs_Test_Agnostic_500 <- t(ktspStatsTestUnRes_500$comparisons)
KTSP_STATs_Test_Agnostic_500[KTSP_STATs_Test_Agnostic_500 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_500, KTSP_STATs_Test_Agnostic_500, file = "./Objs/KTSP/KTSP_STATs_Agnostic_Combined_500.rda")

