################################################################################
### Mohamed Omar
### 10/4/2019
### Goal: Creating unrestricted K-TSP classifier
### Cross-study validation : JHU out
###############################################################################

rm(list = ls())

setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

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
load("./Objs/LOO/MetastasisData_JHUOut.rda")
load("./Objs/Correlation/RGenes.rda")

##########################################################
## Quantile normalize the datasets
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")[RGenes, ]
usedTestMat <- testMat[RGenes, ]

usedTestMat <- t(scale(t(usedTestMat), center = F, scale = TRUE))

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

###########################################################################
### TRAINING using all expressed genes

## Set Feature number and max K
featN <- nrow(usedTrainMat) # the same as in the mechanistic classifier 
ktsp <- c(3:25)  # the same as in the mechanistic classifier


###########
### Train a classifier using the default filter function
ktspPredictorUnRes_500Feat <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = ktsp, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 500)
ktspPredictorUnRes_500Feat

#Agnostic_KTSP <- cbind(ktspPredictorUnRes$TSPs, ktspPredictorUnRes$score)
#colnames(Agnostic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Agnostic_KTSP, type = "latex"), file = "./Objs/KTSP/Agnostic.tex")

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_500Feat <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_500Feat, CombineFunc = sum)
summary(ktspStatsTrainUnRes_500Feat$statistics)

## Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes_500Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

## Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes_500Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: reorder the levels ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain_500Feat <- roc(usedTrainGroup, ktspStatsTrainUnRes_500Feat$statistics, plot = FALSE, ci = TRUE, print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")
ROCTrain_500Feat

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes_500Feat <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_500Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionTrain_500Feat <- confusionMatrix(usedTrainPredictionUnRes_500Feat, usedTrainGroup, positive = "Mets")
confusionTrain_500Feat

MCC_Agnostic_Train_500Feat <- mltools::mcc(pred = usedTrainPredictionUnRes_500Feat, actuals = usedTrainGroup)
MCC_Agnostic_Train_500Feat

# Put the performance metrics together
TrainPerf_500Feat <- data.frame("Training" = c(ROCTrain_500Feat$ci, confusionTrain_500Feat$overall["Accuracy"], confusionTrain_500Feat$byClass["Balanced Accuracy"], confusionTrain_500Feat$byClass["Sensitivity"], confusionTrain_500Feat$byClass["Specificity"], MCC_Agnostic_Train_500Feat))
TrainPerf_500Feat[1:3, ] <- TrainPerf_500Feat[c(2,1,3), ]
rownames(TrainPerf_500Feat) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#################################################################################
###############################################################################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_500Feat <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_500Feat, CombineFunc = sum)
summary(ktspStatsTestUnRes_500Feat$statistics)

# plot the curve
ROCTest_500Feat <- roc(usedTestGroup, ktspStatsTestUnRes_500Feat$statistics, plot = F, ci = T, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")
ROCTest_500Feat

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes_500Feat <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes_500Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionTest_500Feat <- confusionMatrix(usedTestPredictionUnRes_500Feat, usedTestGroup, positive = "Mets", mode = "everything")
confusionTest_500Feat

MCC_Agnostic_Test_500Feat <- mltools::mcc(pred = usedTestPredictionUnRes_500Feat, actuals = usedTestGroup)
MCC_Agnostic_Test_500Feat

## Save the performance metrics
TestPerf_500Feat <- data.frame("Testing" = c(ROCTest_500Feat$ci, confusionTest_500Feat$overall["Accuracy"], confusionTest_500Feat$byClass["Balanced Accuracy"], confusionTest_500Feat$byClass["Sensitivity"], confusionTest_500Feat$byClass["Specificity"], MCC_Agnostic_Test_500Feat))
TestPerf_500Feat[1:3, ] <- TestPerf_500Feat[c(2,1,3), ]
rownames(TestPerf_500Feat) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
JHU_Out_AgnosticPerformance_500Feat <- cbind(TrainPerf_500Feat, TestPerf_500Feat)

# Save
save(JHU_Out_AgnosticPerformance_500Feat, file = "./Objs/KTSP/JHU_Out_AgnosticPerformance_500Feat.rda")
