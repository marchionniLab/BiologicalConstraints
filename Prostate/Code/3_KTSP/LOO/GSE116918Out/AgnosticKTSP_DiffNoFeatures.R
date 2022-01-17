################################################################################
### Mohamed Omar
### 10/4/2019
### Goal: Creating unrestricted K-TSP classifier
### Cross-study validation : GSE116918 out
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
load("./Objs/LOO/MetastasisData_GSE116918Out.rda")
load("./Objs/Correlation/RGenes.rda")

##########################################################
## Quantile normalize the datasets
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")[RGenes, ]
usedTestMat <- testMat[RGenes, ]

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
ktspPredictorUnRes_200Feat <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = ktsp, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 200)
ktspPredictorUnRes_200Feat

#Agnostic_KTSP <- cbind(ktspPredictorUnRes$TSPs, ktspPredictorUnRes$score)
#colnames(Agnostic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Agnostic_KTSP, type = "latex"), file = "./Objs/KTSP/Agnostic.tex")

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_200Feat <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_200Feat, CombineFunc = sum)
summary(ktspStatsTrainUnRes_200Feat$statistics)

## Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes_200Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

## Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes_200Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: reorder the levels ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain_200Feat <- roc(usedTrainGroup, ktspStatsTrainUnRes_200Feat$statistics, plot = FALSE, ci = TRUE, print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")
ROCTrain_200Feat

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes_200Feat <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_200Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionTrain_200Feat <- confusionMatrix(usedTrainPredictionUnRes_200Feat, usedTrainGroup, positive = "Mets")
confusionTrain_200Feat

MCC_Agnostic_Train_200Feat <- mltools::mcc(pred = usedTrainPredictionUnRes_200Feat, actuals = usedTrainGroup)
MCC_Agnostic_Train_200Feat

# Put the performance metrics together
TrainPerf_200Feat <- data.frame("Training" = c(ROCTrain_200Feat$ci, confusionTrain_200Feat$overall["Accuracy"], confusionTrain_200Feat$byClass["Balanced Accuracy"], confusionTrain_200Feat$byClass["Sensitivity"], confusionTrain_200Feat$byClass["Specificity"], MCC_Agnostic_Train_200Feat))
TrainPerf_200Feat[1:3, ] <- TrainPerf_200Feat[c(2,1,3), ]
rownames(TrainPerf_200Feat) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#################################################################################
###############################################################################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_200Feat <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_200Feat, CombineFunc = sum)
summary(ktspStatsTestUnRes_200Feat$statistics)

# plot the curve
ROCTest_200Feat <- roc(usedTestGroup, ktspStatsTestUnRes_200Feat$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", ci = T, col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")
ROCTest_200Feat

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes_200Feat <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes_200Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionTest_200Feat <- confusionMatrix(usedTestPredictionUnRes_200Feat, usedTestGroup, positive = "Mets", mode = "everything")
confusionTest_200Feat

MCC_Agnostic_Test_200Feat <- mltools::mcc(pred = usedTestPredictionUnRes_200Feat, actuals = usedTestGroup)
MCC_Agnostic_Test_200Feat

## Save the performance metrics
TestPerf_200Feat <- data.frame("Testing" = c(ROCTest_200Feat$ci, confusionTest_200Feat$overall["Accuracy"], confusionTest_200Feat$byClass["Balanced Accuracy"], confusionTest_200Feat$byClass["Sensitivity"], confusionTest_200Feat$byClass["Specificity"], MCC_Agnostic_Test_200Feat))
TestPerf_200Feat[1:3, ] <- TestPerf_200Feat[c(2,1,3), ]
rownames(TestPerf_200Feat) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
GSE116918_Out_AgnosticPerformance_200Feat <- cbind(TrainPerf_200Feat, TestPerf_200Feat)

# Save
save(GSE116918_Out_AgnosticPerformance_200Feat, file = "./Objs/KTSP/GSE116918_Out_AgnosticPerformance_200Feat.rda")

