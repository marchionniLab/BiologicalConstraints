################################################################################
### Mohamed Omar
### 10/4/2019
### Goal: Creating unrestricted K-TSP classifier
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
ktsp <- c(3:25)  # the same as in the mechanistic classifier


###########
### Agnostic KTSP All Features
set.seed(333)

ktspPredictorUnRes_AllFeat <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = ktsp, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo = featN)
ktspPredictorUnRes_AllFeat

#Agnostic_KTSP <- cbind(ktspPredictorUnRes$TSPs, ktspPredictorUnRes$score)
#colnames(Agnostic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Agnostic_KTSP, type = "latex"), file = "./Objs/KTSP/Agnostic.tex")

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_AllFeat <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_AllFeat, CombineFunc = sum)
summary(ktspStatsTrainUnRes_AllFeat$statistics)

## Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes_AllFeat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

## Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes_AllFeat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: reorder the levels ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainUnRes_AllFeat$statistics, plot = TRUE, print.auc = TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes_AllFeat <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_AllFeat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionUnRes_AllFeat, usedTrainGroup, positive = "Mets")

MCC_Train <- mltools::mcc(pred = usedTrainPredictionUnRes_AllFeat, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, ConfusionTrain$overall["Accuracy"], ConfusionTrain$byClass["Balanced Accuracy"], ConfusionTrain$byClass["Sensitivity"], ConfusionTrain$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")


################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_AllFeat <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_AllFeat, CombineFunc = sum)
summary(ktspStatsTestUnRes_AllFeat$statistics)

# plot the curve
ROCTest <- roc(usedTestGroup, ktspStatsTestUnRes_AllFeat$statistics, plot = F, print.thres = thr, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes_AllFeat <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes_AllFeat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionUnRes_AllFeat, usedTestGroup, positive = "Mets", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionUnRes_AllFeat, actuals = usedTestGroup)
MCC_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
AgnosticKTSP_AllFeat_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(AgnosticKTSP_AllFeat_Perf, file = "./Objs/KTSP/AgnosticKTSP_AllFeat_Perf.rda")

####################################################################################
#####################################################################################
#################################################################################3
### Agnostic KTSP top 50 Features
set.seed(333)

ktspPredictorUnRes_50Feat <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                              krange = ktsp, 
                                              FilterFunc = SWAP.Filter.Wilcoxon, featureNo = 50)
ktspPredictorUnRes_50Feat

#Agnostic_KTSP <- cbind(ktspPredictorUnRes$TSPs, ktspPredictorUnRes$score)
#colnames(Agnostic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Agnostic_KTSP, type = "latex"), file = "./Objs/KTSP/Agnostic.tex")

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_50Feat <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_50Feat, CombineFunc = sum)
summary(ktspStatsTrainUnRes_50Feat$statistics)

## Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes_50Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

## Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes_50Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: reorder the levels ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainUnRes_50Feat$statistics, plot = TRUE, print.auc = TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes_50Feat <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_50Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionUnRes_50Feat, usedTrainGroup, positive = "Mets")

MCC_Train <- mltools::mcc(pred = usedTrainPredictionUnRes_50Feat, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, ConfusionTrain$overall["Accuracy"], ConfusionTrain$byClass["Balanced Accuracy"], ConfusionTrain$byClass["Sensitivity"], ConfusionTrain$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")


###################


#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_50Feat <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_50Feat, CombineFunc = sum)
summary(ktspStatsTestUnRes_50Feat$statistics)

# plot the curve
ROCTest <- roc(usedTestGroup, ktspStatsTestUnRes_50Feat$statistics, plot = F, print.thres = thr, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes_50Feat <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes_50Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionUnRes_50Feat, usedTestGroup, positive = "Mets", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionUnRes_50Feat, actuals = usedTestGroup)
MCC_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
AgnosticKTSP_50Feat_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(AgnosticKTSP_50Feat_Perf, file = "./Objs/KTSP/AgnosticKTSP_50Feat_Perf.rda")


####################################################################################
#####################################################################################
#################################################################################3
### Agnostic KTSP top 100 Features
set.seed(333)

ktspPredictorUnRes_100Feat <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                             krange = ktsp, 
                                             FilterFunc = SWAP.Filter.Wilcoxon, featureNo = 100)
ktspPredictorUnRes_100Feat

#Agnostic_KTSP <- cbind(ktspPredictorUnRes$TSPs, ktspPredictorUnRes$score)
#colnames(Agnostic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Agnostic_KTSP, type = "latex"), file = "./Objs/KTSP/Agnostic.tex")

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_100Feat <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_100Feat, CombineFunc = sum)
summary(ktspStatsTrainUnRes_100Feat$statistics)

## Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes_100Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

## Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes_100Feat$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: reorder the levels ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainUnRes_100Feat$statistics, plot = TRUE, print.auc = TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes_100Feat <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_100Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionUnRes_100Feat, usedTrainGroup, positive = "Mets")

MCC_Train <- mltools::mcc(pred = usedTrainPredictionUnRes_100Feat, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, ConfusionTrain$overall["Accuracy"], ConfusionTrain$byClass["Balanced Accuracy"], ConfusionTrain$byClass["Sensitivity"], ConfusionTrain$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")


###################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_100Feat <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_100Feat, CombineFunc = sum)
summary(ktspStatsTestUnRes_100Feat$statistics)

# plot the curve
ROCTest <- roc(usedTestGroup, ktspStatsTestUnRes_100Feat$statistics, plot = F, print.thres = thr, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes_100Feat <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes_100Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionUnRes_100Feat, usedTestGroup, positive = "Mets", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionUnRes_100Feat, actuals = usedTestGroup)
MCC_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
AgnosticKTSP_100Feat_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(AgnosticKTSP_100Feat_Perf, file = "./Objs/KTSP/AgnosticKTSP_100Feat_Perf.rda")

####################################################################################
#####################################################################################
#################################################################################3
### Agnostic KTSP top 200 Features
set.seed(333)

ktspPredictorUnRes_200Feat <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                              krange = ktsp, 
                                              FilterFunc = SWAP.Filter.Wilcoxon, featureNo = 200)
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
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainUnRes_200Feat$statistics, plot = TRUE, print.auc = TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes_200Feat <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_200Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionUnRes_200Feat, usedTrainGroup, positive = "Mets")

MCC_Train <- mltools::mcc(pred = usedTrainPredictionUnRes_200Feat, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, ConfusionTrain$overall["Accuracy"], ConfusionTrain$byClass["Balanced Accuracy"], ConfusionTrain$byClass["Sensitivity"], ConfusionTrain$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")


###################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_200Feat <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_200Feat, CombineFunc = sum)
summary(ktspStatsTestUnRes_200Feat$statistics)

# plot the curve
ROCTest <- roc(usedTestGroup, ktspStatsTestUnRes_200Feat$statistics, plot = F, print.thres = thr, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes_200Feat <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes_200Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionUnRes_200Feat, usedTestGroup, positive = "Mets", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionUnRes_200Feat, actuals = usedTestGroup)
MCC_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
AgnosticKTSP_200Feat_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(AgnosticKTSP_200Feat_Perf, file = "./Objs/KTSP/AgnosticKTSP_200Feat_Perf.rda")

####################################################################################
#####################################################################################
#################################################################################3
### Agnostic KTSP top 500 Features
set.seed(333)

ktspPredictorUnRes_500Feat <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                              krange = ktsp, 
                                              FilterFunc = SWAP.Filter.Wilcoxon, featureNo = 500)
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
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainUnRes_500Feat$statistics, plot = TRUE, print.auc = TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes_500Feat <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_500Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionUnRes_500Feat, usedTrainGroup, positive = "Mets")

MCC_Train <- mltools::mcc(pred = usedTrainPredictionUnRes_500Feat, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, ConfusionTrain$overall["Accuracy"], ConfusionTrain$byClass["Balanced Accuracy"], ConfusionTrain$byClass["Sensitivity"], ConfusionTrain$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")


###################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_500Feat <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_500Feat, CombineFunc = sum)
summary(ktspStatsTestUnRes_500Feat$statistics)

# plot the curve
ROCTest <- roc(usedTestGroup, ktspStatsTestUnRes_500Feat$statistics, plot = F, print.thres = thr, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes_500Feat <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes_500Feat, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionUnRes_500Feat, usedTestGroup, positive = "Mets", mode = "everything")
ConfusionTest

MCC_Test <- mltools::mcc(pred = usedTestPredictionUnRes_500Feat, actuals = usedTestGroup)
MCC_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
AgnosticKTSP_500Feat_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(AgnosticKTSP_500Feat_Perf, file = "./Objs/KTSP/AgnosticKTSP_500Feat_Perf.rda")

