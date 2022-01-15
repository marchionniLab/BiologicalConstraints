
###### 
# Clean Work space
rm(list = ls())

############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(enrichR)
library(mltools)
library(xtable)

###########################################################################
### Load expression and phenotype data
load("./Objs/ChemoDataNew_GSE103668Out.rda")

############################################################################
## Load the selected genes *(Notch pairs)
load("./Objs/NotchPairs.rda")
load("./Objs/MycPairs.rda")

myTSPs <- rbind(NotchPairs, MycPairs)
colnames(myTSPs) <- c("BadGene", "GoodGene")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")
usedTestMat <- testMat

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))

usedTrainMat <- usedTrainMat[keepGns, ]
usedTestMat <- usedTestMat[keepGns, ]

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #7
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

ktspPredictorRes

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("Sensitive", "Resistant"), direction = "<" ), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("Sensitive", "Resistant"), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, print.thres.adj=c(0.01,1.25), print.auc=TRUE, ci = T, print.auc.col="black", levels = c("Sensitive", "Resistant"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")
ROCTrain

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Resistant", mode = "everything")

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = usedTrainGroup)
MCC_Mechanistic_Train

confusionTrain <- confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Resistant", mode = "everything")
confusionTrain

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, confusionTrain$overall["Accuracy"], confusionTrain$byClass["Balanced Accuracy"], confusionTrain$byClass["Sensitivity"], confusionTrain$byClass["Specificity"], MCC_Mechanistic_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
ROCTest <- roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("Sensitive", "Resistant"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROCTest

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionTest <- confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Resistant", mode = "everything")
confusionTest

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Mechanistic_Test


TestPerf <- data.frame("Testing" = c(ROCTest$ci, confusionTest$overall["Accuracy"], confusionTest$byClass["Balanced Accuracy"], confusionTest$byClass["Sensitivity"], confusionTest$byClass["Specificity"], MCC_Mechanistic_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
GSE103668_Out_MechPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(GSE103668_Out_MechPerformance, file = "./Objs/KTSP/GSE103668_Out_MechPerformance.rda")
