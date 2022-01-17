###################################################################################
### Mohamed Omar
### 5/5/2019
### ### Goal: Creating the restricted ktsp classifier.
### Androgen genes
#################################################################################

###### 
# Clean Work space
rm(list = ls())
# Set work directory
setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

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
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

############################################################################
## Load the selected genes
Coverings <- load("/Volumes/Macintosh/Research/Genes/coverings.core_union.rda")

ProstateTargetUnion <- list.target$Prostate$union

myTSPs <- t(combn(ProstateTargetUnion, 2))

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))
#keepGns_TF_MiR <- keepGns
#save(keepGns_TF_MiR, file = "./Objs/KTSP/KeepGns_TF_MiR.rda")

#usedTrainMat <- usedTrainMat[keepGns, ]
#usedTestMat <- usedTestMat[keepGns, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

#print(xtable(myTSPs, type = "latex"), file = "./Objs/KTSP/Restricted_Pairs.tex")
#write.csv(myTSPs, file = "./Objs/KTSP/Restricted_Pairs.csv")

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #20
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

ktspPredictorRes

#save(ktspPredictorRes, file = "./Objs/KTSP/MechKTSP94Pairs.rda")

#Mechanistic_KTSP <- cbind(ktspPredictorRes$TSPs, ktspPredictorRes$score)
#colnames(Mechanistic_KTSP) <- c("gene1", "gene2", "score")

#write.csv(Mechanistic_KTSP, file = "./Objs/KTSP/MechanisticKTSP.csv")
#print(xtable(Mechanistic_KTSP, type = "latex"), file = "./Objs/KTSP/Mechanistic.tex")

## Save the mechanistic Pairs
#MechanisticKTSP_Pairs <- c("COX7B>EZH2", "RAD21>ARHGEF11", "GTF2B>SF3B4", "TP53>MTMR4", "BTG2>CHD4", "MZF1>PITX1", "NCOR1>TRIP13")
#save(MechanisticKTSP_Pairs, file = "./Objs/KTSP/MechanisticKTSP_Pairs.rda")

### Check how many TSP
#xx <- SWAP.Filter.Wilcoxon(usedTrainGroup,usedTrainMat,featureNo=featNo)
#dim(myTSPs[myTSPs[,1] %in% xx & myTSPs[,2] %in% xx ,])
#choose(78,2)

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Mets", mode = "everything")

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = usedTrainGroup)
MCC_Mechanistic_Train

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Mets", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Mechanistic_Test

## OR (calculate accuracy, sens, spec and AUC)
SWAP.GetKTSP.PredictionStats(usedTestPredictionRes, usedTestGroup, decision_values = ktspStatsTestRes$statistics)

### basic plots of various metrics
library(precrec)
precrec_obj <- evalmod(scores = ktspStatsTestRes$statistics, labels = usedTestGroup, mode="basic")
autoplot(precrec_obj)   

