###### 
# Clean Work space
rm(list = ls())


setwd("/Volumes/Macintosh/Research/Projects/ICB_Response")

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
library(preprocessCore)
######################################################
load("./Objs/ProcessedDatasets.rda")


# Get the Immune checkpoint genes 
ICgenes_Inhibitor <- c("BTLA", "CD200", "CD200R1", "CD274", 
                       "CD276", "CD80", "CD86", "CEACAM1", "CTLA4", "HAVCR2",
                       "IDO1", "LAG3", "PDCD1", "PDCD1LG2",
                       "PVR", "TIGIT") 

ICgenes_Activator <- c("CD27", "CD28", "CD40",
                       "ICOS", "IL2RB", 
                       "TNFRSF14", "TNFRSF18", "TNFRSF4", 
                       "TNFRSF9","TNFSF14", "TNFSF4", "TNFSF9")

IC_All <- c(ICgenes_Inhibitor, ICgenes_Activator)

keepGns <- intersect(IC_All, rownames(Expr_NP))

#Expr_NP <- Expr_NP[keepGns, ]
#Melanoma_Mat <- Melanoma_Mat[keepGns, ]

myTSPs <- expand.grid(ICgenes_Inhibitor, ICgenes_Activator)
myTSPs <- as.matrix(myTSPs)
#myTSPs <- t(combn(keepGns, 2))

keepGns2 <- c("PDCD1", "CD274", "CTLA4", "CD28", "CD80", "CD86")

myTSPs <- myTSPs[myTSPs[,1] %in% keepGns2 | myTSPs[,2] %in% keepGns2 , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #12
featNo <- nrow(Expr_NP)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  Expr_NP, Pheno_NP$SR, krange=ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, 
  RestrictedPairs = myTSPs, disjoint = F)

ktspPredictorRes




############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = Expr_NP, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

#KTSP_STATs_Train_Mechanistic <- t(ktspStatsTrainRes$comparisons)
#KTSP_STATs_Train_Mechanistic[KTSP_STATs_Train_Mechanistic == FALSE] <- 0

### Threshold
thr <- coords(roc(Pheno_NP$SR, ktspStatsTrainRes$statistics, levels = c("no", "yes"), direction = "<"), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(Pheno_NP$SR, ktspStatsTrainRes$statistics, levels = c("no", "yes"), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(Pheno_NP$SR, ktspStatsTrainRes$statistics, plot = T, print.thres="all", print.auc=TRUE, print.auc.col="black", levels = c("no", "yes"), direction = "<", col="black", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(Expr_NP, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, Pheno_NP$SR, positive = "yes", mode = "everything")

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = Pheno_NP$SR)
MCC_Mechanistic_Train


###########################################3
# Test on all Melanoma
ktspPredictorRes$labels <- c("NR", "R")

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Melanoma_Mat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Melanoma_Group, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Melanoma_Mat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Melanoma_Group, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Melanoma_Group)
MCC_Mechanistic_Test

