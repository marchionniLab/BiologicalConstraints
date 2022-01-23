

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

######################################################
## Load the Processed datasets
load("./Objs/ProcessedDatasets.rda")


# Get the Immune checkpoint genes 
# ICgenes_Inhibitor <- c("ADORA2A", "BTLA", "VSIR", "CD200", "CD200R1", "CD274", 
#                        "CD276", "CD80", "CD86", "CEACAM1", "CTLA4", "GAL3", "TIM3",
#                        "IDO1", "KIR3DL1", "LAG3", "LAIR1", "PDCD1", "PDCD1LG2",
#                        "PVR", "PVRL2", "TIGIT", "VTCN1") 
# 
# ICgenes_Activator <- c("CD266", "CD27", "CD28", "CD40", "CD40LG", "CD70LG",
#                        "DR3", "HAVCR1", "ICOS", "ICOSL", "IL2RB", "NAIL", 
#                        "SLAMF1", "TIM2", "TNFRSF14", "TNFRSF18", "TNFRSF4", 
#                        "TNFRSF9","TNFSF14", "TNFSF18", "TNFSF4", "TNFSF9")

# Get the Immune checkpoint genes 
ICgenes_Inhibitor <- c("BTLA", "CD200", "CD200R1", "CD274", 
                       "CD276", "CD80", "CD86", "CEACAM1", "CTLA4", "HAVCR1",
                       "IDO1", "LAG3", "PDCD1", "PDCD1LG2",
                       "PVR", "TIGIT") 

ICgenes_Activator <- c("CD27", "CD28", "CD40",
                       "ICOS", "IL2RB", 
                       "TNFRSF14", "TNFRSF18", "TNFRSF4", 
                       "TNFRSF9","TNFSF14", "TNFSF4", "TNFSF9")


IC_All <- c(ICgenes_Inhibitor, ICgenes_Activator)

keepGns <- intersect(IC_All, rownames(Expr_NP))

myTSPs <- expand.grid(ICgenes_Inhibitor, ICgenes_Activator)
myTSPs <- as.matrix(myTSPs)
#myTSPs <- t(combn(keepGns, 2))

keepGns2 <- c("PDCD1", "CD274", "CTLA4", "CD28", "CD80", "CD86")

myTSPs <- myTSPs[myTSPs[,1] %in% keepGns2 | myTSPs[,2] %in% keepGns2 , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #8
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

#########################################################################
#########################################################################
### Testing

# Test on Chen_Pre_aCTLA4

## Compute the sum and find the best threshold
ktspPredictorRes$labels <- c("NR", "R")

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_Chen_Pre_aCTLA4, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_Chen_Pre_aCTLA4$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_Chen_Pre_aCTLA4, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_Chen_Pre_aCTLA4$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_Chen_Pre_aCTLA4$Response)
MCC_Mechanistic_Test

###########################################3
# Test on Chen_Pre_aPD1

## Compute the sum and find the best threshold

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_Chen_Pre_aPD1, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_Chen_Pre_aPD1$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_Chen_Pre_aPD1, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_Chen_Pre_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_Chen_Pre_aPD1$Response)
MCC_Mechanistic_Test


###########################################3
# Test on Chen_On_aPD1

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_Chen_On_aPD1, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_Chen_On_aPD1$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_Chen_On_aPD1, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_Chen_On_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_Chen_On_aPD1$Response)
MCC_Mechanistic_Test

###########################################3
# Test on PRAT

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_PRAT, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_PRAT$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_PRAT, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_PRAT$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_PRAT$Response)
MCC_Mechanistic_Test



###########################################3
# Test on HUGO

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_HUGO), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_HUGO$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_HUGO), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_HUGO$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_HUGO$Response)
MCC_Mechanistic_Test

###########################################3
# Test on TCGA

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_TCGA), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_TCGA$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_TCGA), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_TCGA$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_TCGA$Response)
MCC_Mechanistic_Test

####################################################
# Test on RIAZ_Pre_aPD1

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_RIAZ_Pre_aPD1), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_RIAZ_Pre_aPD1$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_RIAZ_Pre_aPD1), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_RIAZ_Pre_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_RIAZ_Pre_aPD1$Response)
MCC_Mechanistic_Test

####################################################
# Test on RIAZ_On_aPD1

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_RIAZ_On_aPD1), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_RIAZ_On_aPD1$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_RIAZ_On_aPD1), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_RIAZ_On_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_RIAZ_On_aPD1$Response)
MCC_Mechanistic_Test

####################################################
# Test on MGH_On_aPD1

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_MGH_On_aPD1), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_MGH_On_aPD1$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_MGH_On_aPD1), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_MGH_On_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_MGH_On_aPD1$Response)
MCC_Mechanistic_Test

####################################################
# Test on MGH_On_aCTLA4

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_MGH_On_aCTLA4), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_MGH_On_aCTLA4$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_MGH_On_aCTLA4), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_MGH_On_aCTLA4$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_MGH_On_aCTLA4$Response)
MCC_Mechanistic_Test

####################################################
# Test on VanAllen

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_VanAllen), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_VanAllen$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_VanAllen), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_VanAllen$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_VanAllen$Response)
MCC_Mechanistic_Test



###########################################3
# Test on all Melanoma

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




