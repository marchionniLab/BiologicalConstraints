###################################################################################
### Mohamed Omar
### April 2, 2020
### ### Goal: Creating the restricted pairs
### Adhesion genes
#################################################################################

###### 
# Clean Work space
rm(list = ls())
# Set work directory
#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

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
#load("./Objs/myTSPs.rda")
load("./Objs/LipogenesisPairs.rda")

lipogenesisPairs[, c(1,2)] <- lipogenesisPairs[, c(2,1)]
#colnames(lipogenesisPairs) <- colnames(myTSPs)

myTSPs <- lipogenesisPairs

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))
#keepGns_TF_MiR <- keepGns
#save(keepGns_TF_MiR, file = "./Objs/KTSP/KeepGns_TF_MiR.rda")

# restricted to the biological genes
usedTrainMat <- usedTrainMat[keepGns, ]
usedTestMat <- usedTestMat[keepGns, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k   [[500 without restrictedPairs -- 57 with restrictedPairs ]]
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=100,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs, disjoint = F)

ktspPredictorRes

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

KTSP_STATs_Train_Mechanistic <- t(ktspStatsTrainRes$comparisons)
KTSP_STATs_Train_Mechanistic[KTSP_STATs_Train_Mechanistic == FALSE] <- 0

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)


KTSP_STATs_Test_Mechanistic <- t(ktspStatsTestRes$comparisons)
KTSP_STATs_Test_Mechanistic[KTSP_STATs_Test_Mechanistic == FALSE] <- 0

save(KTSP_STATs_Train_Mechanistic, KTSP_STATs_Test_Mechanistic, file = "./Objs/KTSP/KTSP_STATs_Mechanistic_Lipogenesis.rda")

