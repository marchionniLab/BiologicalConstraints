############################################################################


rm(list = ls())

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/BreastChemo")

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
library(mltools)
library(xtable)
library(dplyr)
library(Ckmeans.1d.dp)


## Load data
load("./Objs/ChemoDataNew.rda")

## Load the selected genes (TFs-Targets)
load("./Objs/NotchPairs.rda")
load("./Objs/MycPairs.rda")


myTSPs <- rbind(NotchPairs, MycPairs)
colnames(myTSPs) <- c("BadGene", "GoodGene")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))

usedTrainMat <- usedTrainMat[keepGns, ]
UsedTestMat <- usedTestMat[keepGns, ]

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=241, featureNo= featNo, 
  FilterFunc = SWAP.Filter.Wilcoxon, RestrictedPairs = myTSPs)

ktspPredictorRes

# # Filter the mechanistic KTSP based on biology
# mech_keep <- ktspPredictorRes$TSPs[,1] %in% unique(myTSPs[,"BadGene"]) & ktspPredictorRes$TSPs[,2] %in% unique(myTSPs[,"GoodGene"])
# table(mech_keep)
# 
# ### Subset
# ktspPredictorRes$score <- ktspPredictorRes$score[mech_keep]
# ktspPredictorRes$TSPs <- ktspPredictorRes$TSPs[mech_keep, ]
# ktspPredictorRes$tieVote <- droplevels(ktspPredictorRes$tieVote[mech_keep])

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



save(KTSP_STATs_Train_Mechanistic, KTSP_STATs_Test_Mechanistic, file = "./Objs/KTSP/TNBC_KTSP_STATs_Mechanistic_NotchAndMYC2.rda")

