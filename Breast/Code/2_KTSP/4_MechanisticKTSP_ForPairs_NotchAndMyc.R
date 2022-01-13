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
library(enrichR)
library(mltools)
library(xtable)
library(pdp)
library(DiagrammeR)
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

#Mechanistic_KTSP <- cbind(ktspPredictorRes$TSPs, ktspPredictorRes$score)
#colnames(Mechanistic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Mechanistic_KTSP, type = "latex"), file = "./Objs/KTSP/Mechanistic.tex")

## Save the mechanistic Pairs
#MechanisticKTSP_Pairs <- c("ZKSCAN1>ZW10", "HES1>NR3C1", "ZNF160>ZMIZ1", "KLHL9>PML", "NUDC>RELA", "SNCA>ZBTB7A", "ALDH7A1>NFIC", "SMC3>RNF44", "UBE2N>SUZ12", "ARL2>MAX", "EGR1>TSTA3")
#save(MechanisticKTSP_Pairs, file = "./Objs/KTSP/MechanisticKTSP_Pairs.rda")


###########################################################################
### Check consistency with biology
# keep <- ktspPredictorRes$TSPs[,1] %in% myTSPs[,1]
# table(keep)
# #
# # ###Subset
# ktspPredictorRes$name <- paste(sum(keep), "TSPs", sep="")
# ktspPredictorRes$TSPs <- ktspPredictorRes$TSPs[ keep, ]
# ktspPredictorRes$score <- ktspPredictorRes$score[ keep ]
# #
# # ### Visualize the classifier
# ktspPredictorRes


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

