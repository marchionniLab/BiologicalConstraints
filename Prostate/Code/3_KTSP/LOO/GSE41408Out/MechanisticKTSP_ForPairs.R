###################################################################################
### Mohamed Omar
### 25/11/2019
### ### Goal: Creating the restricted ktsp classifier.
### TF-Target genes
## Leaving 1 dataset out (GSE41408)
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
load("./Objs/LOO/MetastasisData_GSE41408Out.rda")
load("./Objs/Correlation/RGenes.rda")

############################################################################
## Load the selected genes
Genes1 <- read.delim("./geneset.txt")
Genes1 <- as.matrix(Genes1)
Genes1 <- Genes1[-1,]

Genes2 <- read.delim("./geneset2.txt")
Genes2 <- as.matrix(Genes2)
Genes2 <- Genes2[-1,]

Genes3 <- read.delim("./geneset3.txt")
Genes3 <- as.matrix(Genes3)
Genes3 <- Genes3[-1,]

Genes <- c(Genes1,Genes2, Genes3)
Genes <- Genes[!duplicated(Genes)]

myTSPs <- t(combn(Genes,2))

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")[RGenes,]
usedTestMat <- testMat[RGenes, ]

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
ktsp <- c(3:25) #70
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=100,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

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

save(KTSP_STATs_Train_Mechanistic, KTSP_STATs_Test_Mechanistic, file = "./Objs/KTSP/LOO/KTSP_STATs_Mechanistic_GSE41408Out.rda")


