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
## load mechanistic pairs
load("./Objs/StromaPairs.rda")
StromaPairs[,c(1,2)] <- StromaPairs[,c(2,1)]

load("./Objs/LipogenesisPairs.rda")
lipogenesisPairs[,c(1,2)] <- lipogenesisPairs[,c(2,1)]

# Make mechanistic pairs

#BadGenes <- c("FDFT1", "FASN", "SREBF1", "DHCR7", "LSS", "HMGCS1", "DHCR24", "SCD", "FADS2", "SQLE", "IDI1", "CYP51A1", "SC5D", "HMGCR", "FDPS", "MVD", "MVK", "NSDHL", "ACSS2", "TM7SF2", "HSD17B7", "MSMO1", "FADS3", "ACSL1", "SREBF2", "ACACA", "FADS1", "ELOVL5", "ACSL3", "ME1", "HSD17B12", "ELOVL7", "FAXDC2", "ELOVL4", "ACSL4", "PTPLAD1", "PTPLB", "ACOT4")

#GoodGenes <- c("ELOVL3", "ACLY", "TECR", "ACOT1", "PTPLA", "OLAH", "SCD5", "MLYCD", "ACSL6", "ACOT2", "ACSBG2", "ACSBG1", "FADS6", "ELOVL2", "ACSL5", "IDI2", "EBP", "ACAT2", "ACOT7", "ELOVL1", "ELOVL6", "PMVK", "MCAT", "OXSM", "ME2")


#Lipogenesis <- expand.grid(BadGenes, GoodGenes)
#Lipogenesis <- as.matrix(Lipogenesis)
#Lipogenesis[,c(1,2)] <- Lipogenesis[, c(2,1)]
#colnames(Lipogenesis) <- c("GoodGene", "BadGene")

myTSPs <- rbind(StromaPairs, lipogenesisPairs)
#myTSPs <- lipogenesisPairs
colnames(myTSPs) <- c("GoodGene", "BadGene")

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
  usedTrainMat, usedTrainGroup, krange=614,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs, disjoint = F)

ktspPredictorRes

### Check consistency with biology
keepTest <- ktspPredictorRes$TSPs[,1] %in% myTSPs[,"GoodGene"] & ktspPredictorRes$TSPs[,2] %in% myTSPs[,"BadGene"]

table(keepTest)

### Subset
ktspPredictorRes$score <- ktspPredictorRes$score[keepTest]
ktspPredictorRes$TSPs <- ktspPredictorRes$TSPs[keepTest, ]
#ktspPredictorRes$tieVote  <- ktspPredictorRes$tieVote[!is.na(ktspPredictorRes$tieVote)][keepTest]
ktspPredictorRes$tieVote <- droplevels(ktspPredictorRes$tieVote[keepTest])

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

save(KTSP_STATs_Train_Mechanistic, KTSP_STATs_Test_Mechanistic, file = "./Objs/KTSP/KTSP_STATs_Mechanistic_StromaANDlipogenesis.rda")

