###############################################################################
### Mohamed Omar
### 14/07/2019
### Goal : Creating the Agnostic K-TSP classifier for prostate cancer metastasis
#################################################################################

rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

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

###################################################################
## Load Metastasis Data
load("./Objs/MetastasisData_New.rda")


## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")
#boxplot(usedTrainMat, outline = FALSE)

usedTestMat <- normalizeBetweenArrays(testMat, method = "quantile")
#boxplot(usedTestMat, outline = FALSE)

#
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup




### Train a classifier using the default filter function
featN <- nrow(usedTrainMat)
ktsp <- c(3:25) #20

ktspPredictorUnRes <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = 500, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo = 1000)

ktspPredictorUnRes

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainUnRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTrainUnRes$statistics)

KTSP_STATs_Train_Agnostic <- t(ktspStatsTrainUnRes$comparisons)
KTSP_STATs_Train_Agnostic[KTSP_STATs_Train_Agnostic == FALSE] <- 0


#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestUnRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTestUnRes$statistics)

KTSP_STATs_Test_Agnostic <- t(ktspStatsTestUnRes$comparisons)
KTSP_STATs_Test_Agnostic[KTSP_STATs_Test_Agnostic == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic, KTSP_STATs_Test_Agnostic, file = "./Objs/KTSP_STATs_Agnostic.rda")
