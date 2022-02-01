
rm(list = ls())


####################### 
##Load required packages
library(switchBox)
library(Biobase)
library(limma)
library(pROC)
library(caret)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(plotROC)
library(xtable)
library(mltools)

#####################################################################
### Load data
load("./Objs/LOO/MetastasisData_GSE116918Out.rda")
load("./Objs/Correlation/RGenes.rda")

##########################################################
## Quantile normalize the datasets
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")[RGenes, ]
usedTestMat <- testMat[RGenes, ]

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

###########################################################################
### TRAINING using all expressed genes

## Set Feature number and max K
featN <- nrow(usedTrainMat) # the same as in the mechanistic classifier 
#ktsp <- c(3:25)  # the same as in the mechanistic classifier


###########
### Train a classifier using the default filter function
ktspPredictorUnRes <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = 50, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= featN)
ktspPredictorUnRes

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTrainUnRes$statistics)


KTSP_STATs_Train_Agnostic <- t(ktspStatsTrainUnRes$comparisons)
KTSP_STATs_Train_Agnostic[KTSP_STATs_Train_Agnostic == FALSE] <- 0

#################################################################################
###############################################################################
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTestUnRes$statistics)

KTSP_STATs_Test_Agnostic <- t(ktspStatsTestUnRes$comparisons)
KTSP_STATs_Test_Agnostic[KTSP_STATs_Test_Agnostic == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic, KTSP_STATs_Test_Agnostic, file = "./Objs/KTSP/LOO/KTSP_STATs_Agnostic_GSE116918Out.rda")

