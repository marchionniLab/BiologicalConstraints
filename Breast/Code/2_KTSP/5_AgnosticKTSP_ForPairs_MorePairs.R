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
library(DiagrammeR)
library(dplyr)
library(Ckmeans.1d.dp)

## Load data
load("./Objs/ChemoDataNew.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

featN <- nrow(usedTrainMat) # the same as in the mechanistic classifier 
#ktsp <- c(3:25)  # the same as in the mechanistic classifier


#############################################################################
### Train a classifier using the top 50 features > 25 pairs
ktspPredictorUnRes_25 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup,
                                      krange = 25,
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 50)
ktspPredictorUnRes_25


#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_25 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_25, CombineFunc = sum)
summary(ktspStatsTrainUnRes_25$statistics)


KTSP_STATs_Train_Agnostic_25 <- t(ktspStatsTrainUnRes_25$comparisons)
KTSP_STATs_Train_Agnostic_25[KTSP_STATs_Train_Agnostic_25 == FALSE] <- 0


###############################################################################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_25 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_25, CombineFunc = sum)
summary(ktspStatsTestUnRes_25$statistics)

KTSP_STATs_Test_Agnostic_25 <- t(ktspStatsTestUnRes_25$comparisons)
KTSP_STATs_Test_Agnostic_25[KTSP_STATs_Test_Agnostic_25 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_25, KTSP_STATs_Test_Agnostic_25, file = "./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_25.rda")

#############################################################################
###############################################################################

### Train a classifier using the top 168 features > 84 pairs
ktspPredictorUnRes_50 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                         krange = 50, 
                                         FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 100)
ktspPredictorUnRes_50


#####################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_50 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_50, CombineFunc = sum)
summary(ktspStatsTrainUnRes_50$statistics)


KTSP_STATs_Train_Agnostic_50 <- t(ktspStatsTrainUnRes_50$comparisons)
KTSP_STATs_Train_Agnostic_50[KTSP_STATs_Train_Agnostic_50 == FALSE] <- 0


######################################
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_50 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_50, CombineFunc = sum)
summary(ktspStatsTestUnRes_50$statistics)

KTSP_STATs_Test_Agnostic_50 <- t(ktspStatsTestUnRes_50$comparisons)
KTSP_STATs_Test_Agnostic_50[KTSP_STATs_Test_Agnostic_50 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_50, KTSP_STATs_Test_Agnostic_50, file = "./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_50.rda")

#############################################################################
###############################################################################

### Train a classifier using the top 200 features > 100 pairs
ktspPredictorUnRes_100 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                         krange = 100, 
                                         FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 200)
ktspPredictorUnRes_100


#############################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_100 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_100, CombineFunc = sum)
summary(ktspStatsTrainUnRes_100$statistics)


KTSP_STATs_Train_Agnostic_100 <- t(ktspStatsTrainUnRes_100$comparisons)
KTSP_STATs_Train_Agnostic_100[KTSP_STATs_Train_Agnostic_100 == FALSE] <- 0


################################
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_100 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_100, CombineFunc = sum)
summary(ktspStatsTestUnRes_100$statistics)

KTSP_STATs_Test_Agnostic_100 <- t(ktspStatsTestUnRes_100$comparisons)
KTSP_STATs_Test_Agnostic_100[KTSP_STATs_Test_Agnostic_100 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_100, KTSP_STATs_Test_Agnostic_100, file = "./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_100.rda")


#############################################################################
###############################################################################

### Train a classifier using the top 500 features > 250 pairs
ktspPredictorUnRes_250 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                          krange = 250, 
                                          FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 500)
ktspPredictorUnRes_250


###################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_250 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_250, CombineFunc = sum)
summary(ktspStatsTrainUnRes_250$statistics)


KTSP_STATs_Train_Agnostic_250 <- t(ktspStatsTrainUnRes_250$comparisons)
KTSP_STATs_Train_Agnostic_250[KTSP_STATs_Train_Agnostic_250 == FALSE] <- 0


####################################
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_250 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_250, CombineFunc = sum)
summary(ktspStatsTestUnRes_250$statistics)

KTSP_STATs_Test_Agnostic_250 <- t(ktspStatsTestUnRes_250$comparisons)
KTSP_STATs_Test_Agnostic_250[KTSP_STATs_Test_Agnostic_250 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_250, KTSP_STATs_Test_Agnostic_250, file = "./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_250.rda")



