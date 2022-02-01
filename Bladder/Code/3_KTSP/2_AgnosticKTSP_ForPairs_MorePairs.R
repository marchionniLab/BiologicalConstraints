################################################################################
### Mohamed Omar
### 10/4/2019
### Goal: Creating unrestricted K-TSP classifier
###############################################################################

rm(list = ls()) 

#setwd("/Volumes/Macintosh/Research/Projects/Bladder")

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
load("./Objs/progressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")

##########################################################
## Quantile normalize the datasets
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#####
### TRAINING using all expressed genes

## Set Feature number and max K
# featN <- nrow(usedTrainMat) # the same as in the mechanistic classifier 
#ktsp <- c(3:25)  # the same as in the mechanistic classifier
# ALL AVAILABLE PAIRS = 37 (as in the mechanistic)

###################################################################
### Train a classifier using the top 74 features > 37 pairs
# 
set.seed(333)

ktspPredictorUnRes_37 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = 37, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 74)
ktspPredictorUnRes_37

#######
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_37 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_37, CombineFunc = sum)
summary(ktspStatsTrainUnRes_37$statistics)


KTSP_STATs_Train_Agnostic_37 <- t(ktspStatsTrainUnRes_37$comparisons)
KTSP_STATs_Train_Agnostic_37[KTSP_STATs_Train_Agnostic_37 == FALSE] <- 0


##########
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_37 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_37, CombineFunc = sum)
summary(ktspStatsTestUnRes_37$statistics)


KTSP_STATs_Test_Agnostic_37 <- t(ktspStatsTestUnRes_37$comparisons)
KTSP_STATs_Test_Agnostic_37[KTSP_STATs_Test_Agnostic_37 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_37, KTSP_STATs_Test_Agnostic_37, file = "./Objs/KTSP/KTSP_STATs_Agnostic_37.rda")


#############################################################################
##############################################
###########
### Train a classifier using the top 100 features > 50 pairs
# 
set.seed(333)

ktspPredictorUnRes_50 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                         krange = 50, 
                                         FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 100)
ktspPredictorUnRes_50

##########
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_50 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_50, CombineFunc = sum)
summary(ktspStatsTrainUnRes_50$statistics)


KTSP_STATs_Train_Agnostic_50 <- t(ktspStatsTrainUnRes_50$comparisons)
KTSP_STATs_Train_Agnostic_50[KTSP_STATs_Train_Agnostic_50 == FALSE] <- 0


##########
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_50 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_50, CombineFunc = sum)
summary(ktspStatsTestUnRes_50$statistics)


KTSP_STATs_Test_Agnostic_50 <- t(ktspStatsTestUnRes_50$comparisons)
KTSP_STATs_Test_Agnostic_50[KTSP_STATs_Test_Agnostic_50 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_50, KTSP_STATs_Test_Agnostic_50, file = "./Objs/KTSP/KTSP_STATs_Agnostic_50.rda")

#############################################################################
##############################################
###########
### Train a classifier using the top 200 features > 100 pairs

set.seed(333)

ktspPredictorUnRes_100 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                          krange = 100, 
                                          FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 200)
ktspPredictorUnRes_100

############
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_100 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_100, CombineFunc = sum)
summary(ktspStatsTrainUnRes_100$statistics)


KTSP_STATs_Train_Agnostic_100 <- t(ktspStatsTrainUnRes_100$comparisons)
KTSP_STATs_Train_Agnostic_100[KTSP_STATs_Train_Agnostic_100 == FALSE] <- 0


##########
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_100 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_100, CombineFunc = sum)
summary(ktspStatsTestUnRes_100$statistics)


KTSP_STATs_Test_Agnostic_100 <- t(ktspStatsTestUnRes_100$comparisons)
KTSP_STATs_Test_Agnostic_100[KTSP_STATs_Test_Agnostic_100 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_100, KTSP_STATs_Test_Agnostic_100, file = "./Objs/KTSP/KTSP_STATs_Agnostic_100.rda")

#############################################################################
##############################################
###########
### Train a classifier using the top 500 features > 250 pairs

set.seed(333)

ktspPredictorUnRes_250 <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                          krange = 250, 
                                          FilterFunc = SWAP.Filter.Wilcoxon, featureNo= 500)
ktspPredictorUnRes_250

##############
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes_250 <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes_250, CombineFunc = sum)
summary(ktspStatsTrainUnRes_250$statistics)


KTSP_STATs_Train_Agnostic_250 <- t(ktspStatsTrainUnRes_250$comparisons)
KTSP_STATs_Train_Agnostic_250[KTSP_STATs_Train_Agnostic_250 == FALSE] <- 0


##########
#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes_250 <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes_250, CombineFunc = sum)
summary(ktspStatsTestUnRes_250$statistics)


KTSP_STATs_Test_Agnostic_250 <- t(ktspStatsTestUnRes_250$comparisons)
KTSP_STATs_Test_Agnostic_250[KTSP_STATs_Test_Agnostic_250 == FALSE] <- 0

save(KTSP_STATs_Train_Agnostic_250, KTSP_STATs_Test_Agnostic_250, file = "./Objs/KTSP/KTSP_STATs_Agnostic_250.rda")

