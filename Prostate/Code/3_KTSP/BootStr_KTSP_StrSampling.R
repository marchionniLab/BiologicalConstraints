
rm(list = ls())

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
library(boot)
library(patchwork)

### Load expression and phenotype data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

###############################
## Load the selected genes
Genes1 <- read.delim("./objs/GO_Adhesion.txt")
Genes1 <- as.matrix(Genes1)
Genes1 <- Genes1[-1,]

Genes2 <- read.delim("./objs/GO_Activation.txt")
Genes2 <- as.matrix(Genes2)
Genes2 <- Genes2[-1,]

Genes3 <- read.delim("./objs/GO_O2Response.txt")
Genes3 <- as.matrix(Genes3)
Genes3 <- Genes3[-1,]

Genes <- c(Genes1,Genes2, Genes3)
Genes <- Genes[!duplicated(Genes)]

myTSPs <- t(combn(Genes,2))

#################################
### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

#########
### TRAINING using restricted pairs
#######

### Set Feature number and max k
ktsp <- c(1:50) #12
featNo <- nrow(usedTrainMat)

# Combine the expression matrix with the phenotype in 1 data frame
Data <- cbind(t(usedTrainMat), usedTrainGroup)
Data <- as.data.frame(Data)
Data$usedTrainGroup <- as.factor(Data$usedTrainGroup)
levels(Data[, "usedTrainGroup"]) <- c("No_Mets", "Mets")

######################################################
## Mechanistic and agnostic using top 50 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=50)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}

set.seed(333)
bootobject_50 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15)

######################################################
## Mechanistic and agnostic using top 100 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=100)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_100 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15) 

######################################################
## Mechanistic and agnostic using top 200 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=200)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_200 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15) 

######################################################
## Mechanistic and agnostic using top 500 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=500)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_500 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 15) 

############################################################
## Save all
save(bootobject_50, bootobject_100, bootobject_200, bootobject_500, file = "./Objs/KTSP/bootobjectKTSP_AdhesionActivationO2response.rda")

#load("./Objs/KTSP/bootobjectKTSP_AdhesionActivationO2response.rda")


##############################################################
### Work with boot object 50  
All_50 <- bootobject_50$t
colnames(All_50) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_50 <- All_50[,"Diff_Agnostic"]
range(Diff_Agnostic_50)
quantile(Diff_Agnostic_50, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_50 <- All_50[,"Diff_Mechanistic"]
range(Diff_Mechanistic_50)
quantile(Diff_Mechanistic_50, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_50 <- data.frame(AUC = All_50[, "AUC_Train_Mech"])
AgnosticAUCTrain_50 <- data.frame(AUC = All_50[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_50$modelType <- "Mechanistic"
AgnosticAUCTrain_50$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_50 <- rbind(MechanisticAUCTrain_50, AgnosticAUCTrain_50)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_50 <- data.frame(AUC = All_50[, "AUC_Test_Mech"])
AgnosticAUCTest_50 <- data.frame(AUC = All_50[, "AUC_Test_Agnostic"])

MechanisticAUCTest_50$modelType <- "Mechanistic"
AgnosticAUCTest_50$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_50 <- rbind(MechanisticAUCTest_50, AgnosticAUCTest_50)


## Save the AUCs in the training and testing data
ModelCompareAUCTrain_50$data_type <- "Training"
ModelCompareAUCTest_50$data_type <- "Testing"

ModelCompareAUCTrain_50$NofFeatAgn <- "50_Genes"
ModelCompareAUCTest_50$NofFeatAgn <- "50_Genes"
# 
# ###########################################################################
## Save for the main figure
ModelCompare_KTSP <- rbind(ModelCompareAUCTrain_50, ModelCompareAUCTest_50)
ModelCompare_KTSP$algorithm <- "KTSP"
save(ModelCompare_KTSP, file = "./Objs/KTSP/ModelCompare_KTSP.rda")

##############################################################
### Work with boot object 100  
All_100 <- bootobject_100$t
colnames(All_100) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_100 <- All_100[,"Diff_Agnostic"]
range(Diff_Agnostic_100)
quantile(Diff_Agnostic_100, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_100 <- All_100[,"Diff_Mechanistic"]
range(Diff_Mechanistic_100)
quantile(Diff_Mechanistic_100, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_100 <- data.frame(AUC = All_100[, "AUC_Train_Mech"])
AgnosticAUCTrain_100 <- data.frame(AUC = All_100[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_100$modelType <- "Mechanistic"
AgnosticAUCTrain_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_100 <- rbind(MechanisticAUCTrain_100, AgnosticAUCTrain_100)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_100 <- data.frame(AUC = All_100[, "AUC_Test_Mech"])
AgnosticAUCTest_100 <- data.frame(AUC = All_100[, "AUC_Test_Agnostic"])

MechanisticAUCTest_100$modelType <- "Mechanistic"
AgnosticAUCTest_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_100 <- rbind(MechanisticAUCTest_100, AgnosticAUCTest_100)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_100$data_type <- "Training"
ModelCompareAUCTest_100$data_type <- "Testing"

ModelCompareAUCTrain_100$NofFeatAgn <- "100_Genes"
ModelCompareAUCTest_100$NofFeatAgn <- "100_Genes"

# ## Save for the main figure
#ModelCompare_KTSP <- rbind(ModelCompareAUCTrain_100, ModelCompareAUCTest_100)
#ModelCompare_KTSP$algorithm <- "KTSP"
#save(ModelCompare_KTSP, file = "./Objs/KTSP/ModelCompare_KTSP.rda")


###############################
##############################################################
### Work with boot object 200  
All_200 <- bootobject_200$t
colnames(All_200) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_200 <- All_200[,"Diff_Agnostic"]
range(Diff_Agnostic_200)
quantile(Diff_Agnostic_200, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_200 <- All_200[,"Diff_Mechanistic"]
range(Diff_Mechanistic_200)
quantile(Diff_Mechanistic_200, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_200 <- data.frame(AUC = All_200[, "AUC_Train_Mech"])
AgnosticAUCTrain_200 <- data.frame(AUC = All_200[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_200$modelType <- "Mechanistic"
AgnosticAUCTrain_200$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_200 <- rbind(MechanisticAUCTrain_200, AgnosticAUCTrain_200)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_200 <- data.frame(AUC = All_200[, "AUC_Test_Mech"])
AgnosticAUCTest_200 <- data.frame(AUC = All_200[, "AUC_Test_Agnostic"])

MechanisticAUCTest_200$modelType <- "Mechanistic"
AgnosticAUCTest_200$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_200 <- rbind(MechanisticAUCTest_200, AgnosticAUCTest_200)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_200$data_type <- "Training"
ModelCompareAUCTest_200$data_type <- "Testing"

ModelCompareAUCTrain_200$NofFeatAgn <- "200_Genes"
ModelCompareAUCTest_200$NofFeatAgn <- "200_Genes"

###############################
##############################################################
### Work with boot object 500  
All_500 <- bootobject_500$t
colnames(All_500) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_500 <- All_500[,"Diff_Agnostic"]
range(Diff_Agnostic_500)
quantile(Diff_Agnostic_500, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_500 <- All_500[,"Diff_Mechanistic"]
range(Diff_Mechanistic_500)
quantile(Diff_Mechanistic_500, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_500 <- data.frame(AUC = All_500[, "AUC_Train_Mech"])
AgnosticAUCTrain_500 <- data.frame(AUC = All_500[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_500$modelType <- "Mechanistic"
AgnosticAUCTrain_500$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_500 <- rbind(MechanisticAUCTrain_500, AgnosticAUCTrain_500)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_500 <- data.frame(AUC = All_500[, "AUC_Test_Mech"])
AgnosticAUCTest_500 <- data.frame(AUC = All_500[, "AUC_Test_Agnostic"])

MechanisticAUCTest_500$modelType <- "Mechanistic"
AgnosticAUCTest_500$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_500 <- rbind(MechanisticAUCTest_500, AgnosticAUCTest_500)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_500$data_type <- "Training"
ModelCompareAUCTest_500$data_type <- "Testing"

ModelCompareAUCTrain_500$NofFeatAgn <- "500_Genes"
ModelCompareAUCTest_500$NofFeatAgn <- "500_Genes"

####################################################################################
####################################################################################
####################################################################################
# Combine all together in one dataframe

ModelCompare_KTSP_DiffNoFeat <- rbind(ModelCompareAUCTrain_50,
                                      ModelCompareAUCTest_50,
                                      ModelCompareAUCTrain_100,
                                      ModelCompareAUCTest_100,
                                      ModelCompareAUCTrain_200,
                                      ModelCompareAUCTest_200,
                                      ModelCompareAUCTrain_500,
                                      ModelCompareAUCTest_500
)

save(ModelCompare_KTSP_DiffNoFeat, file = "./Objs/KTSP/ModelCompare_KTSP_DiffNoFeat.rda")

####################################################################################
####################################################################################