
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

#save(myTSPs, file = './objs/metastasisPairs.rda')
colnames(myTSPs) <- c('gene1', 'gene2')
write_csv(as.data.frame(myTSPs), file = './objs/metastasis_pairs.csv')

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
  # get the pairs
  rownames(KTSP_Train_Agnostic$TSPs) <- gsub(',', '-', rownames(KTSP_Train_Agnostic$TSPs))
  rownames(KTSP_Train_Mech$TSPs) <- gsub(',', '-', rownames(KTSP_Train_Mech$TSPs))
  pairs_agnostic <- paste(rownames(KTSP_Train_Agnostic$TSPs), collapse = ',')
  pairs_Mech <- paste(rownames(KTSP_Train_Mech$TSPs), collapse = ',')
  # get the invidiual genes
  gene1_agnostic <- paste(KTSP_Train_Agnostic$TSPs[,1], collapse = ',')
  gene2_agnostic <- paste(KTSP_Train_Agnostic$TSPs[,2], collapse = ',')
  gene1_Mech <- paste(KTSP_Train_Mech$TSPs[,1], collapse = ',')
  gene2_Mech <- paste(KTSP_Train_Mech$TSPs[,2], collapse = ',')
  ##
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
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic, pairs_agnostic, pairs_Mech, gene1_agnostic, gene2_agnostic, gene1_Mech, gene2_Mech))
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

load("./Objs/KTSP/bootobjectKTSP_AdhesionActivationO2response.rda")


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

################
## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_50 <- data.frame(NofPairs = All_50[, "N_Pairs_Mech"])
Agnostic_NofPairs_50 <- data.frame(NofPairs = All_50[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_50$modelType <- "Mechanistic"
Agnostic_NofPairs_50$modelType <- "Agnostic"

ModelCompare_NofPairs_50 <- rbind(Mechanistic_NofPairs_50, Agnostic_NofPairs_50)

ModelCompare_NofPairs_50$NofFeatAgn <- "50 Genes"

# ###########################################################################
## Save for the main figure
ModelCompare_KTSP <- rbind(ModelCompareAUCTrain_50, ModelCompareAUCTest_50)
ModelCompare_KTSP$algorithm <- "KTSP"
save(ModelCompare_KTSP, file = "./Objs/KTSP/ModelCompare_KTSP.rda")

##############################################################
### Work with boot object 100  
All_100 <- bootobject_100$t
colnames(All_100) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic", "pairs_agnostic", "pairs_Mech", "gene1_agnostic", "gene2_agnostic", "gene1_Mech", "gene2_Mech")

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

################
## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_100 <- data.frame(NofPairs = All_100[, "N_Pairs_Mech"])
Agnostic_NofPairs_100 <- data.frame(NofPairs = All_100[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_100$modelType <- "Mechanistic"
Agnostic_NofPairs_100$modelType <- "Agnostic"

ModelCompare_NofPairs_100 <- rbind(Mechanistic_NofPairs_100, Agnostic_NofPairs_100)

ModelCompare_NofPairs_100$NofFeatAgn <- "100 Genes"

# ## Save for the main figure
#ModelCompare_KTSP <- rbind(ModelCompareAUCTrain_100, ModelCompareAUCTest_100)
#ModelCompare_KTSP$algorithm <- "KTSP"
#save(ModelCompare_KTSP, file = "./Objs/KTSP/ModelCompare_KTSP.rda")

###############
# get the common pairs (repeated)
pairs_agnostic <- strsplit(All_100[, 'pairs_agnostic'], ',')
pairs_mech <- strsplit(All_100[, 'pairs_Mech'], ',')

pairs_agnostic_df <- plyr::ldply(pairs_agnostic, rbind)
pairs_mech_df <- plyr::ldply(pairs_mech, rbind)

find_rep <- function(dat, feature){
  # NAs create problems in the function so we substitute that with "unknown"
  dat[is.na(dat)] <- "unknown"
  rep_rows <- sum(apply(dat, 1,  function(x) any(x == feature)))
  names(rep_rows) <- feature
  as.data.frame(rep_rows)
}

features_agnostic <- na.omit(unique(as.vector(as.matrix(pairs_agnostic_df))))
list_results_agnostic <- lapply(features_agnostic, find_rep, dat = pairs_agnostic_df)
sum_result_agnostic <- do.call(rbind, list_results_agnostic)
sum_result_agnostic$feature <- rownames(sum_result_agnostic)
sum_result_agnostic <- as.data.frame(sum_result_agnostic[order(sum_result_agnostic$rep_rows, decreasing = T), ])

features_mech <- na.omit(unique(as.vector(as.matrix(pairs_mech_df))))
list_results_mech <- lapply(features_mech, find_rep, dat = pairs_mech_df)
sum_result_mech <- do.call(rbind, list_results_mech)
sum_result_mech$feature <- rownames(sum_result_mech)
sum_result_mech <- as.data.frame(sum_result_mech[order(sum_result_mech$rep_rows, decreasing = T), ])

###
# bar plots
MechFreq <- ggplot(data=sum_result_mech[c(1:20), ], aes(x=rep_rows, y=reorder(feature, rep_rows))) +
  geom_col(width=0.5) + 
  #coord_cartesian(xlim = c(1, 600))
  scale_x_continuous(limits = c(1,550), n.breaks =10, oob = scales::squish) +
  labs(y = "Pair", x = "Frequency", title = "Mechanistic") +
  theme(plot.title = element_text(size=10, hjust=0.5), 
        axis.text.y = element_text(size=8), 
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8)
  )

AgnFreq <- ggplot(data=sum_result_agnostic[c(1:20), ], aes(x=rep_rows, y=reorder(feature, rep_rows))) +
  geom_col(width=0.5) + 
  #coord_cartesian(xlim = c(1, 600))
  scale_x_continuous(limits = c(1,350), n.breaks =10, oob = scales::squish) +
  labs(y = "Pair", x = "Frequency", title = "Agnostic") +
  theme(plot.title = element_text(size=10, hjust=0.5), 
        axis.text.y = element_text(size=8), 
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8)
  )


tiff(filename = "./Figs/prostate_frequency.tiff", width = 3000, height = 2000, res = 300)
((MechFreq + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10))) | 
    (AgnFreq + plot_layout(tag_level = "new") & theme(plot.tag = element_text(size = 10))) 
) +
  #plot_layout(widths = c(0.4, 1)) + 
  plot_annotation(
    title = 'The 20 most frequent gene pairs returned by the k-TSPs model',
    tag_levels = c('A'),
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )
dev.off()

#########
## get the individual genes for pairs of interest

# agnostic
agnostic_indvGns_good <- All_100 %>%
  as.data.frame() %>%
  dplyr::select(pairs_agnostic, gene1_agnostic, gene2_agnostic) %>%
  #mutate(tmp = strsplit(as.character(pairs_Mech),',')) %>%
  #unnest(tmp) %>%
  group_by(pairs_agnostic) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp) %>%
  #pivot_longer(cols = c(4:13), values_to = 'pair', values_drop_na = TRUE) %>%
  #select(-name) %>%
  #mutate(tmp_gene1 = strsplit(as.character(gene1_Mech),',')) %>%
  #group_by(pair) %>%
  #unnest(tmp_gene1, keep_empty = T) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp_gene1) %>%
  separate_rows(c(pairs_agnostic, gene1_agnostic, gene2_agnostic), sep = ',') %>%
  #select(-pairs_Mech) %>%
  dplyr::rename(gene1=gene1_agnostic, gene2=gene2_agnostic)

agnostic_genes_notGood <- sum_result_agnostic %>%
  dplyr::mutate(tmp = strsplit(as.character(feature),'-')) %>%
  dplyr::mutate(gene1 = map_chr(tmp, 1),
                gene2 = map_chr(tmp, 2)) %>%
  column_to_rownames(var = 'feature') %>%
  select(-tmp) %>%
  relocate(rep_rows, .after = gene2)

agnostic_genes_freq <- data.frame(rep_rows = agnostic_genes_notGood$rep_rows, 
                                  pairs_agnostic = rownames(agnostic_genes_notGood),
                                  row.names = rownames(agnostic_genes_notGood))


# merge
agnostic_indvGns_good_unique <- agnostic_indvGns_good[!duplicated(agnostic_indvGns_good$pairs_agnostic), ]
agnostic_indvGns_good_clean <- merge(x = agnostic_indvGns_good_unique, y = agnostic_genes_freq, 
                                     by = 'pairs_agnostic', suffixes = colnames(agnostic_indvGns_good_unique))

agnostic_indvGns_good_clean <- agnostic_indvGns_good_clean[order(agnostic_indvGns_good_clean$rep_rows, decreasing = T), ] 
agnostic_indvGns_good_clean$pairs_agnostic <- NULL

#agnostic_indvGns_good_clean <- filter(agnostic_indvGns_good_clean, !grepl("///", gene1, ignore.case = TRUE))
#agnostic_indvGns_good_clean <- filter(agnostic_indvGns_good_clean, !grepl("///", gene2, ignore.case = TRUE))

agnostic_indvGns_good_clean <- agnostic_indvGns_good_clean[agnostic_indvGns_good_clean$rep_rows>10, ]

#######
# mechanistic
mech_indvGns_good <- All_100 %>%
  as.data.frame() %>%
  dplyr::select(pairs_Mech, gene1_Mech, gene2_Mech) %>%
  #mutate(tmp = strsplit(as.character(pairs_Mech),',')) %>%
  #unnest(tmp) %>%
  group_by(pairs_Mech) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp) %>%
  #pivot_longer(cols = c(4:13), values_to = 'pair', values_drop_na = TRUE) %>%
  #select(-name) %>%
  #mutate(tmp_gene1 = strsplit(as.character(gene1_Mech),',')) %>%
  #group_by(pair) %>%
  #unnest(tmp_gene1, keep_empty = T) %>%
  #mutate(row = row_number()) %>%
  #spread(row, tmp_gene1) %>%
  separate_rows(c(pairs_Mech, gene1_Mech, gene2_Mech), sep = ',') %>%
  #select(-pairs_Mech) %>%
  dplyr::rename(gene1=gene1_Mech, gene2=gene2_Mech)


mech_genes_notGood <- sum_result_mech %>%
  dplyr::mutate(tmp = strsplit(as.character(feature),'-')) %>%
  dplyr::mutate(gene1 = map_chr(tmp, 1),
                gene2 = map_chr(tmp, 2)) %>%
  column_to_rownames(var = 'feature') %>%
  select(-tmp) %>%
  relocate(rep_rows, .after = gene2)

mech_genes_freq <- data.frame(rep_rows = mech_genes_notGood$rep_rows, 
                              pairs_Mech = rownames(mech_genes_notGood),
                              row.names = rownames(mech_genes_notGood))

# merge
mech_indvGns_good_unique <- mech_indvGns_good[!duplicated(mech_indvGns_good$pairs_Mech), ]
mech_indvGns_good_clean <- merge(x = mech_indvGns_good_unique, y = mech_genes_freq, 
                                 by = 'pairs_Mech', suffixes = colnames(mech_indvGns_good_unique))

mech_indvGns_good_clean <- mech_indvGns_good_clean[order(mech_indvGns_good_clean$rep_rows, decreasing = T), ] 
mech_indvGns_good_clean$pairs_Mech <- NULL

mech_indvGns_good_clean <- mech_indvGns_good_clean[mech_indvGns_good_clean$rep_rows>10, ]

#############################
## plot

# weighted bipartite network
library(igraph)

#mech_genes_net <- as.network(as.matrix(mech_indvGns), 
#                             ignore.eval = F, 
#                             names.eval = "rep_rows")

#ggnet2(mech_genes_net, label = T, label.size = 3)

##########################
#network_mech <- graph_from_data_frame(d=mech_indvGns, directed = T)
#network_mech <- graph_from_edgelist(as.matrix(mech_indvGns), directed = T) 
#deg_mech <- degree(network_mech, mode="all", normalized = T)

#rescale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}
#deg_mech2 <- rescale(degree(network_mech), 1, 1000, 1,10)

#l <- layout_nicely(network_mech, dim = 2)

# plot(network_mech, 
#      vertex.size=log(degree(network_mech)), 
#      vertex.label.dist=0.1,
#      edge.arrow.width=0.3, 
#      edge.arrow.size=0.2, 
#      edge.width = 0.5,
#      arrow.size = 0.5, 
#      arrow.width = 0.5,
#      #label.cex=0.05,
#      vertex.label.cex=0.6,
#      layout = layout_with_lgl)


network_mech <- graph_from_data_frame(as.matrix(mech_indvGns_good_clean), directed = T)
E(network_mech)$width <- log(mech_indvGns_good_clean$rep_rows)+1
#V(network_mech)$color <- ifelse(names(V(network_mech)) %in% mech_indvGns_good_clean$gene1, 'green', 'red')

set.seed(333)
tiff(filename = 'figs/prostate_mech.tiff', width =  3000, height = 2000, res = 200)
l_mech <- layout_nicely(network_mech)
l_mech <- norm_coords(l_mech, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5) #default -- scaled
plot(network_mech, 
     vertex.size=degree(network_mech), 
     vertex.label.dist=0,
     edge.arrow.width=0.1, 
     edge.arrow.size=0.1, 
     #edge.width = 0.5,
     arrow.size = 0.5, 
     arrow.width = 0.5,
     #label.cex=0.05,
     vertex.label.cex=0.4,
     rescale=F,
     layout = l_mech)
dev.off()

network_agnostic <- graph_from_data_frame(as.matrix(agnostic_indvGns_good_clean), directed = F) 
E(network_agnostic)$width <- log(agnostic_indvGns_good_clean$rep_rows)
tiff(filename = 'figs/prostate_agnostic.tiff', width =  3000, height = 2000, res = 200)
l_agnostic <- layout_nicely(network_agnostic)
l_agnostic <- norm_coords(l_agnostic, ymin=-1, ymax=1, xmin=-1.5, xmax=1.5) #default -- scaled
plot(network_agnostic, 
     vertex.size=degree(network_agnostic), 
     vertex.label.dist=0,
     edge.arrow.width=0.1, 
     edge.arrow.size=0.1, 
     #edge.width = 0.5,
     arrow.size = 0.5, 
     arrow.width = 0.5,
     #label.cex=0.05,
     vertex.label.cex=0.4,
     layout = l_agnostic
)
dev.off()


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

################
## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_200 <- data.frame(NofPairs = All_200[, "N_Pairs_Mech"])
Agnostic_NofPairs_200 <- data.frame(NofPairs = All_200[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_200$modelType <- "Mechanistic"
Agnostic_NofPairs_200$modelType <- "Agnostic"

ModelCompare_NofPairs_200 <- rbind(Mechanistic_NofPairs_200, Agnostic_NofPairs_200)

ModelCompare_NofPairs_200$NofFeatAgn <- "200 Genes"

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

################
## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Mech"])
Agnostic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_500$modelType <- "Mechanistic"
Agnostic_NofPairs_500$modelType <- "Agnostic"

ModelCompare_NofPairs_500 <- rbind(Mechanistic_NofPairs_500, Agnostic_NofPairs_500)

ModelCompare_NofPairs_500$NofFeatAgn <- "500 Genes"

####################################################################################
# combine the number of pairs distribution

ModelCompare_KTSP_Npairs <- rbind(ModelCompare_NofPairs_50,
                                  ModelCompare_NofPairs_100,
                                  ModelCompare_NofPairs_200,
                                  ModelCompare_NofPairs_500
)

save(ModelCompare_KTSP_Npairs, file = "./Objs/KTSP/ModelCompare_KTSP_Npairs.rda")


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