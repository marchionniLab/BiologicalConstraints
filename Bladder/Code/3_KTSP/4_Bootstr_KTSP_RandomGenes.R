
## Run it with the agnostic then the mechanistic
# The desired outcome is to repeat the training process n times and calculate the AUC in the testing data for each n repeat
# Take the difference between the AUC in the test data using mechanistic and agnostic approach
# AUC_Test_Mech - AUC_Test_Agnostic
# Then calculate the Confidence interval of this difference


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
l1 = load("../../Objs/progressionDataGood2.rda")
l2 = load("../../Objs/Correlation/RGenes.rda")

#########
## Load the selected genes *(TF-MiR from Lotte)
TF_MiR <- load("../../../../Genes/allTSPs.rda")

#TF_MiR <- load("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/USER/Mohamed/MechanisticModels/Genes/TF_MiR.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

if(FALSE){

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
levels(Data[, "usedTrainGroup"]) <- c("NoProgression", "Progression")


######################################################
## Mechanistic and agnostic using 74 random genes

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  
  sel.genes = sample(rownames(ExprMat), 74)
  
  # Finally the function
  
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat[sel.genes, ], phenoGroup = Pheno, krange=ktsp, FilterFunc = NULL)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_74 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R=1000, parallel = "multicore", ncpus = 12) 

######################################################
## Mechanistic and agnostic using top 100 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  
  sel.genes = sample(rownames(ExprMat), 100)
  
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat[sel.genes, ], phenoGroup = Pheno, krange=ktsp, FilterFunc = NULL)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_100 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 12) 

######################################################
## Mechanistic and agnostic using top 200 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  
  sel.genes = sample(rownames(ExprMat), 200)
  
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat[sel.genes, ], phenoGroup = Pheno, krange=ktsp, FilterFunc = NULL)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_200 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R= 1000, parallel = "multicore", ncpus = 12) 

######################################################
## Mechanistic and agnostic using top 500 DEGs

SWAP.Train.KTSPStrap <- function(data, indices) {
  d = data[indices, ]
  Pheno <- d$usedTrainGroup
  d$usedTrainGroup <- NULL
  ExprMat = t(data.matrix(d, rownames.force = T)) # transpose the matrix so that rows are genes and columns are samples
  
  sel.genes = sample(rownames(ExprMat), 500)
  
  # Finally the function
  KTSP_Train_Agnostic <- SWAP.Train.KTSP(inputMat=ExprMat[sel.genes, ], phenoGroup = Pheno, krange=ktsp, FilterFunc = NULL)
  KTSP_Train_Mech <- SWAP.Train.KTSP(inputMat=ExprMat, phenoGroup = Pheno, krange=ktsp, FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)  
  N_Pairs_Agnostic <- nrow(KTSP_Train_Agnostic$TSPs)
  N_Pairs_Mech <- nrow(KTSP_Train_Mech$TSPs)
  ktspStatsTrainAgnostic <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTestAgnostic <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Agnostic, CombineFunc = sum)
  ktspStatsTrainMech <- SWAP.KTSP.Statistics(inputMat = ExprMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ktspStatsTestMech <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = KTSP_Train_Mech, CombineFunc = sum)
  ROCTrainAgnostic <- roc(Pheno, ktspStatsTrainAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, ktspStatsTestAgnostic$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTrainMech <- roc(Pheno, ktspStatsTrainMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  ROCTestMech <- roc(usedTestGroup, ktspStatsTestMech$statistics, plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", levels = c("NoProgression", "Progression"), direction = "<", col="blue", lwd=2, grid=TRUE)
  AUC_Train_Agnostic <- ROCTrainAgnostic$auc
  AUC_Test_Agnostic <- ROCTestAgnostic$auc
  AUC_Train_Mech <- ROCTrainMech$auc
  AUC_Test_Mech <- ROCTestMech$auc
  Diff_Agnostic <- AUC_Train_Agnostic - AUC_Test_Agnostic
  Diff_Mechanistic <- AUC_Train_Mech - AUC_Test_Mech
  return(c(N_Pairs_Agnostic, AUC_Train_Agnostic, AUC_Test_Agnostic, N_Pairs_Mech, AUC_Train_Mech, AUC_Test_Mech, Diff_Agnostic, Diff_Mechanistic))
}


set.seed(333)
bootobject_500 <- boot(data= Data, statistic= SWAP.Train.KTSPStrap, R=1000, parallel = "multicore", ncpus=12) 

############################################################
## Save all
save(bootobject_74, bootobject_100, bootobject_200, bootobject_500, file = "../../Objs/KTSP/bootobjectKTSP_RAND.rda")


} # if-false

load("../../Objs/KTSP/bootobjectKTSP_RAND.rda")



##############################################################
### Work with boot object 50  
All_74 <- bootobject_74$t
colnames(All_74) <- c("N_Pairs_Agnostic", "AUC_Train_Agnostic", "AUC_Test_Agnostic", "N_Pairs_Mech", "AUC_Train_Mech", "AUC_Test_Mech", "Diff_Agnostic", "Diff_Mechanistic")

## Calculate the difference and CI of the difference (Training data)
Diff_Agnostic_74 <- All_74[,"Diff_Agnostic"]
range(Diff_Agnostic_74)
quantile(Diff_Agnostic_74, c(0.025, 0.975))

## Calculate the difference and CI of the difference (Testing data)
Diff_Mechanistic_74 <- All_74[,"Diff_Mechanistic"]
range(Diff_Mechanistic_74)
quantile(Diff_Mechanistic_74, c(0.025, 0.975))

## Calculate the difference and CI of the difference
#Diff <- MechKTSP_AUC_Test - AgnosticKTSP_AUC_Test
#quantile(Diff, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain_74 <- data.frame(AUC = All_74[, "AUC_Train_Mech"])
AgnosticAUCTrain_74 <- data.frame(AUC = All_74[, "AUC_Train_Agnostic"])

MechanisticAUCTrain_74$modelType <- "Mechanistic"
AgnosticAUCTrain_74$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_74 <- rbind(MechanisticAUCTrain_74, AgnosticAUCTrain_74)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest_74 <- data.frame(AUC = All_74[, "AUC_Test_Mech"])
AgnosticAUCTest_74 <- data.frame(AUC = All_74[, "AUC_Test_Agnostic"])

MechanisticAUCTest_74$modelType <- "Mechanistic"
AgnosticAUCTest_74$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_74 <- rbind(MechanisticAUCTest_74, AgnosticAUCTest_74)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_74$data_type <- "Training"
ModelCompareAUCTest_74$data_type <- "Testing"

ModelCompareAUCTrain_74$NofFeatAgn <- "74_Genes"
ModelCompareAUCTest_74$NofFeatAgn <- "74_Genes"

###########################################################################
## Save for the main figure
ModelCompare_KTSP <- rbind(ModelCompareAUCTrain_74, ModelCompareAUCTest_74)
ModelCompare_KTSP$algorithm <- "KTSP"
save(ModelCompare_KTSP, file = "../../Objs/KTSP/ModelCompare_KTSP_RAND.rda")
###########################################################################





#### ---------------------- plotting -------------------------------


if(FALSE){


######################
## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_74 <- data.frame(NofPairs = All_74[, "N_Pairs_Mech"])
Agnostic_NofPairs_74 <- data.frame(NofPairs = All_74[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_74$modelType <- "Mech"
Agnostic_NofPairs_74$modelType <- "Agnostic"

ModelCompare_NofPairs_74 <- rbind(Mechanistic_NofPairs_74, Agnostic_NofPairs_74)

##################################
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=7)
)

DiffHist_Agnostic_74 <- ggplot(as.data.frame(Diff_Agnostic_74), aes(Diff_Agnostic_74, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

DiffHist_Mech_74 <- ggplot(as.data.frame(Diff_Mechanistic_74), aes(Diff_Mechanistic_74, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme

AUC_Train_DistrHist_74 <- ggplot(ModelCompareAUCTrain_74, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic (74 random genes) and mechanistic KTSP models in the training data") + My_Theme

AUC_Test_DistrHist_74 <- ggplot(ModelCompareAUCTest_74, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic KTSP models in the testing data") + My_Theme

NofPairs_DistrHist_74 <- ggplot(ModelCompare_NofPairs_74, aes(NofPairs, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(1, 50)) +
  labs(title="N of pairs distribution of the agnostic (74 random genes) and mechanistic KTSP models") + My_Theme

png("../../Figs/KTSP/KTSP_BS_RAND_AUC_74.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_74 / AUC_Test_DistrHist_74) | (DiffHist_Agnostic_74 / DiffHist_Mech_74 | NofPairs_DistrHist_74) + plot_layout(widths = c(1, 1.5))
dev.off()

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


## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_100 <- data.frame(NofPairs = All_100[, "N_Pairs_Mech"])
Agnostic_NofPairs_100 <- data.frame(NofPairs = All_100[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_100$modelType <- "Mech"
Agnostic_NofPairs_100$modelType <- "Agnostic"

ModelCompare_NofPairs_100 <- rbind(Mechanistic_NofPairs_100, Agnostic_NofPairs_100)

##################################
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=7)
)

DiffHist_Agnostic_100 <- ggplot(as.data.frame(Diff_Agnostic_100), aes(Diff_Agnostic_100, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

DiffHist_Mech_100 <- ggplot(as.data.frame(Diff_Mechanistic_100), aes(Diff_Mechanistic_100, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme

AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUCTrain_100, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic (100 random genes) and mechanistic KTSP models in the training data") + My_Theme

AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUCTest_100, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic (100 random genes) and mechanistic KTSP models in the testing data") + My_Theme

NofPairs_DistrHist_100 <- ggplot(ModelCompare_NofPairs_100, aes(NofPairs, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(1, 50)) +
  labs(title="N of pairs distribution of the agnostic and mechanistic KTSP models") + My_Theme

png("../../Figs/KTSP/KTSP_BS_RAND_AUC_100.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (DiffHist_Agnostic_100 / DiffHist_Mech_100 | NofPairs_DistrHist_100) + plot_layout(widths = c(1, 1.5))
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

## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_200 <- data.frame(NofPairs = All_200[, "N_Pairs_Mech"])
Agnostic_NofPairs_200 <- data.frame(NofPairs = All_200[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_200$modelType <- "Mech"
Agnostic_NofPairs_200$modelType <- "Agnostic"

ModelCompare_NofPairs_200 <- rbind(Mechanistic_NofPairs_200, Agnostic_NofPairs_200)

##################################
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=7)
)

My_Theme2 = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=15)
)

DiffHist_Agnostic_200 <- ggplot(as.data.frame(Diff_Agnostic_200), aes(Diff_Agnostic_200, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

DiffHist_Mech_200 <- ggplot(as.data.frame(Diff_Mechanistic_200), aes(Diff_Mechanistic_200, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme

AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUCTrain_200, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic (200 random genes) and mechanistic KTSP models in the training data") + My_Theme

AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUCTest_200, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic (200 random genes) and mechanistic KTSP models in the testing data") + My_Theme

NofPairs_DistrHist_200 <- ggplot(ModelCompare_NofPairs_200, aes(NofPairs, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(1, 50)) +
  labs(title="N of pairs distribution of the agnostic and mechanistic KTSP models") + My_Theme

png("../../Figs/KTSP/KTSP_BS_RAND_AUC_200.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (DiffHist_Agnostic_200 / DiffHist_Mech_200 | NofPairs_DistrHist_200) + plot_layout(widths = c(1, 1.5))
dev.off()

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

ModelCompare_KTSP_DiffNoFeat <- rbind(ModelCompareAUCTrain_74,
                                    ModelCompareAUCTest_74,
                                    ModelCompareAUCTrain_100,
                                    ModelCompareAUCTest_100,
                                    ModelCompareAUCTrain_200,
                                    ModelCompareAUCTest_200,
                                    ModelCompareAUCTrain_500,
                                    ModelCompareAUCTest_500
)

save(ModelCompare_KTSP_DiffNoFeat, file = "../../Objs/KTSP/ModelCompare_KTSP_RAND_DiffNoFeat.rda")

####################################################################################
####################################################################################
####################################################################################


## Plot the distributions of the N of Pairs from both methods
Mechanistic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Mech"])
Agnostic_NofPairs_500 <- data.frame(NofPairs = All_500[, "N_Pairs_Agnostic"])

Mechanistic_NofPairs_500$modelType <- "Mech"
Agnostic_NofPairs_500$modelType <- "Agnostic"

ModelCompare_NofPairs_500 <- rbind(Mechanistic_NofPairs_500, Agnostic_NofPairs_500)

##################################
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=7)
)

DiffHist_Agnostic_500 <- ggplot(as.data.frame(Diff_Agnostic_500), aes(Diff_Agnostic_500, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

DiffHist_Mech_500 <- ggplot(as.data.frame(Diff_Mechanistic_500), aes(Diff_Mechanistic_500, fill = "red")) +
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0, 0.4)) +
  labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme

AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUCTrain_500, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic (500 random genes) and mechanistic KTSP models in the training data") + My_Theme

AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUCTest_500, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic (500 random genes) and mechanistic KTSP models in the testing data") + My_Theme

NofPairs_DistrHist_500 <- ggplot(ModelCompare_NofPairs_500, aes(NofPairs, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(1, 50)) +
  labs(title="N of pairs distribution of the agnostic and mechanistic KTSP models") + My_Theme

png("../../Figs/KTSP/KTSP_BS_RAND_AUC_500.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) | (DiffHist_Agnostic_500 / DiffHist_Mech_500 | NofPairs_DistrHist_500) + plot_layout(widths = c(1, 1.5))
dev.off()


###########################################
###########################################
My_Theme2 = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=9)
)

#############################
AUC_Train_DistrHist_74 <- ggplot(ModelCompareAUCTrain_74, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (74 random genes) vs mechanistic KTSP models in the training data") + My_Theme2

AUC_Test_DistrHist_74 <- ggplot(ModelCompareAUCTest_74, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (74 random genes) vs mechanistic KTSP models in the testing data") + My_Theme2
#############################
AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUCTrain_100, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (100 random genes) vs mechanistic KTSP models in the training data") + My_Theme2

AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUCTest_100, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (100 random genes) vs mechanistic KTSP models in the testing data") + My_Theme2
##############################
AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUCTrain_200, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (200 random genes) vs mechanistic KTSP models in the training data") + My_Theme2

AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUCTest_200, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (200 random genes) vs mechanistic KTSP models in the testing data") + My_Theme2
###############################
AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUCTrain_500, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (500 random genes) vs mechanistic KTSP models in the training data") + My_Theme2

AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUCTest_500, aes(AUC, fill = modelType)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (500 random genes) vs mechanistic KTSP models in the testing data") + My_Theme2

# All in one figure
png("../../Figs/BootStrap_Diff_No_features/KTSP_BS_RAND_AUC_DiffFeatures.png", width = 3000, height = 1500, res = 150)
(AUC_Train_DistrHist_74 / AUC_Test_DistrHist_74) | (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) + plot_layout(widths = c(1, 1))
dev.off()



}






