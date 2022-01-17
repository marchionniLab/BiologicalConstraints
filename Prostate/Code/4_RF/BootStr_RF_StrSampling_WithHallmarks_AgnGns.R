###################################################
## RF
# Mech

rm(list = ls())

## Load necessary packages
library(caret)
library(switchBox)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(randomForest)
library(patchwork)
library(boot)



############
## Mechanistic 


## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Hallmarks.rda")
load("./Objs/MetastasisDataGood.rda")


usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


predictor_data_Train_Mech <- t(KTSP_STATs_Train_Mechanistic)
predictor_data_Test_Mech <- t(KTSP_STATs_Test_Mechanistic)

DataMech_Train <- cbind(predictor_data_Train_Mech, usedTrainGroup)
DataMech_Train <- as.data.frame(DataMech_Train)
DataMech_Train$usedTrainGroup <- as.factor(DataMech_Train$usedTrainGroup)
levels(DataMech_Train[, "usedTrainGroup"]) <- c("No_Mets", "Mets")

names(DataMech_Train) <- make.names(names(DataMech_Train))

colnames(predictor_data_Train_Mech) <- make.names(colnames(predictor_data_Train_Mech))
colnames(predictor_data_Test_Mech) <- make.names(colnames(predictor_data_Test_Mech))

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Mech <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Mech <- Importance_Mech[order(Importance_Mech, decreasing = TRUE), ]
  Importance_Mech <- Importance_Mech[Importance_Mech > 0]
  N_ImportanVariables <- length(Importance_Mech)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Mech, type = "vote")
  ROCTrainMech <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestMech <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainMech$auc, ROCTestMech$auc, N_ImportanVariables))
}

set.seed(333)
bootobjectMech <- boot(data= DataMech_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

AUCs_RF_Mech <- bootobjectMech$t
colnames(AUCs_RF_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

###################################################
## RF
# Agnostic top 74 DEGs


## Load the data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")



### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 50)

## Subset the expression matrix to the top DE genes only
usedTrainMat <- usedTrainMat[TopDEgenes, ]
usedTestMat <- usedTestMat[TopDEgenes, ]
########################


predictor_data_Train_Agnostic <- t(usedTrainMat)
predictor_data_Test_Agnostic <- t(usedTestMat)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("No_Mets", "Mets")

#names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

#colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
#colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

## Finally we run the RF algorithm. 

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_50 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

################################################################################
# Agnostic top 100 DEGs

## Load the data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")



### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 100)

## Subset the expression matrix to the top DE genes only
usedTrainMat <- usedTrainMat[TopDEgenes, ]
usedTestMat <- usedTestMat[TopDEgenes, ]
########################


predictor_data_Train_Agnostic <- t(usedTrainMat)
predictor_data_Test_Agnostic <- t(usedTestMat)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("No_Mets", "Mets")

#names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

#colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
#colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

## Finally we run the RF algorithm. 

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_100 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

################################################################################
# Agnostic top 200 DEGs

## Load the data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 200)

## Subset the expression matrix to the top DE genes only
usedTrainMat <- usedTrainMat[TopDEgenes, ]
usedTestMat <- usedTestMat[TopDEgenes, ]
########################


predictor_data_Train_Agnostic <- t(usedTrainMat)
predictor_data_Test_Agnostic <- t(usedTestMat)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("No_Mets", "Mets")

#names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

#colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
#colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

## Finally we run the RF algorithm. 

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_200 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 


################################################################################
# Agnostic top 500 DEGs

## Load the data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 500)

## Subset the expression matrix to the top DE genes only
usedTrainMat <- usedTrainMat[TopDEgenes, ]
usedTestMat <- usedTestMat[TopDEgenes, ]
########################


predictor_data_Train_Agnostic <- t(usedTrainMat)
predictor_data_Test_Agnostic <- t(usedTestMat)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("No_Mets", "Mets")

#names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

#colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
#colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

## Finally we run the RF algorithm. 

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_500 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

################################################################################

## Save all bootobjects
save(bootobjectMech, bootobjectAgnostic_50, bootobjectAgnostic_100, bootobjectAgnostic_200, bootobjectAgnostic_500, file = "./Objs/RF/RFBootObjects_Hallmarks.rda")

## load
load("./Objs/RF/RFBootObjects_Hallmarks.rda")

##################################################################################
##################################################################################


## Work with Agnostic bootobject 74 vs mechanistic

AUCs_RF_Agnostic_50 <- bootobjectAgnostic_50$t
colnames(AUCs_RF_Agnostic_50) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


AUCs_RF_Mech <- bootobjectMech$t
colnames(AUCs_RF_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_50 <- AUCs_RF_Agnostic_50[, "AUC_Train"] - AUCs_RF_Agnostic_50[, "AUC_Test"]
quantile(DiffAgnostic_50, c(0.025, 0.975))

## Calculate the difference and CI of the difference between mechanistic training and mechanistic testing
DiffMech <- AUCs_RF_Mech[, "AUC_Train"] - AUCs_RF_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUC_Train <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Train"])
AgnosticAUC_Train_50 <- data.frame(AUC = AUCs_RF_Agnostic_50[, "AUC_Train"])

MechanisticAUC_Train$modelType <- "Mech"
AgnosticAUC_Train_50$modelType <- "Agnostic"

ModelCompareAUC_Train_50 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_50)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUC_Test <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Test"])
AgnosticAUC_Test_50 <- data.frame(AUC = AUCs_RF_Agnostic_50[, "AUC_Test"])

MechanisticAUC_Test$modelType <- "Mech"
AgnosticAUC_Test_50$modelType <- "Agnostic"

ModelCompareAUC_Test_50 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_50)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_50$data_type <- "Training"
ModelCompareAUC_Test_50$data_type <- "Testing"

ModelCompare_RF <- rbind(ModelCompareAUC_Train_50, ModelCompareAUC_Test_50)
ModelCompare_RF$algorithm <- "RF"
save(ModelCompare_RF, file = "./Objs/RF/ModelCompare_RF.rda")

##########
## The N of important features (output)
# ImportanceAgnostic <- data.frame(NofOutputFeatures = AUCs_RF_Agnostic[, "N_ImportanVariables"])
# ImportanceMech <- data.frame(NofOutputFeatures = AUCs_RF_Mech[, "N_ImportanVariables"])
# 
# ImportanceAgnostic$modelType <- "Agnostic"
# ImportanceMech$modelType <- "Mech"
# 
# ModelCompare_ImportantFeatures <- rbind(ImportanceAgnostic, ImportanceMech)

############
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=8)
)

DiffHist_Agnostic_50 <- ggplot(as.data.frame(DiffAgnostic_50), aes(DiffAgnostic_50, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0.10, 0.50)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0.10, 0.50)) +
  labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme


AUC_Train_DistrHist_50 <- ggplot(ModelCompareAUC_Train_50, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the training data") + My_Theme

AUC_Test_DistrHist_50 <- ggplot(ModelCompareAUC_Test_50, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the testing data") + My_Theme

# BarPlotImportance <- ggplot(data=ModelCompare_ImportantFeatures, aes(x=NofOutputFeatures, fill = modelType)) +
#   geom_bar(stat="bin") +
#   scale_x_continuous(limits = c(30, 80)) +
#   labs(title="Distribution of the number of important features (output)") + 
#   #facet_grid(~modelType, scale='free_x') +
#   My_Theme

png("./Figs/RF/RF_BS_AUC_50.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_50 / AUC_Test_DistrHist_50) | (DiffHist_Agnostic_50 / DiffHist_Mech) + plot_layout(widths = c(1, 1)) #| BarPlotImportance
dev.off()

###################################################################################3
###################################################################################
## Work with Agnostic bootobject 100 vs mechanistic

AUCs_RF_Agnostic_100 <- bootobjectAgnostic_100$t
colnames(AUCs_RF_Agnostic_100) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_100 <- AUCs_RF_Agnostic_100[, "AUC_Train"] - AUCs_RF_Agnostic_100[, "AUC_Test"]
quantile(DiffAgnostic_100, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_100 <- data.frame(AUC = AUCs_RF_Agnostic_100[, "AUC_Train"])

AgnosticAUC_Train_100$modelType <- "Agnostic"

ModelCompareAUC_Train_100 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_100)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_100 <- data.frame(AUC = AUCs_RF_Agnostic_100[, "AUC_Test"])

AgnosticAUC_Test_100$modelType <- "Agnostic"

ModelCompareAUC_Test_100 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_100)

############
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=8)
)

DiffHist_Agnostic_100 <- ggplot(as.data.frame(DiffAgnostic_100), aes(DiffAgnostic_100, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0.10, 0.50)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUC_Train_100, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the training data") + My_Theme

AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUC_Test_100, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the testing data") + My_Theme

png("./Figs/RF/RF_BS_AUC_100.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (DiffHist_Agnostic_100 / DiffHist_Mech) + plot_layout(widths = c(1, 1)) #| BarPlotImportance
dev.off()

###################################################################################3
###################################################################################
## Work with Agnostic bootobject 200 vs mechanistic

AUCs_RF_Agnostic_200 <- bootobjectAgnostic_200$t
colnames(AUCs_RF_Agnostic_200) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_200 <- AUCs_RF_Agnostic_200[, "AUC_Train"] - AUCs_RF_Agnostic_200[, "AUC_Test"]
quantile(DiffAgnostic_200, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_200 <- data.frame(AUC = AUCs_RF_Agnostic_200[, "AUC_Train"])

AgnosticAUC_Train_200$modelType <- "Agnostic"

ModelCompareAUC_Train_200 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_200)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_200 <- data.frame(AUC = AUCs_RF_Agnostic_200[, "AUC_Test"])

AgnosticAUC_Test_200$modelType <- "Agnostic"

ModelCompareAUC_Test_200 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_200)

############
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=8)
)

DiffHist_Agnostic_200 <- ggplot(as.data.frame(DiffAgnostic_200), aes(DiffAgnostic_200, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0.10, 0.50)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUC_Train_200, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the training data") + My_Theme

AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUC_Test_200, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the testing data") + My_Theme

png("./Figs/RF/RF_BS_AUC_200.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (DiffHist_Agnostic_200 / DiffHist_Mech) + plot_layout(widths = c(1, 1)) #| BarPlotImportance
dev.off()

###################################################################################3
###################################################################################
## Work with Agnostic bootobject 500 vs mechanistic

AUCs_RF_Agnostic_500 <- bootobjectAgnostic_500$t
colnames(AUCs_RF_Agnostic_500) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_500 <- AUCs_RF_Agnostic_500[, "AUC_Train"] - AUCs_RF_Agnostic_500[, "AUC_Test"]
quantile(DiffAgnostic_500, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_500 <- data.frame(AUC = AUCs_RF_Agnostic_500[, "AUC_Train"])

AgnosticAUC_Train_500$modelType <- "Agnostic"

ModelCompareAUC_Train_500 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_500)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_500 <- data.frame(AUC = AUCs_RF_Agnostic_500[, "AUC_Test"])

AgnosticAUC_Test_500$modelType <- "Agnostic"

ModelCompareAUC_Test_500 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_500)

############
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 5),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 5),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=8)
)

DiffHist_Agnostic_500 <- ggplot(as.data.frame(DiffAgnostic_500), aes(DiffAgnostic_500, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0.10, 0.50)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUC_Train_500, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the training data") + My_Theme

AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUC_Test_500, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the testing data") + My_Theme

png("./Figs/RF/RF_BS_AUC_500.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) | (DiffHist_Agnostic_500 / DiffHist_Mech) + plot_layout(widths = c(1, 1)) #| BarPlotImportance
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
AUC_Train_DistrHist_50 <- ggplot(ModelCompareAUC_Train_50, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 50 DEGs) vs mechanistic RF models in the training data") + My_Theme2

AUC_Test_DistrHist_50 <- ggplot(ModelCompareAUC_Test_50, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 50 DEGs) vs mechanistic RF models in the testing data") + My_Theme2
#############################
AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUC_Train_100, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 100 DEGs) vs mechanistic RF models in the training data") + My_Theme2

AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUC_Test_100, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 100 DEGs) vs mechanistic RF models in the testing data") + My_Theme2
##############################
AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUC_Train_200, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 200 DEGs) vs mechanistic RF models in the training data") + My_Theme2

AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUC_Test_200, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 200 DEGs) vs mechanistic RF models in the testing data") + My_Theme2
###############################
AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUC_Train_500, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 500 DEGs) vs mechanistic RF models in the training data") + My_Theme2

AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUC_Test_500, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.5, 1)) +
  labs(title="Agnostic (top 500 DEGs) vs mechanistic RF models in the testing data") + My_Theme2


# All in one figure
png("./Figs/RF/RF_BS_AUC_DiffFeatures.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist_50 / AUC_Test_DistrHist_50) | (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) + plot_layout(widths = c(1, 1)) 
dev.off()


