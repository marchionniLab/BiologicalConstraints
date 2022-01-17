###################################################
## SVM

rm(list = ls())

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(mltools)
library(parallel)
library(doParallel)

###################################################
## SVM
# Mechanistic
## Load data
load("./Objs/KTSP/TNBC_KTSP_STATs_Mechanistic_NotchAndMyc2_100.rda")
load("./Objs/ChemoDataNew.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Mechanistic)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Mechanistic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Mechanistic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Mechanistic) <- make.names(names(Data_train_Mechanistic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 3, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Mech <- varImp(SVM, scale = TRUE)
  Importance_Mech <- Importance_Mech$importance
  Importance_Mech <- Importance_Mech[order(Importance_Mech$Resistant, decreasing = TRUE),]
  Importance_Mech <- Importance_Mech[Importance_Mech$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Mech)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainMechanistic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestMechanistic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainMechanistic$auc, ROCTestMechanistic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectMech <- boot(data= Data_train_Mechanistic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

#########################################################################################
#########################################################################################

# Agnostic

# 25 pairs (from 50 DEGs)

## Load data
load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_25.rda")
load("./Objs/ChemoDataNew.rda")


####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_25)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical
Testing <- t(KTSP_STATs_Test_Agnostic_25)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Agnostic) <- make.names(names(Data_train_Agnostic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 5, scale = 0.01, C = 0.25)
# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_25 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15)

#########################################################################################
#########################################################################################

# Agnostic

# 50 pairs (from 100 DEGs)

## Load data
load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_50.rda")
load("./Objs/ChemoDataNew.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_50)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic_50)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Agnostic) <- make.names(names(Data_train_Agnostic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 1, scale = 0.01, C = 0.25)
# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_50 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

#########################################################################################
#########################################################################################

# Agnostic

# 100 pairs (from 200 DEGs)

## Load data
load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_100.rda")
load("./Objs/ChemoDataNew.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_100)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic_100)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Agnostic) <- make.names(names(Data_train_Agnostic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 1, scale = 0.01, C = 0.25)
# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_100 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

#########################################################################################
#########################################################################################

# Agnostic

# 250 pairs (from 500 DEGs)

## Load data
load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_250.rda")
load("./Objs/ChemoDataNew.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_250)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic_250)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Agnostic) <- make.names(names(Data_train_Agnostic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 1, scale = 0.01, C = 0.25)
# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Resistant, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Resistant > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_250 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
########################################################################################
## Save all Objects
save(bootobjectMech, bootobjectAgnostic_25, bootobjectAgnostic_50, bootobjectAgnostic_100, bootobjectAgnostic_250, file= "./Objs/SVM/SVMBootObjects_NotchAndMyc_AgnPairs_mech100pairs.rda")

## Load
#load("./Objs/SVM/SVMBootObjects_NotchAndMyc_AgnPairs.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 25 vs mechanistic

AUCs_SVM_Agnostic_25 <- bootobjectAgnostic_25$t
colnames(AUCs_SVM_Agnostic_25) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_25 <- AUCs_SVM_Agnostic_25[, "AUC_Train"] - AUCs_SVM_Agnostic_25[, "AUC_Test"]
quantile(DiffAgnostic_25, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain_25 <- data.frame(AUC = AUCs_SVM_Agnostic_25[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mechanistic"
AgnosticAUCTrain_25$modelType <- "Agnostic_Pairs"

ModelCompareAUCTrain_25 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_25)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_25 <- data.frame(AUC = AUCs_SVM_Agnostic_25[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_25$modelType <- "Agnostic_Pairs"

ModelCompareAUCTest_25 <- rbind(MechanisticAUCTest, AgnosticAUCTest_25)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_25$data_type <- "Training"
ModelCompareAUCTest_25$data_type <- "Testing"

ModelCompareAUCTrain_25$NofFeatAgn <- "25_Pairs"
ModelCompareAUCTest_25$NofFeatAgn <- "25_Pairs"

#############################################################################
## Save for the main figure
ModelCompare_SVM <- rbind(ModelCompareAUCTrain_25, ModelCompareAUCTest_25)
ModelCompare_SVM$algorithm <- "SVM"
save(ModelCompare_SVM, file = "./Objs/SVM/ModelCompare_SVM_AgnPairs_mech100pairs.rda")
 
# #########################################################################
# ############################################################
# # Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_50_mech100pairs.rda")

# Combine
ModelCompareAUC_Train_25_50 <- rbind(ModelCompareAUCTrain_25, ModelCompareAUCTrain_50)
ModelCompareAUC_Test_25_50 <- rbind(ModelCompareAUCTest_25, ModelCompareAUCTest_50)

########################################################################################
########################################################################################

## Work with Agnostic bootobject 50 vs mechanistic

AUCs_SVM_Agnostic_50 <- bootobjectAgnostic_50$t
colnames(AUCs_SVM_Agnostic_50) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_50 <- AUCs_SVM_Agnostic_50[, "AUC_Train"] - AUCs_SVM_Agnostic_50[, "AUC_Test"]
quantile(DiffAgnostic_50, c(0.025, 0.975))

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain_50 <- data.frame(AUC = AUCs_SVM_Agnostic_50[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mechanistic"
AgnosticAUCTrain_50$modelType <- "Agnostic_Pairs"

ModelCompareAUCTrain_50 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_50)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_50 <- data.frame(AUC = AUCs_SVM_Agnostic_50[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_50$modelType <- "Agnostic_Pairs"

ModelCompareAUCTest_50 <- rbind(MechanisticAUCTest, AgnosticAUCTest_50)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_50$data_type <- "Training"
ModelCompareAUCTest_50$data_type <- "Testing"

ModelCompareAUCTrain_50$NofFeatAgn <- "50_Pairs"
ModelCompareAUCTest_50$NofFeatAgn <- "50_Pairs"

#############################################################################
## Save for the main figure
#ModelCompare_SVM <- rbind(ModelCompareAUCTrain_50, ModelCompareAUCTest_50)
#ModelCompare_SVM$algorithm <- "SVM"
#save(ModelCompare_SVM, file = "./Objs/SVM/ModelCompare_SVM_AgnPairs.rda")

#########################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_100_mech100pairs.rda")

# Combine
ModelCompareAUC_Train_50_100 <- rbind(ModelCompareAUCTrain_50, ModelCompareAUCTrain_100)
ModelCompareAUC_Test_50_100 <- rbind(ModelCompareAUCTest_50, ModelCompareAUCTest_100)

########################################################################################
########################################################################################

## Work with Agnostic bootobject 100 vs mechanistic

AUCs_SVM_Agnostic_100 <- bootobjectAgnostic_100$t
colnames(AUCs_SVM_Agnostic_100) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_100 <- AUCs_SVM_Agnostic_100[, "AUC_Train"] - AUCs_SVM_Agnostic_100[, "AUC_Test"]
quantile(DiffAgnostic_100, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUCTrain_100 <- data.frame(AUC = AUCs_SVM_Agnostic_100[, "AUC_Train"])

AgnosticAUCTrain_100$modelType <- "Agnostic_Pairs"

ModelCompareAUCTrain_100 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_100)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_100 <- data.frame(AUC = AUCs_SVM_Agnostic_100[, "AUC_Test"])

AgnosticAUCTest_100$modelType <- "Agnostic_Pairs"

ModelCompareAUCTest_100 <- rbind(MechanisticAUCTest, AgnosticAUCTest_100)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_100$data_type <- "Training"
ModelCompareAUCTest_100$data_type <- "Testing"

ModelCompareAUCTrain_100$NofFeatAgn <- "100_Pairs"
ModelCompareAUCTest_100$NofFeatAgn <- "100_Pairs"

#########################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_200_mech100pairs.rda")

# Combine
ModelCompareAUC_Train_100_200 <- rbind(ModelCompareAUCTrain_100, ModelCompareAUCTrain_200)
ModelCompareAUC_Test_100_200 <- rbind(ModelCompareAUCTest_100, ModelCompareAUCTest_200)

########################################################################################
########################################################################################

## Work with Agnostic bootobject 500 vs mechanistic

AUCs_SVM_Agnostic_250 <- bootobjectAgnostic_250$t
colnames(AUCs_SVM_Agnostic_250) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_250 <- AUCs_SVM_Agnostic_250[, "AUC_Train"] - AUCs_SVM_Agnostic_250[, "AUC_Test"]
quantile(DiffAgnostic_250, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUCTrain_250 <- data.frame(AUC = AUCs_SVM_Agnostic_250[, "AUC_Train"])

AgnosticAUCTrain_250$modelType <- "Agnostic_Pairs"

ModelCompareAUCTrain_250 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_250)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_250 <- data.frame(AUC = AUCs_SVM_Agnostic_250[, "AUC_Test"])

AgnosticAUCTest_250$modelType <- "Agnostic_Pairs"

ModelCompareAUCTest_250 <- rbind(MechanisticAUCTest, AgnosticAUCTest_250)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_250$data_type <- "Training"
ModelCompareAUCTest_250$data_type <- "Testing"

ModelCompareAUCTrain_250$NofFeatAgn <- "250_Pairs"
ModelCompareAUCTest_250$NofFeatAgn <- "250_Pairs"

#########################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_500_mech100pairs.rda")

# Combine
ModelCompareAUC_Train_250_500 <- rbind(ModelCompareAUCTrain_250, ModelCompareAUCTrain_500)
ModelCompareAUC_Test_250_500 <- rbind(ModelCompareAUCTest_250, ModelCompareAUCTest_500)

####################################################################################
####################################################################################
####################################################################################
# Combine all together in one dataframe

ModelCompare_SVM_DiffNoFeat <- rbind(ModelCompareAUC_Train_25_50,
                                     ModelCompareAUC_Test_25_50,
                                     ModelCompareAUC_Train_50_100,
                                     ModelCompareAUC_Test_50_100,
                                     ModelCompareAUC_Train_100_200,
                                     ModelCompareAUC_Test_100_200,
                                     ModelCompareAUC_Train_250_500,
                                     ModelCompareAUC_Test_250_500
)

save(ModelCompare_SVM_DiffNoFeat, file = "./Objs/SVM/ModelCompare_SVM_DiffNoFeat_mech100pairs.rda")

####################################################################################
####################################################################################
####################################################################################

