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
library(boot)
library(patchwork)

##################################################
## SVM
# Mechanistic
## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")

####
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

names(Data_train_Mechanistic) <- make.names(names(Data_train_Mechanistic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

###########################################################################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_mech <- train(usedTrainGroup~., data=Data_train_Mechanistic, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly_mech

################################################
# Use the best parameters in the bootstrap

Grid_mech <- expand.grid(degree = 1, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_mech, metric = "ROC")
  Importance_Mech <- varImp(SVM, scale = TRUE)
  Importance_Mech <- Importance_Mech$importance
  Importance_Mech <- Importance_Mech[order(Importance_Mech$Mets, decreasing = TRUE),]
  Importance_Mech <- Importance_Mech[Importance_Mech$Mets > 0, ]
  N_ImportanVariables <- nrow(Importance_Mech)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainMechanistic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestMechanistic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainMechanistic$auc, ROCTestMechanistic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectMech <- boot(data= Data_train_Mechanistic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

save(bootobjectMech, file= "./Objs/SVM/SVM_MechBootObject_AdhesionActivationO2response.rda")

###################################################################################
### Agnostic

### top 50 DEGs

## Load data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")


### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 50)

## Subset the expression matrix to the top DE genes only
usedTrainMat_tmp <- usedTrainMat[TopDEgenes, ]

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
Training_tmp <- t(usedTrainMat_tmp)

## Making sure that sample names are identical in both Training and usedTrainGroup
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Training_tmp <- as.data.frame(Training_tmp)

Data_train_Agnostic <- cbind(Training, usedTrainGroup)
Data_train_Agnostic_tmp <- cbind(Training_tmp, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical
Testing <- t(usedTestMat)
all(rownames(Testing) == names(usedTestGroup))


#########################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_agnostic50 <- train(usedTrainGroup~., data=Data_train_Agnostic_tmp, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly_agnostic50

###########################
# Use the best parameters in the bootstrap
Grid_agn50 <- expand.grid(degree = 3, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  # get the top 50 DEGs
  Top50genes <- SWAP.Filter.Wilcoxon(phenoGroup = d[,"usedTrainGroup"], inputMat = as.matrix(t(d[,!colnames(d) == "usedTrainGroup"])), featureNo = 50)
  
  # subset the data to those top genes
  d <- d[, c(Top50genes, 'usedTrainGroup')]
  
  # train the model
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_agn50, metric = "ROC")
  
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Mets, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Mets > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_50 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15)

########################################################################################
### Agnostic

### top 100 DEGs

## Load data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")


### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]


####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


all(names(usedTrainGroup) == colnames(usedTrainMat))
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 100)

## Subset the expression matrix to the top DE genes only
usedTrainMat_tmp <- usedTrainMat[TopDEgenes, ]

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training_tmp <- t(usedTrainMat_tmp)

## Combining the expression matrix and the phenotype in one data frame
Training_tmp <- as.data.frame(Training_tmp)
Data_train_Agnostic_tmp <- cbind(Training_tmp, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
all(rownames(Testing) == names(usedTestGroup))

#########################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_agnostic100 <- train(usedTrainGroup~., data=Data_train_Agnostic_tmp, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly_agnostic100

###########################
# Use the best parameters in the bootstrap
Grid_agn100 <- expand.grid(degree = 3, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  # get the top 100 DEGs
  Top100genes <- SWAP.Filter.Wilcoxon(phenoGroup = d[,"usedTrainGroup"], inputMat = as.matrix(t(d[,!colnames(d) == "usedTrainGroup"])), featureNo = 100)
  
  # subset the data to those top genes
  d <- d[, c(Top100genes, 'usedTrainGroup')]
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_agn100, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Mets, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Mets > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_100 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
### Agnostic

### top 200 DEGs

## Load data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")


### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]


####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


all(names(usedTrainGroup) == colnames(usedTrainMat))
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 200)

## Subset the expression matrix to the top DE genes only
usedTrainMat_tmp <- usedTrainMat[TopDEgenes, ]

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat_tmp)

## Combining the expression matrix and the phenotype in one data frame
Training_tmp <- as.data.frame(Training_tmp)
Data_train_Agnostic_tmp <- cbind(Training_tmp, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
all(rownames(Testing) == names(usedTestGroup))

#########################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_agnostic200 <- train(usedTrainGroup~., data=Data_train_Agnostic, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly_agnostic200

###########################
# Use the best parameters in the bootstrap
Grid_agn200 <- expand.grid(degree = 2, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  # get the top 200 DEGs
  Top200genes <- SWAP.Filter.Wilcoxon(phenoGroup = d[,"usedTrainGroup"], inputMat = as.matrix(t(d[,!colnames(d) == "usedTrainGroup"])), featureNo = 200)
  
  # subset the data to those top genes
  d <- d[, c(Top200genes, 'usedTrainGroup')]
  
  # train the model
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_agn200, metric = "ROC")
  
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Mets, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Mets > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_200 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
### Agnostic

### top 500 DEGs

## Load data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")


### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]


####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


all(names(usedTrainGroup) == colnames(usedTrainMat))
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 500)

## Subset the expression matrix to the top DE genes only
usedTrainMat_tmp <- usedTrainMat[TopDEgenes, ]

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training_tmp <- t(usedTrainMat_tmp)

## Combining the expression matrix and the phenotype in one data frame
Training_tmp <- as.data.frame(Training_tmp)
Data_train_Agnostic_tmp <- cbind(Training_tmp, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
all(rownames(Testing) == names(usedTestGroup))

#########################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_agnostic500 <- train(usedTrainGroup~., data=Data_train_Agnostic, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly_agnostic500

###########################
# Use the best parameters in the bootstrap
Grid_agn500 <- expand.grid(degree = 3, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  # get the top 500 DEGs
  Top500genes <- SWAP.Filter.Wilcoxon(phenoGroup = d[,"usedTrainGroup"], inputMat = as.matrix(t(d[,!colnames(d) == "usedTrainGroup"])), featureNo = 500)
  
  # subset the data to those top genes
  d <- d[, c(Top500genes, 'usedTrainGroup')]
  
  # train the model
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_agn500, metric = "ROC")
  
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Mets, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Mets > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_500 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
########################################################################################
## Save all Objects
save(bootobjectMech, bootobjectAgnostic_50, bootobjectAgnostic_100, bootobjectAgnostic_200, bootobjectAgnostic_500, file= "./Objs/SVM/SVMBootObjects_AdhesionActivationO2response_new.rda")

## Load
load("./Objs/SVM/SVMBootObjects_AdhesionActivationO2response_new.rda")

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
AgnosticAUCTrain_50$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_50 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_50)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_50 <- data.frame(AUC = AUCs_SVM_Agnostic_50[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_50$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_50 <- rbind(MechanisticAUCTest, AgnosticAUCTest_50)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_50$data_type <- "Training"
ModelCompareAUCTest_50$data_type <- "Testing"

ModelCompareAUCTrain_50$NofFeatAgn <- "50_Genes"
ModelCompareAUCTest_50$NofFeatAgn <- "50_Genes"

save(ModelCompareAUCTrain_50, ModelCompareAUCTest_50, file = "./Objs/SVM/ModelCompareAUC_50_new.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 100 vs mechanistic

AUCs_SVM_Agnostic_100 <- bootobjectAgnostic_100$t
colnames(AUCs_SVM_Agnostic_100) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_100 <- AUCs_SVM_Agnostic_100[, "AUC_Train"] - AUCs_SVM_Agnostic_100[, "AUC_Test"]
quantile(DiffAgnostic_100, c(0.025, 0.975))

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain_100 <- data.frame(AUC = AUCs_SVM_Agnostic_100[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mechanistic"
AgnosticAUCTrain_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_100 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_100)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_100 <- data.frame(AUC = AUCs_SVM_Agnostic_100[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_100 <- rbind(MechanisticAUCTest, AgnosticAUCTest_100)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_100$data_type <- "Training"
ModelCompareAUCTest_100$data_type <- "Testing"

ModelCompareAUCTrain_100$NofFeatAgn <- "100_Genes"
ModelCompareAUCTest_100$NofFeatAgn <- "100_Genes"

save(ModelCompareAUCTrain_100, ModelCompareAUCTest_100, file = "./Objs/SVM/ModelCompareAUC_100_new.rda")

###########################################################################3
## Save for the main figure
ModelCompare_SVM <- rbind(ModelCompareAUCTrain_100, ModelCompareAUCTest_100)
ModelCompare_SVM$algorithm <- "SVM"
save(ModelCompare_SVM, file = "./Objs/SVM/ModelCompare_SVM_new.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 200 vs mechanistic

AUCs_SVM_Agnostic_200 <- bootobjectAgnostic_200$t
colnames(AUCs_SVM_Agnostic_200) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_200 <- AUCs_SVM_Agnostic_200[, "AUC_Train"] - AUCs_SVM_Agnostic_200[, "AUC_Test"]
quantile(DiffAgnostic_200, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUCTrain_200 <- data.frame(AUC = AUCs_SVM_Agnostic_200[, "AUC_Train"])

AgnosticAUCTrain_200$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_200 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_200)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_200 <- data.frame(AUC = AUCs_SVM_Agnostic_200[, "AUC_Test"])

AgnosticAUCTest_200$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_200 <- rbind(MechanisticAUCTest, AgnosticAUCTest_200)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_200$data_type <- "Training"
ModelCompareAUCTest_200$data_type <- "Testing"

ModelCompareAUCTrain_200$NofFeatAgn <- "200_Genes"
ModelCompareAUCTest_200$NofFeatAgn <- "200_Genes"

save(ModelCompareAUCTrain_200, ModelCompareAUCTest_200, file = "./Objs/SVM/ModelCompareAUC_200_new.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 500 vs mechanistic

AUCs_SVM_Agnostic_500 <- bootobjectAgnostic_500$t
colnames(AUCs_SVM_Agnostic_500) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_500 <- AUCs_SVM_Agnostic_500[, "AUC_Train"] - AUCs_SVM_Agnostic_500[, "AUC_Test"]
quantile(DiffAgnostic_500, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUCTrain_500 <- data.frame(AUC = AUCs_SVM_Agnostic_500[, "AUC_Train"])

AgnosticAUCTrain_500$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_500 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_500)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_500 <- data.frame(AUC = AUCs_SVM_Agnostic_500[, "AUC_Test"])

AgnosticAUCTest_500$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_500 <- rbind(MechanisticAUCTest, AgnosticAUCTest_500)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_500$data_type <- "Training"
ModelCompareAUCTest_500$data_type <- "Testing"

ModelCompareAUCTrain_500$NofFeatAgn <- "500_Genes"
ModelCompareAUCTest_500$NofFeatAgn <- "500_Genes"

save(ModelCompareAUCTrain_500, ModelCompareAUCTest_500, file = "./Objs/SVM/ModelCompareAUC_500_new.rda")

