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
load("./Objs/KTSP/KTSP_STATs_Mechanistic.rda")
load("./Objs/ProgressionDataGood2.rda")


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

Grid_mech <- expand.grid(degree = 3, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_mech, metric = "ROC")
  Importance_Mech <- varImp(SVM, scale = TRUE)
  Importance_Mech <- Importance_Mech$importance
  Importance_Mech <- Importance_Mech[order(Importance_Mech$Progression, decreasing = TRUE),]
  Importance_Mech <- Importance_Mech[Importance_Mech$Progression > 0, ]
  N_ImportanVariables <- nrow(Importance_Mech)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainMechanistic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestMechanistic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainMechanistic$auc, ROCTestMechanistic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectMech <- boot(data= Data_train_Mechanistic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

save(bootobjectMech, file= "./Objs/SVM/SVM_MechBootObject.rda")

###################################################################################
### Agnostic

### top 74 DEGs

## Load data
load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")


### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 74)

## Subset the expression matrix to the top DE genes only
usedTrainMat_tmp <- usedTrainMat[TopDEgenes, ]
#usedTestMat_tmp <- usedTestMat[TopDEgenes, ]

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
Training_tmp <- t(usedTrainMat_tmp)

#names_train <- c(as.vector(rownames(usedTrainMat)))
#colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
#names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Training_tmp <- as.data.frame(Training_tmp)

Data_train_Agnostic <- cbind(Training, usedTrainGroup)
Data_train_Agnostic_tmp <- cbind(Training_tmp, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
#names_Test <- c(as.vector(rownames(usedTestMat)))
#colnames(Testing) <- names_Test
#names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

#########################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_agnostic74 <- train(usedTrainGroup~., data=Data_train_Agnostic_tmp, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly_agnostic74

###########################
# Use the best parameters in the bootstrap
Grid_agn74 <- expand.grid(degree = 3, scale = 10, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  # get the top 74 DEGs
  Top74genes <- SWAP.Filter.Wilcoxon(phenoGroup = d[,"usedTrainGroup"], inputMat = as.matrix(t(d[,!colnames(d) == "usedTrainGroup"])), featureNo = 74)
  
  # subset the data to those top genes
  d <- d[, c(Top74genes, 'usedTrainGroup')]
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_agn74, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Progression, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Progression > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_74 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
### Agnostic

### top 100 DEGs

load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")
 
 
### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

#names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 100)

## Subset the expression matrix to the top DE genes only
usedTrainMat_tmp <- usedTrainMat[TopDEgenes, ]
usedTestMat_tmp <- usedTestMat[TopDEgenes, ]

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
#Training <- t(usedTrainMat) # we already have this from above

Training_tmp <- t(usedTrainMat_tmp)


## Making sure that sample names are identical in both Training and usedTrainGroup
#all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
#Training <- as.data.frame(Training)
Training_tmp <- as.data.frame(Training_tmp)

#Data_train_Agnostic <- cbind(Training, usedTrainGroup)
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
Grid_agn100 <- expand.grid(degree = 3, scale = 0.1, C = 0.25)

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
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Progression, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Progression > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_100 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
### Agnostic

### top 200 DEGs

load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")


### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = 200)

## Subset the expression matrix to the top DE genes only
usedTrainMat_tmp <- usedTrainMat[TopDEgenes, ]

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training_tmp <- t(usedTrainMat)

## Combining the expression matrix and the phenotype in one data frame
Training_tmp <- as.data.frame(Training_tmp)
Data_train_Agnostic_tmp <- cbind(Training_tmp, usedTrainGroup)

#########################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_agnostic200 <- train(usedTrainGroup~., data=Data_train_Agnostic_tmp, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly_agnostic200

###########################
# Use the best parameters in the bootstrap
Grid_agn200 <- expand.grid(degree = 3, scale = 0.1, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  
  # get the top 200 DEGs
  Top200genes <- SWAP.Filter.Wilcoxon(phenoGroup = d[,"usedTrainGroup"], inputMat = as.matrix(t(d[,!colnames(d) == "usedTrainGroup"])), featureNo = 200)
  
  # subset the data to those top genes
  d <- d[, c(Top200genes, 'usedTrainGroup')]
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_agn200, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Progression, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Progression > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}

set.seed(333)
bootobjectAgnostic_200 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
### Agnostic

### top 500 DEGs

load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

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

#########################
## Get the best parameters

control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

# 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly_agnostic500 <- train(usedTrainGroup~., data=Data_train_Agnostic_tmp, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
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
  
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid_agn500, metric = "ROC")
  Importance_Agnostic <- varImp(SVM, scale = TRUE)
  Importance_Agnostic <- Importance_Agnostic$importance
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Progression, decreasing = TRUE),]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Progression > 0, ]
  N_ImportanVariables <- nrow(Importance_Agnostic)
  Training <- d
  PhenoTrain <- d$usedTrainGroup
  Training$usedTrainGroup <- NULL
  train_preds <- predict(SVM, newdata = Training, type = "prob")
  test_preds <- predict(SVM, newdata = Testing, type = "prob")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("NoProgression", "Progression"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}

set.seed(333)
bootobjectAgnostic_500 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
########################################################################################
## Save all Objects
save(bootobjectMech, bootobjectAgnostic_74, bootobjectAgnostic_100, bootobjectAgnostic_200, bootobjectAgnostic_500, file= "./Objs/SVM/SVMBootObjects_new.rda")

## Load
load("./Objs/SVM/SVMBootObjects_new.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 74 vs mechanistic

AUCs_SVM_Agnostic_74 <- bootobjectAgnostic_74$t
colnames(AUCs_SVM_Agnostic_74) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_74 <- AUCs_SVM_Agnostic_74[, "AUC_Train"] - AUCs_SVM_Agnostic_74[, "AUC_Test"]
quantile(DiffAgnostic_74, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain_74 <- data.frame(AUC = AUCs_SVM_Agnostic_74[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mechanistic"
AgnosticAUCTrain_74$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_74 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_74)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_74 <- data.frame(AUC = AUCs_SVM_Agnostic_74[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_74$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_74 <- rbind(MechanisticAUCTest, AgnosticAUCTest_74)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_74$data_type <- "Training"
ModelCompareAUCTest_74$data_type <- "Testing"

ModelCompareAUCTrain_74$NofFeatAgn <- "74_Genes"
ModelCompareAUCTest_74$NofFeatAgn <- "74_Genes"

save(ModelCompareAUCTrain_74, ModelCompareAUCTest_74, file = "./Objs/SVM/ModelCompareAUC_74_new.rda")

############################################################################
# Save for the main figure
ModelCompare_SVM <- rbind(ModelCompareAUCTrain_74, ModelCompareAUCTest_74)
ModelCompare_SVM$algorithm <- "SVM"
save(ModelCompare_SVM, file = "./Objs/SVM/ModelCompare_SVM_new.rda")

############################################################################

######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# DiffHist_Agnostic_74 <- ggplot(as.data.frame(DiffAgnostic_74), aes(DiffAgnostic_74, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_74 <- ggplot(ModelCompareAUCTrain_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_74 <- ggplot(ModelCompareAUCTest_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the testing data") + My_Theme
# 
# # BarPlotImportance <- ggplot(data=ModelCompare_ImportantFeatures, aes(x=NofOutputFeatures, fill = modelType)) +
# #   geom_bar(stat="bin") +
# #   scale_x_continuous(limits = c(30, 80)) +
# #   labs(title="Distribution of the number of important features (output)") + 
# #   facet_grid(~modelType, scale='free_x') +
# #   My_Theme
# 
# png("./Figs/SVM/SVM_BS_AUC_74.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_74 / AUC_Test_DistrHist_74) | (DiffHist_Agnostic_74 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

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

AgnosticAUCTrain_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTrain_100 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_100)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_100 <- data.frame(AUC = AUCs_SVM_Agnostic_100[, "AUC_Test"])

AgnosticAUCTest_100$modelType <- "Agnostic_DEGs"

ModelCompareAUCTest_100 <- rbind(MechanisticAUCTest, AgnosticAUCTest_100)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_100$data_type <- "Training"
ModelCompareAUCTest_100$data_type <- "Testing"

ModelCompareAUCTrain_100$NofFeatAgn <- "100_Genes"
ModelCompareAUCTest_100$NofFeatAgn <- "100_Genes"

save(ModelCompareAUCTrain_100, ModelCompareAUCTest_100, file = "./Objs/SVM/ModelCompareAUC_100_new.rda")

######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# DiffHist_Agnostic_100 <- ggplot(as.data.frame(DiffAgnostic_100), aes(DiffAgnostic_100, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUCTrain_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUCTest_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the testing data") + My_Theme
# 
# # BarPlotImportance <- ggplot(data=ModelCompare_ImportantFeatures, aes(x=NofOutputFeatures, fill = modelType)) +
# #   geom_bar(stat="bin") +
# #   scale_x_continuous(limits = c(30, 80)) +
# #   labs(title="Distribution of the number of important features (output)") + 
# #   facet_grid(~modelType, scale='free_x') +
# #   My_Theme
# 
# png("./Figs/SVM/SVM_BS_AUC_100.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (DiffHist_Agnostic_100 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

########################################################################################
########################################################################################

## Work with Agnostic bootobject 200 vs mechanistic

AUCs_SVM_Agnostic_200 <- bootobjectAgnostic_200$t
colnames(AUCs_SVM_Agnostic_200) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_200 <- AUCs_SVM_Agnostic_200[, "AUC_Train"] - AUCs_SVM_Agnostic_200[, "AUC_Test"]
quantile(DiffAgnostic_200, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


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

######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# DiffHist_Agnostic_200 <- ggplot(as.data.frame(DiffAgnostic_200), aes(DiffAgnostic_200, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUCTrain_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUCTest_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the testing data") + My_Theme
# 
# 
# png("./Figs/SVM/SVM_BS_AUC_200.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (DiffHist_Agnostic_200 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

########################################################################################
########################################################################################

## Work with Agnostic bootobject 500 vs mechanistic

AUCs_SVM_Agnostic_500 <- bootobjectAgnostic_500$t
colnames(AUCs_SVM_Agnostic_500) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_500 <- AUCs_SVM_Agnostic_500[, "AUC_Train"] - AUCs_SVM_Agnostic_500[, "AUC_Test"]
quantile(DiffAgnostic_500, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


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

######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# DiffHist_Agnostic_500 <- ggplot(as.data.frame(DiffAgnostic_500), aes(DiffAgnostic_500, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUCTrain_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUCTest_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the testing data") + My_Theme
# 
# # BarPlotImportance <- ggplot(data=ModelCompare_ImportantFeatures, aes(x=NofOutputFeatures, fill = modelType)) +
# #   geom_bar(stat="bin") +
# #   scale_x_continuous(limits = c(30, 80)) +
# #   labs(title="Distribution of the number of important features (output)") + 
# #   facet_grid(~modelType, scale='free_x') +
# #   My_Theme
# 
# png("./Figs/SVM/SVM_BS_AUC_500.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) | (DiffHist_Agnostic_500 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()
# 
# ###################################################################################
# ###########################################
# ###########################################
# My_Theme2 = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=9)
# )
# 
# #############################
# AUC_Train_DistrHist_74 <- ggplot(ModelCompareAUCTrain_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 74 DEGs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_74 <- ggplot(ModelCompareAUCTest_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 74 DEGs) vs mechanistic SVM models in the testing data") + My_Theme2
# #############################
# AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUCTrain_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 100 DEGs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUCTest_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 100 DEGs) vs mechanistic SVM models in the testing data") + My_Theme2
# ##############################
# AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUCTrain_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 200 DEGs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUCTest_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 200 DEGs) vs mechanistic SVM models in the testing data") + My_Theme2
# ###############################
# AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUCTrain_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 500 DEGs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUCTest_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (top 500 DEGs) vs mechanistic SVM models in the testing data") + My_Theme2
# 
# 
# # All in one figure
# png("./Figs/BootStrap_Diff_No_features/SVM_BS_AUC_DiffFeatures.png", width = 3000, height = 1500, res = 150)
# (AUC_Train_DistrHist_74 / AUC_Test_DistrHist_74) | (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) + plot_layout(widths = c(1, 1)) 
# dev.off()
# 
# 
# 
# 
# 
