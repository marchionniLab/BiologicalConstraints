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

###################################################

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


###############################################################################

# Agnostic

# 37 pairs (from 74 features)

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_37.rda")
load("./Objs/ProgressionDataGood2.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_37)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic_37)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Agnostic) <- make.names(names(Data_train_Agnostic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 2, scale = 0.001, C = 1)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
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
bootobjectAgnostic_37 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

###############################################################################

# Agnostic

# 50 pairs (from 100 features)

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_50.rda")
load("./Objs/ProgressionDataGood2.rda")


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

Grid <- expand.grid(degree = 2, scale = 0.001, C = 1)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
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
bootobjectAgnostic_50 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

###############################################################################

# Agnostic

# 100 pairs (from 200 features)

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_100.rda")
load("./Objs/ProgressionDataGood2.rda")


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

Grid <- expand.grid(degree = 2, scale = 0.001, C = 1)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
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


###############################################################################

# Agnostic

# 250 pairs (from 500 features)

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_250.rda")
load("./Objs/ProgressionDataGood2.rda")


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

Grid <- expand.grid(degree = 2, scale = 0.001, C = 1)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
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
bootobjectAgnostic_250 <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

########################################################################################
########################################################################################
## Save all Objects
save(bootobjectMech, bootobjectAgnostic_37, bootobjectAgnostic_50, bootobjectAgnostic_100, bootobjectAgnostic_250, file= "./Objs/SVM/SVMBootObjects_AgnPairs.rda")

## Load
load("./Objs/SVM/SVMBootObjects_AgnPairs.rda")

########################################################################################
########################################################################################

## Work with Agnostic bootobject 37 vs mechanistic

AUCs_SVM_Agnostic_37 <- bootobjectAgnostic_37$t
colnames(AUCs_SVM_Agnostic_37) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_37 <- AUCs_SVM_Agnostic_37[, "AUC_Train"] - AUCs_SVM_Agnostic_37[, "AUC_Test"]
quantile(DiffAgnostic_37, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain_37 <- data.frame(AUC = AUCs_SVM_Agnostic_37[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mechanistic"
AgnosticAUCTrain_37$modelType <- "Agnostic_Pairs"

ModelCompareAUCTrain_37 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_37)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest_37 <- data.frame(AUC = AUCs_SVM_Agnostic_37[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mechanistic"
AgnosticAUCTest_37$modelType <- "Agnostic_Pairs"

ModelCompareAUCTest_37 <- rbind(MechanisticAUCTest, AgnosticAUCTest_37)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_37$data_type <- "Training"
ModelCompareAUCTest_37$data_type <- "Testing"

ModelCompareAUCTrain_37$NofFeatAgn <- "37_Pairs"
ModelCompareAUCTest_37$NofFeatAgn <- "37_Pairs"


#########################################################################
## Save for the main figure
ModelCompare_SVM <- rbind(ModelCompareAUCTrain_37, ModelCompareAUCTest_37)
ModelCompare_SVM$algorithm <- "SVM"
save(ModelCompare_SVM, file = "./Objs/SVM/ModelCompare_SVM_AgnPairs.rda")
#########################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_74.rda")

# Combine
ModelCompareAUC_Train_37_74 <- rbind(ModelCompareAUCTrain_37, ModelCompareAUCTrain_74)
ModelCompareAUC_Test_37_74 <- rbind(ModelCompareAUCTest_37, ModelCompareAUCTest_74)


##############################################
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# DiffHist_Agnostic_37 <- ggplot(as.data.frame(DiffAgnostic_37), aes(DiffAgnostic_37, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_37 <- ggplot(ModelCompareAUCTrain_37, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_37 <- ggplot(ModelCompareAUCTest_37, aes(AUC, fill = modelType)) + 
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
# png("./Figs/SVM/SVM_BS_AUC_37_AgnPairs.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_37 / AUC_Test_DistrHist_37) | (DiffHist_Agnostic_37 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

########################################################################################
########################################################################################

## Work with Agnostic bootobject 100 vs mechanistic

AUCs_SVM_Agnostic_50 <- bootobjectAgnostic_50$t
colnames(AUCs_SVM_Agnostic_50) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_50 <- AUCs_SVM_Agnostic_50[, "AUC_Train"] - AUCs_SVM_Agnostic_50[, "AUC_Test"]
quantile(DiffAgnostic_50, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUCTrain_50 <- data.frame(AUC = AUCs_SVM_Agnostic_50[, "AUC_Train"])

AgnosticAUCTrain_50$modelType <- "Agnostic_Pairs"

ModelCompareAUCTrain_50 <- rbind(MechanisticAUCTrain, AgnosticAUCTrain_50)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUCTest_50 <- data.frame(AUC = AUCs_SVM_Agnostic_50[, "AUC_Test"])

AgnosticAUCTest_50$modelType <- "Agnostic_Pairs"

ModelCompareAUCTest_50 <- rbind(MechanisticAUCTest, AgnosticAUCTest_50)

## Save the AUCs in the training and testing data
ModelCompareAUCTrain_50$data_type <- "Training"
ModelCompareAUCTest_50$data_type <- "Testing"

ModelCompareAUCTrain_50$NofFeatAgn <- "50_Pairs"
ModelCompareAUCTest_50$NofFeatAgn <- "50_Pairs"

############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_100.rda")

# Combine
ModelCompareAUC_Train_50_100 <- rbind(ModelCompareAUCTrain_50, ModelCompareAUCTrain_100)
ModelCompareAUC_Test_50_100 <- rbind(ModelCompareAUCTest_50, ModelCompareAUCTest_100)

##############################################
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
# DiffHist_Agnostic_50 <- ggplot(as.data.frame(DiffAgnostic_50), aes(DiffAgnostic_50, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_50 <- ggplot(ModelCompareAUCTrain_50, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_50 <- ggplot(ModelCompareAUCTest_50, aes(AUC, fill = modelType)) + 
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
# png("./Figs/SVM/SVM_BS_AUC_50_AgnPairs.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_50 / AUC_Test_DistrHist_50) | (DiffHist_Agnostic_50 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

########################################################################################
########################################################################################

## Work with Agnostic bootobject 200 vs mechanistic

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

############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_200.rda")

# Combine
ModelCompareAUC_Train_100_200 <- rbind(ModelCompareAUCTrain_100, ModelCompareAUCTrain_200)
ModelCompareAUC_Test_100_200 <- rbind(ModelCompareAUCTest_100, ModelCompareAUCTest_200)

##############################################
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
# 
# png("./Figs/SVM/SVM_BS_AUC_100_AgnPairs.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (DiffHist_Agnostic_100 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

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

############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/SVM/ModelCompareAUC_500.rda")

# Combine
ModelCompareAUC_Train_250_500 <- rbind(ModelCompareAUCTrain_250, ModelCompareAUCTrain_500)
ModelCompareAUC_Test_250_500 <- rbind(ModelCompareAUCTest_250, ModelCompareAUCTest_500)

##############################################

####################################################################################
####################################################################################
####################################################################################
# Combine all together in one dataframe

ModelCompare_SVM_DiffNoFeat <- rbind(ModelCompareAUC_Train_37_74,
                                     ModelCompareAUC_Test_37_74,
                                     ModelCompareAUC_Train_50_100,
                                     ModelCompareAUC_Test_50_100,
                                     ModelCompareAUC_Train_100_200,
                                     ModelCompareAUC_Test_100_200,
                                     ModelCompareAUC_Train_250_500,
                                     ModelCompareAUC_Test_250_500
)

save(ModelCompare_SVM_DiffNoFeat, file = "./Objs/SVM/ModelCompare_SVM_DiffNoFeat.rda")

####################################################################################
####################################################################################
####################################################################################

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
# DiffHist_Agnostic_250 <- ggplot(as.data.frame(DiffAgnostic_250), aes(DiffAgnostic_250, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# 
# DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
#   geom_histogram(bins = 25) +
#   scale_x_continuous(limits = c(0, 0.5)) +
#   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_250 <- ggplot(ModelCompareAUCTrain_250, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_250 <- ggplot(ModelCompareAUCTest_250, aes(AUC, fill = modelType)) + 
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
# png("./Figs/SVM/SVM_BS_AUC_250_AgnPairs.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_250 / AUC_Test_DistrHist_250) | (DiffHist_Agnostic_250 / DiffHist_Mech) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
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
# AUC_Train_DistrHist_37 <- ggplot(ModelCompareAUCTrain_37, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (37 pairs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_37 <- ggplot(ModelCompareAUCTest_37, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (37 pairs) vs mechanistic SVM models in the testing data") + My_Theme2
# 
# #############################
# AUC_Train_DistrHist_50 <- ggplot(ModelCompareAUCTrain_50, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (50 pairs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_50 <- ggplot(ModelCompareAUCTest_50, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (50 pairs) vs mechanistic SVM models in the testing data") + My_Theme2
# 
# ##############################
# AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUCTrain_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (100 pairs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUCTest_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (100 pairs) vs mechanistic SVM models in the testing data") + My_Theme2
# 
# ###############################
# AUC_Train_DistrHist_250 <- ggplot(ModelCompareAUCTrain_250, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5, adjust=0.01) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (250 pairs) vs mechanistic SVM models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_250 <- ggplot(ModelCompareAUCTest_250, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (250 pairs) vs mechanistic SVM models in the testing data") + My_Theme2
# 
# 
# # All in one figure
# png("./Figs/BootStrap_Diff_No_features/GenePairs/SVM_BS_AUC_DiffPairs.png", width = 3000, height = 1500, res = 150)
# (AUC_Train_DistrHist_37 / AUC_Test_DistrHist_37) | (AUC_Train_DistrHist_50 / AUC_Test_DistrHist_50) | (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (AUC_Train_DistrHist_250 / AUC_Test_DistrHist_250) + plot_layout(widths = c(1, 1)) 
# dev.off()
# 
