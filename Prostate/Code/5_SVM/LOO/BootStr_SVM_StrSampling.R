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


# Agnostic

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train_Agnostic <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
#control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

names(Data_train_Agnostic) <- make.names(names(Data_train_Agnostic))

colnames(Training) <- make.names(colnames(Training))
colnames(Testing) <- make.names(colnames(Testing))

Grid <- expand.grid(degree = 3, scale = 0.01, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
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
bootobjectAgnostic <- boot(data= Data_train_Agnostic, statistic= SVM_Strap, R= 1000) 

AUCs_SVM_Agnostic <- bootobjectAgnostic$t
colnames(AUCs_SVM_Agnostic) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

###################################################
## SVM
# Mechanistic
## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")


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

Grid <- expand.grid(degree = 3, scale = 0.001, C = 0.25)

# The function for bootstraping
SVM_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  SVM <- train(usedTrainGroup~., data=d, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
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
bootobjectMech <- boot(data= Data_train_Mechanistic, statistic= SVM_Strap, R= 1000) 

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

save(bootobjectAgnostic, bootobjectMech, file= "./Objs/SVM/SVMBootObjects.rda")

load("./Objs/SVM/SVMBootObjects.rda")

AUCs_SVM_Agnostic <- bootobjectAgnostic$t
colnames(AUCs_SVM_Agnostic) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_SVM_Mech <- bootobjectMech$t
colnames(AUCs_SVM_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffTrain <- AUCs_SVM_Mech[, "AUC_Train"] - AUCs_SVM_Agnostic[, "AUC_Train"]
quantile(DiffTrain, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
DiffTest <- AUCs_SVM_Mech[, "AUC_Test"] - AUCs_SVM_Agnostic[, "AUC_Test"]
quantile(DiffTest, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUCTrain <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Train"])
AgnosticAUCTrain <- data.frame(AUC = AUCs_SVM_Agnostic[, "AUC_Train"])

MechanisticAUCTrain$modelType <- "Mech"
AgnosticAUCTrain$modelType <- "Agnostic"

ModelCompareAUCTrain <- rbind(MechanisticAUCTrain, AgnosticAUCTrain)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUCTest <- data.frame(AUC = AUCs_SVM_Mech[, "AUC_Test"])
AgnosticAUCTest <- data.frame(AUC = AUCs_SVM_Agnostic[, "AUC_Test"])

MechanisticAUCTest$modelType <- "Mech"
AgnosticAUCTest$modelType <- "Agnostic"

ModelCompareAUCTest <- rbind(MechanisticAUCTest, AgnosticAUCTest)

##########
## The N of important features (output)
ImportanceAgnostic <- data.frame(NofOutputFeatures = AUCs_SVM_Agnostic[, "N_ImportanVariables"])
ImportanceMech <- data.frame(NofOutputFeatures = AUCs_SVM_Mech[, "N_ImportanVariables"])

ImportanceAgnostic$modelType <- "Agnostic"
ImportanceMech$modelType <- "Mech"

ModelCompare_ImportantFeatures <- rbind(ImportanceAgnostic, ImportanceMech)

######
## Plots

My_Theme = theme(
  axis.title.x = element_text(size = 7),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 7),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=9)
)

DiffHistTrain <- ggplot(as.data.frame(DiffTrain), aes(DiffTrain, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(-0.25, 0.25)) +
  labs(title="Histogram of the difference in the training data") + My_Theme

DiffHistTest <- ggplot(as.data.frame(DiffTest), aes(DiffTest, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(-0.25, 0.25)) +
  labs(title="Histogram of the difference in the testing data") + My_Theme

AUC_Train_DistrHist <- ggplot(ModelCompareAUCTrain, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.6, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic SVM models in the training data") + My_Theme

AUC_Test_DistrHist <- ggplot(ModelCompareAUCTest, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.6, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic SVM models in the testing data") + My_Theme

BarPlotImportance <- ggplot(data=ModelCompare_ImportantFeatures, aes(x=NofOutputFeatures, color = modelType, fill = modelType)) +
  geom_bar(stat="bin") +
  scale_x_continuous(limits = c(90, 105)) +
  labs(title="Distribution of the number of important features (output)") + 
  facet_grid(~modelType, scale='free_x') +
  My_Theme

png("./Figs/SVM/SVM_BS_AUC_Compare_WithDiff.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist / AUC_Test_DistrHist) | (DiffHistTrain / DiffHistTest) + plot_layout(widths = c(1, 1)) | BarPlotImportance
dev.off()


