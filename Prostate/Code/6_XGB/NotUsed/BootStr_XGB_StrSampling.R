###################################################
## XGB

rm(list = ls())

## Load necessary packages
library(xgboost)
library(pdp)
library(lime)
library(preprocessCore)
library(limma)
library(pROC)
library(caret)
library(DiagrammeR)
library(ggplot2)
library(xgboostExplainer)
library(dplyr)
library(Ckmeans.1d.dp)
library(mltools)
library(patchwork)
library(boot)
# Agnostic

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic)

#####################
## Here, WE divide the training data into "actual training" and another validation
# This validation data will be used in the "WatchList" parameter. It is independent of the testing data.
set.seed(333)
inds <- createDataPartition(usedTrainGroup, p = 0.7, list = F)
Training1 <- Training[inds, ]
Validation <- Training[-inds, ]

usedTrainGroup1 <- usedTrainGroup[inds]
usedValGroup <- usedTrainGroup[-inds]

table(usedTrainGroup1)
table(usedValGroup)



## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

all(rownames(Validation) == names(usedValGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training1)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train_Agnostic <- cbind(Training, usedTrainGroup1)

## The same for validation
Validation <- as.data.frame(Validation)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_val_Agnostic <- cbind(Validation, usedValGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train_Agnostic$usedTrainGroup1)  
levels(Data_train_Agnostic$usedTrainGroup1) <- c(0,1) 
Data_train_Agnostic$usedTrainGroup1 <- factor(Data_train_Agnostic$usedTrainGroup1, levels = c(0,1)) # 0=No,1= Yes 
Train_label <- Data_train_Agnostic$usedTrainGroup1
Train_label <- as.vector(Train_label)
table(Train_label)

## The same for validation
table(Data_val_Agnostic$usedValGroup)  
levels(Data_val_Agnostic$usedValGroup) <- c(0,1) 
Data_val_Agnostic$usedValGroup <- factor(Data_val_Agnostic$usedValGroup, levels = c(0,1)) # 0=No,1= Yes 
Val_label <- Data_val_Agnostic$usedValGroup
Val_label <- as.vector(Val_label)
table(Val_label)


## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test_Agnostic <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_test_Agnostic$usedTestGroup)  
levels(Data_test_Agnostic$usedTestGroup) <- c(0,1)
Data_test_Agnostic$usedTestGroup <- factor(Data_test_Agnostic$usedTestGroup, levels = c(0,1)) #0=No, 1=Yes
Test_label <- Data_test_Agnostic$usedTestGroup
Test_label <- as.vector(Test_label)
table(Test_label)



## Convert to xgb.DMatrix
DataTrain_Agnostic <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataVal_Agnostic <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
DataTest_Agnostic <- xgb.DMatrix(as.matrix(Testing), label = Test_label)

## Creat a watch list
watchlist <- list(train  = DataTrain_Agnostic, test = DataVal_Agnostic)

##########################
# Scale weight (to compensate for un-balanced class sizes)
Mets <- sum(Train_label == 1)
No_Mets <- sum(Train_label == 0)


## Parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 3,             # 3
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,           #0.5 /0.4    # default = 1,   range: (0,1]
  colsample_bytree   = 1,              #1   # default = 1,   range: (0,1]
  colsample_bylevel  = 1,               #1   # default = 1,   range: (0,1]
  lambda             = 1,             # 0.5 /2  # default = 1
  alpha              = 0,           # 0.5 /0   # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

# The function for bootstraping
XGBStrap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = No_Mets/Mets)
  Importance_Agnostic <- xgb.importance(model = XGB)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic$Gain, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic$Gain > 0, ]
  N_ImportanVariables <- length(Importance_Agnostic$Gain)
  train_preds <- predict(XGB, d2)
  test_preds <- predict(XGB, DataTest_Agnostic)
  ROCTrainAgnostic <- roc(Train_label, train_preds, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(Test_label, test_preds, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R= 1000) 

AUCs_XGB_Agnostic <- bootobjectAgnostic$t

colnames(AUCs_XGB_Agnostic) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

###################################################
## XGB
# Mechanistic
## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Mechanistic)

#####################
## Here, WE divide the training data into "actual training" and another validation
# This validation data will be used in the "WatchList" parameter. It is independent of the testing data.
set.seed(333)
inds <- createDataPartition(usedTrainGroup, p = 0.7, list = F)
Training1 <- Training[inds, ]
Validation <- Training[-inds, ]

usedTrainGroup1 <- usedTrainGroup[inds]
usedValGroup <- usedTrainGroup[-inds]

table(usedTrainGroup1)
table(usedValGroup)


## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

all(rownames(Validation) == names(usedValGroup))

## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training1)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train_Mech <- cbind(Training, usedTrainGroup1)

## The same for validation
Validation <- as.data.frame(Validation)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_val_Mech <- cbind(Validation, usedValGroup)


########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Mechanistic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train_Mech$usedTrainGroup1)  
levels(Data_train_Mech$usedTrainGroup1) <- c(0,1) 
Data_train_Mech$usedTrainGroup1 <- factor(Data_train_Mech$usedTrainGroup1, levels = c(0,1)) # 0=No,1= Yes 
Train_label <- Data_train_Mech$usedTrainGroup1
Train_label <- as.vector(Train_label)
table(Train_label)


## The same for validation
table(Data_val_Mech$usedValGroup)  
levels(Data_val_Mech$usedValGroup) <- c(0,1) 
Data_val_Mech$usedValGroup <- factor(Data_val_Mech$usedValGroup, levels = c(0,1)) # 0=No,1= Yes 
Val_label <- Data_val_Mech$usedValGroup
Val_label <- as.vector(Val_label)
table(Val_label)

## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test_Mech <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_test_Mech$usedTestGroup)  
levels(Data_test_Mech$usedTestGroup) <- c(0,1)
Data_test_Mech$usedTestGroup <- factor(Data_test_Mech$usedTestGroup, levels = c(0,1)) #0=No, 1=Yes
Test_label <- Data_test_Mech$usedTestGroup
Test_label <- as.vector(Test_label)
table(Test_label)


#save(Train_label, Test_label, file = "./Objs/XGB/Labels.rda")

## Convert to xgb.DMatrix
DataTrain_Mech <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataVal_Mech <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
DataTest_Mech <- xgb.DMatrix(as.matrix(Testing), label = Test_label)

## Creat a watch list
watchlist <- list(train  = DataTrain_Mech, test = DataVal_Mech)

##########################
# Scale weight (to compensate for un-balanced class sizes)
Mets <- sum(Train_label == 1)
No_Mets <- sum(Train_label == 0)


## Parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.01,           #0.1 /0.01 /0.1/0.01   # default = 0.3, range: [0,1]
  gamma              = 0,             #0  # default = 0,   range: [0,∞]
  max_depth          = 3,             # 3 /6 /3
  min_child_weight   = 1,             #1   # default = 1,   range: [0,∞]
  subsample          = 0.3,           #0.3   # default = 1,   range: (0,1]
  colsample_bytree   = 1,           #1      # default = 1,   range: (0,1]
  colsample_bylevel  = 1,    #1  # default = 1,   range: (0,1]
  lambda             = 2,             #0.5 /1 /0.5  # default = 1
  alpha              = 0,           # 0.5 /1/0/0.1   # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

# The function for bootstraping
XGBStrap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = No_Mets/Mets)
  Importance_Mech <- xgb.importance(model = XGB)
  Importance_Mech <- Importance_Mech[order(Importance_Mech$Gain, decreasing = TRUE), ]
  Importance_Mech <- Importance_Mech[Importance_Mech$Gain > 0, ]
  N_ImportanVariables <- length(Importance_Mech$Gain)
  train_preds <- predict(XGB, d2)
  test_preds <- predict(XGB, DataTest_Mech)
  ROCTrainMech <- roc(Train_label, train_preds, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestMech <- roc(Test_label, test_preds, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainMech$auc, ROCTestMech$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectMech <- boot(data= Data_train_Mech, statistic= XGBStrap, R= 1000) 
AUCs_XG_Mech <- bootobjectMech$t

colnames(AUCs_XG_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

#save(bootobjectAgnostic, bootobjectMech, file = "./Objs/XGB/XGBBootObjects.rda")

load("./Objs/XGB/XGBBootObjects.rda")

AUCs_XGB_Agnostic <- bootobjectAgnostic$t
AUCs_XG_Mech <- bootobjectMech$t
colnames(AUCs_XGB_Agnostic) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")
colnames(AUCs_XG_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
Diff_Train <- AUCs_XG_Mech[,"AUC_Train"] - AUCs_XGB_Agnostic[,"AUC_Train"]
range(Diff_Train)
quantile(Diff_Train, c(0.025, 0.975))
#colnames(Diff_Train) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
Diff_Test <- AUCs_XG_Mech[,"AUC_Test"] - AUCs_XGB_Agnostic[,"AUC_Test"]
range(Diff_Test)
quantile(Diff_Test, c(0.025, 0.975))
#colnames(Diff_Test) <- "Diff"


## Plot the distributions of the AUCs from both methods i the training data
MechanisticTrain <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Train"])
AgnosticTrain <- data.frame(AUC = AUCs_XGB_Agnostic[, "AUC_Train"])

MechanisticTrain$modelType <- "Mech"
AgnosticTrain$modelType <- "Agnostic"

ModelCompareTrain <- rbind(MechanisticTrain, AgnosticTrain)

## Plot the distributions of the AUCs from both methods i the testing data
MechanisticTest <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Test"])
AgnosticTest <- data.frame(AUC = AUCs_XGB_Agnostic[, "AUC_Test"])

MechanisticTest$modelType <- "Mech"
AgnosticTest$modelType <- "Agnostic"

ModelCompareTest <- rbind(MechanisticTest, AgnosticTest)

#######
## Plots

My_Theme = theme(
  axis.title.x = element_text(size = 7),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 7),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=9)
)

DiffHistTrain <- ggplot(as.data.frame(Diff_Train), aes(Diff_Train, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(-0.25, 0.25)) +
  labs(title="The difference in the training data") + My_Theme

DiffHistTest <- ggplot(as.data.frame(Diff_Test), aes(Diff_Test, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(-0.25, 0.25)) +
  labs(title="The difference in the testing data") + My_Theme


AUCDistrHistTrain <- ggplot(ModelCompareTrain, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.6, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic XGB models in the training data") + My_Theme


AUCDistrHistTest <- ggplot(ModelCompareTest, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.6, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic XGB models in the testing data") + My_Theme

png("./Figs/XGB/XGB_BS_AUC_Compare_WithDiff.png", width = 3000, height = 1500, res = 300)
(AUCDistrHistTrain / AUCDistrHistTest) | (DiffHistTrain / DiffHistTest) + plot_layout(widths = c(1, 1)) 
dev.off()



