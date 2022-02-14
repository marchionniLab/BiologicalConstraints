###################################################
## XGB

## with random genes - Wikum


rm(list = ls())

## Load necessary packages
library(xgboost)
library(pdp)
library(lime)
library(preprocessCore)
library(limma)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(mltools)
library(patchwork)
library(boot)
library(switchBox)

##########################################################
## XGB
# Mechanistic

## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic.rda")
load("./Objs/ProgressionDataGood2.rda")


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

## Parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1   # default = 0.3, range: [0,1]
  gamma              = 0,             #0  # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #1   # default = 1,   range: [0,∞]
  subsample          = 0.4,           #0.4     # default = 1,   range: (0,1]
  colsample_bytree   = 1,           #1      # default = 1,   range: (0,1]
  colsample_bylevel  = 1,    #1  # default = 1,   range: (0,1]
  lambda             = 0,             #0  # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

# The function for bootstraping
XGBStrap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  Progression <- sum(Train_label == 1)
  NoProgression <- sum(Train_label == 0)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = NoProgression/Progression)
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
bootobjectMech <- boot(data= Data_train_Mech, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15) 
AUCs_XG_Mech <- bootobjectMech$t
colnames(AUCs_XG_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

#########################################################################3
#########################################################################3

# Agnostic
## 74 random genes

## Load data
load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")

### Quantile normalize
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

Training <- t(usedTrainMat)

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
#names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

#names(usedValGroup) <- rownames(Validation) 
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
Testing <- t(usedTestMat)

#names(usedTestGroup) <- rownames(Testing)
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

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,           # 0.4      # default = 1,   range: (0,1]
  colsample_bytree   = 0.5,           # 0.5      # default = 1,   range: (0,1]
  colsample_bylevel  = 0.3,              # 0.3   # default = 1,   range: (0,1]
  lambda             = 3,             # 3 # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)


# The function for bootstraping
XGBStrap <- function(data, indices) {
  
  sel.genes = sample(setdiff(colnames(data),  "usedTrainGroup1"), 74)
  data = data[, c(sel.genes, "usedTrainGroup1")]
  
  # subset the validation data
  Validation <- Validation[ , c(sel.genes)]
  
  # subset the testing data
  Testing <- Testing[ , c(sel.genes)]
  
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  Progression <- sum(Train_label == 1)
  NoProgression <- sum(Train_label == 0)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  
  DataVal_Agnostic <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
  DataTest_Agnostic <- xgb.DMatrix(as.matrix(Testing), label = Test_label)
  watchlist <- list(train  = d2, test = DataVal_Agnostic)
  
  
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = NoProgression/Progression)
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
bootobjectAgnostic_74 <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R=1000, parallel="multicore", ncpus=15) 

#########################################################################3
#########################################################################3

# Agnostic
## 100 random genes

## Load data
load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")

### Quantile normalize
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

Training <- t(usedTrainMat)

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
#names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

#names(usedValGroup) <- rownames(Validation) 
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
Testing <- t(usedTestMat)

#names(usedTestGroup) <- rownames(Testing)
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

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,           # 0.4      # default = 1,   range: (0,1]
  colsample_bytree   = 0.5,           # 0.5      # default = 1,   range: (0,1]
  colsample_bylevel  = 0.3,              # 0.3   # default = 1,   range: (0,1]
  lambda             = 3,             # 3 # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)


# The function for bootstraping
XGBStrap <- function(data, indices) {
  
  sel.genes = sample(setdiff(colnames(data),  "usedTrainGroup1"), 100)
  data = data[, c(sel.genes, "usedTrainGroup1")]
  
  # subset the validation data
  Validation <- Validation[ , c(sel.genes)]
  
  # subset the testing data
  Testing <- Testing[ , c(sel.genes)]
  
  
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  Progression <- sum(Train_label == 1)
  NoProgression <- sum(Train_label == 0)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  
  DataVal_Agnostic <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
  DataTest_Agnostic <- xgb.DMatrix(as.matrix(Testing), label = Test_label)
  watchlist <- list(train  = d2, test = DataVal_Agnostic)
  
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = NoProgression/Progression)
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
bootobjectAgnostic_100 <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15) 

#########################################################################3
#########################################################################3

# Agnostic
## 200 random genes

## Load data
load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")

### Quantile normalize
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

Training <- t(usedTrainMat)

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
#names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

#names(usedValGroup) <- rownames(Validation) 
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
Testing <- t(usedTestMat)

#names(usedTestGroup) <- rownames(Testing)
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

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,           # 0.4      # default = 1,   range: (0,1]
  colsample_bytree   = 0.5,           # 0.5      # default = 1,   range: (0,1]
  colsample_bylevel  = 0.3,              # 0.3   # default = 1,   range: (0,1]
  lambda             = 3,             # 3 # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)


# The function for bootstraping
XGBStrap <- function(data, indices) {
  
  sel.genes = sample(setdiff(colnames(data),  "usedTrainGroup1"), 100)
  data = data[, c(sel.genes, "usedTrainGroup1")]
  
  # subset the validation data
  Validation <- Validation[ , c(sel.genes)]
  
  # subset the testing data
  Testing <- Testing[ , c(sel.genes)]
  
  
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  Progression <- sum(Train_label == 1)
  NoProgression <- sum(Train_label == 0)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  
  DataVal_Agnostic <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
  DataTest_Agnostic <- xgb.DMatrix(as.matrix(Testing), label = Test_label)
  watchlist <- list(train  = d2, test = DataVal_Agnostic)
  
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = NoProgression/Progression)
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
bootobjectAgnostic_200 <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15) 

#########################################################################3
#########################################################################3

# Agnostic
## 500 Random genes

## Load data
load("./Objs/ProgressionDataGood2.rda")
load("./Objs/Correlation/RGenes.rda")

### Quantile normalize
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

Training <- t(usedTrainMat)

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
#names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

#names(usedValGroup) <- rownames(Validation) 
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
Testing <- t(usedTestMat)

#names(usedTestGroup) <- rownames(Testing)
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

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,           # 0.4      # default = 1,   range: (0,1]
  colsample_bytree   = 0.5,           # 0.5      # default = 1,   range: (0,1]
  colsample_bylevel  = 0.3,              # 0.3   # default = 1,   range: (0,1]
  lambda             = 3,             # 3 # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)


# The function for bootstraping
XGBStrap <- function(data, indices) {
  
  sel.genes = sample(setdiff(colnames(data),  "usedTrainGroup1"), 100)
  data = data[, c(sel.genes, "usedTrainGroup1")]
  
  # subset the validation data
  Validation <- Validation[ , c(sel.genes)]
  
  # subset the testing data
  Testing <- Testing[ , c(sel.genes)]
  
  
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  Progression <- sum(Train_label == 1)
  NoProgression <- sum(Train_label == 0)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  
  DataVal_Agnostic <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
  DataTest_Agnostic <- xgb.DMatrix(as.matrix(Testing), label = Test_label)
  watchlist <- list(train  = d2, test = DataVal_Agnostic)
  
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = NoProgression/Progression)
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
bootobjectAgnostic_500 <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15) 

################################################################################################
################################################################################################

### Save boot objects
save(bootobjectMech, bootobjectAgnostic_74, bootobjectAgnostic_100, bootobjectAgnostic_200, bootobjectAgnostic_500, file = "./Objs/XGB/XGBBootObjects_RAND.rda")

# Load
load("./Objs/XGB/XGBBootObjects_RAND.rda")

################################################################################################
################################################################################################

## Work with Agnostic bootobject 74 vs mechanistic


AUCs_XGB_Agnostic_74 <- bootobjectAgnostic_74$t
colnames(AUCs_XGB_Agnostic_74) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_XG_Mech <- bootobjectMech$t
colnames(AUCs_XG_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_74 <- AUCs_XGB_Agnostic_74[, "AUC_Train"] - AUCs_XGB_Agnostic_74[, "AUC_Test"]
range(DiffAgnostic_74)
quantile(DiffAgnostic_74, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_XG_Mech[, "AUC_Train"] - AUCs_XG_Mech[, "AUC_Test"]
range(DiffMech)
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUC_Train <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Train"])
AgnosticAUC_Train_74 <- data.frame(AUC = AUCs_XGB_Agnostic_74[, "AUC_Train"])

MechanisticAUC_Train$modelType <- "Mechanistic"
AgnosticAUC_Train_74$modelType <- "Random genes"

ModelCompareAUC_Train_74 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_74)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUC_Test <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Test"])
AgnosticAUC_Test_74 <- data.frame(AUC = AUCs_XGB_Agnostic_74[, "AUC_Test"])

MechanisticAUC_Test$modelType <- "Mechanistic"
AgnosticAUC_Test_74$modelType <- "Random genes"

ModelCompareAUC_Test_74 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_74)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_74$data_type <- "Training"
ModelCompareAUC_Test_74$data_type <- "Testing"

ModelCompareAUC_Train_74$NofFeatAgn <- "74_Genes"
ModelCompareAUC_Test_74$NofFeatAgn <- "74_Genes"

save(ModelCompareAUC_Train_74, ModelCompareAUC_Test_74, file = "./Objs/XGB/ModelCompare_RAND_AUC_74.rda")

#########################################################################
## Save for the main figure
ModelCompare_XGB <- rbind(ModelCompareAUC_Train_74, ModelCompareAUC_Test_74)
ModelCompare_XGB$algorithm <- "XGB"
save(ModelCompare_XGB, file = "./Objs/XGB/ModelCompare_XGB_RAND.rda")
#########################################################################


##########
#######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# # DiffHist_Agnostic_74 <- ggplot(as.data.frame(DiffAgnostic_74), aes(DiffAgnostic_74, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# # 
# # DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_74 <- ggplot(ModelCompareAUC_Train_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_74 <- ggplot(ModelCompareAUC_Test_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the testing data") + My_Theme
# 
# 
# png("../../Figs/XGB/XGB_BS_RAND_AUC_74.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_74 / AUC_Test_DistrHist_74) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

################################################################################################
################################################################################################

## Work with Agnostic bootobject 100 vs mechanistic

AUCs_XGB_Agnostic_100 <- bootobjectAgnostic_100$t
colnames(AUCs_XGB_Agnostic_100) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_100 <- AUCs_XGB_Agnostic_100[, "AUC_Train"] - AUCs_XGB_Agnostic_100[, "AUC_Test"]
range(DiffAgnostic_100)
quantile(DiffAgnostic_100, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_100 <- data.frame(AUC = AUCs_XGB_Agnostic_100[, "AUC_Train"])

AgnosticAUC_Train_100$modelType <- "Random genes"

ModelCompareAUC_Train_100 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_100)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_100 <- data.frame(AUC = AUCs_XGB_Agnostic_100[, "AUC_Test"])

AgnosticAUC_Test_100$modelType <- "Random genes"

ModelCompareAUC_Test_100 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_100)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_100$data_type <- "Training"
ModelCompareAUC_Test_100$data_type <- "Testing"

ModelCompareAUC_Train_100$NofFeatAgn <- "100_Genes"
ModelCompareAUC_Test_100$NofFeatAgn <- "100_Genes"

save(ModelCompareAUC_Train_100, ModelCompareAUC_Test_100, file = "./Objs/XGB/ModelCompare_RAND_AUC_100.rda")

#######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# # DiffHist_Agnostic_100 <- ggplot(as.data.frame(DiffAgnostic_100), aes(DiffAgnostic_100, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# # 
# # DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUC_Train_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUC_Test_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the testing data") + My_Theme
# 
# png("../../Figs/XGB/XGB_BS_RAND_AUC_100.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

################################################################################################
################################################################################################

## Work with Agnostic bootobject 200 vs mechanistic

AUCs_XGB_Agnostic_200 <- bootobjectAgnostic_200$t
colnames(AUCs_XGB_Agnostic_200) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_200 <- AUCs_XGB_Agnostic_200[, "AUC_Train"] - AUCs_XGB_Agnostic_200[, "AUC_Test"]
range(DiffAgnostic_200)
quantile(DiffAgnostic_200, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_200 <- data.frame(AUC = AUCs_XGB_Agnostic_200[, "AUC_Train"])

AgnosticAUC_Train_200$modelType <- "Random_genes"

ModelCompareAUC_Train_200 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_200)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_200 <- data.frame(AUC = AUCs_XGB_Agnostic_200[, "AUC_Test"])

AgnosticAUC_Test_200$modelType <- "Random_genes"

ModelCompareAUC_Test_200 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_200)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_200$data_type <- "Training"
ModelCompareAUC_Test_200$data_type <- "Testing"

ModelCompareAUC_Train_200$NofFeatAgn <- "200_Genes"
ModelCompareAUC_Test_200$NofFeatAgn <- "200_Genes"

save(ModelCompareAUC_Train_200, ModelCompareAUC_Test_200, file = "./Objs/XGB/ModelCompare_RAND_AUC_200.rda")

#######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# # DiffHist_Agnostic_200 <- ggplot(as.data.frame(DiffAgnostic_200), aes(DiffAgnostic_200, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# # 
# # DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUC_Train_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUC_Test_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the testing data") + My_Theme
# 
# png("../../Figs/XGB/XGB_BS_AUC_200.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

################################################################################################
################################################################################################

## Work with Agnostic bootobject 500 vs mechanistic

AUCs_XGB_Agnostic_500 <- bootobjectAgnostic_500$t
colnames(AUCs_XGB_Agnostic_500) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_500 <- AUCs_XGB_Agnostic_500[, "AUC_Train"] - AUCs_XGB_Agnostic_500[, "AUC_Test"]
range(DiffAgnostic_500)
quantile(DiffAgnostic_500, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_500 <- data.frame(AUC = AUCs_XGB_Agnostic_500[, "AUC_Train"])

AgnosticAUC_Train_500$modelType <- "Random genes"

ModelCompareAUC_Train_500 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_500)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_500 <- data.frame(AUC = AUCs_XGB_Agnostic_500[, "AUC_Test"])

AgnosticAUC_Test_500$modelType <- "Random genes"

ModelCompareAUC_Test_500 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_500)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_500$data_type <- "Training"
ModelCompareAUC_Test_500$data_type <- "Testing"

ModelCompareAUC_Train_500$NofFeatAgn <- "500_Genes"
ModelCompareAUC_Test_500$NofFeatAgn <- "500_Genes"

save(ModelCompareAUC_Train_500, ModelCompareAUC_Test_500, file = "./Objs/XGB/ModelCompare_RAND_AUC_500.rda")

###############################
## save all

ModelCompare_XGB_RAND_DiffNoFeat <- rbind(ModelCompareAUC_Train_74,
                                          ModelCompareAUC_Test_74,
                                          ModelCompareAUC_Train_100,
                                          ModelCompareAUC_Test_100,
                                          ModelCompareAUC_Train_200,
                                          ModelCompareAUC_Test_200,
                                          ModelCompareAUC_Train_500,
                                          ModelCompareAUC_Test_500
)

save(ModelCompare_XGB_RAND_DiffNoFeat, file = "./Objs/XGB/ModelCompare_XGB_RAND_DiffNoFeat.rda")

#######
## Plots
# My_Theme = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=8)
# )
# 
# # DiffHist_Agnostic_500 <- ggplot(as.data.frame(DiffAgnostic_500), aes(DiffAgnostic_500, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme
# # 
# # DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
# #   geom_histogram(bins = 25) +
# #   scale_x_continuous(limits = c(0, 0.5)) +
# #   labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme
# 
# AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUC_Train_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the training data") + My_Theme
# 
# AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUC_Test_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="AUC distribution of the agnostic and mechanistic XGB models in the testing data") + My_Theme
# 
# png("../../Figs/XGB/XGB_BS_RAND_AUC_500.png", width = 3000, height = 1500, res = 300)
# (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) + plot_layout(widths = c(2, 1)) #| BarPlotImportance
# dev.off()

###################################################################################
###########################################
###########################################
# My_Theme2 = theme(
#   axis.title.x = element_text(size = 5),
#   axis.text.x = element_text(size = 5),
#   axis.title.y = element_text(size = 5),
#   axis.text.y = element_text(size = 5),
#   plot.title = element_text(size=9)
# )
# 
# #############################
# AUC_Train_DistrHist_74 <- ggplot(ModelCompareAUC_Train_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (74 random genes) vs mechanistic XGB models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_74 <- ggplot(ModelCompareAUC_Test_74, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (74 random genes) vs mechanistic XGB models in the testing data") + My_Theme2
# #############################
# AUC_Train_DistrHist_100 <- ggplot(ModelCompareAUC_Train_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (100 random genes) vs mechanistic XGB models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_100 <- ggplot(ModelCompareAUC_Test_100, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (100 random genes) vs mechanistic XGB models in the testing data") + My_Theme2
# ##############################
# AUC_Train_DistrHist_200 <- ggplot(ModelCompareAUC_Train_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (200 random genes) vs mechanistic XGB models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_200 <- ggplot(ModelCompareAUC_Test_200, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (200 random genes) vs mechanistic XGB models in the testing data") + My_Theme2
# ###############################
# AUC_Train_DistrHist_500 <- ggplot(ModelCompareAUC_Train_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (500 random genes) vs mechanistic XGB models in the training data") + My_Theme2
# 
# AUC_Test_DistrHist_500 <- ggplot(ModelCompareAUC_Test_500, aes(AUC, fill = modelType)) + 
#   geom_density(alpha = 0.5) +
#   scale_x_continuous(limits = c(0.5, 1)) +
#   labs(title="Agnostic (500 random genes) vs mechanistic XGB models in the testing data") + My_Theme2
# 
# 
# # All in one figure
# png("../../Figs/BootStrap_Diff_No_features/XGB_BS_RAND_AUC_DiffFeatures.png", width = 3000, height = 1500, res = 150)
# (AUC_Train_DistrHist_74 / AUC_Test_DistrHist_74) | (AUC_Train_DistrHist_100 / AUC_Test_DistrHist_100) | (AUC_Train_DistrHist_200 / AUC_Test_DistrHist_200) | (AUC_Train_DistrHist_500 / AUC_Test_DistrHist_500) + plot_layout(widths = c(1, 1)) 
# dev.off()
# 
