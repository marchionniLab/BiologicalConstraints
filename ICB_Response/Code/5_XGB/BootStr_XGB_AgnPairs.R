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
library(dplyr)
library(Ckmeans.1d.dp)
library(mltools)
library(patchwork)

#########################
## XGB
# Mechanistic
## Load data
# load("./Objs/KTSP/KTSP_STATs_Mechanistic.rda")
# load("./Objs/icbData.rda")
# 
# ### Associated groups
# usedTrainGroup <- mixTrainGroup
# usedTestGroup <- mixTestGroup
# 
# ### Transpose usedTrainMat (making samples as rows instead of columns)
# Training <- t(KTSP_STATs_Train_Mechanistic)
# 
# #####################
# ## Here, WE divide the training data into "actual training" and another validation
# # This validation data will be used in the "WatchList" parameter. It is independent of the testing data.
# set.seed(333)
# inds <- createDataPartition(usedTrainGroup, p = 0.7, list = F)
# Training1 <- Training[inds, ]
# Validation <- Training[-inds, ]
# 
# usedTrainGroup1 <- usedTrainGroup[inds]
# usedValGroup <- usedTrainGroup[-inds]
# 
# table(usedTrainGroup1)
# table(usedValGroup)
# 
# 
# ## Making sure that sample names are identical in both Training and usedTrainGroup
# names(usedTrainGroup1) <- rownames(Training1)
# all(rownames(Training1) == names(usedTrainGroup1))
# 
# all(rownames(Validation) == names(usedValGroup))
# 
# ## Combining the expression matrix and the phenotype in one data frame
# Training <- as.data.frame(Training1)
# #usedTrainGroup <- as.data.frame(usedTrainGroup)
# Data_train_Mech <- cbind(Training, usedTrainGroup1)
# 
# ## The same for validation
# Validation <- as.data.frame(Validation)
# #usedTrainGroup <- as.data.frame(usedTrainGroup)
# Data_val_Mech <- cbind(Validation, usedValGroup)
# 
# 
# ########################################################
# # Transpose usedTestMat and make the sample names identical 
# Testing <- t(KTSP_STATs_Test_Mechanistic)
# 
# names(usedTestGroup) <- rownames(Testing)
# all(rownames(Testing) == names(usedTestGroup))
# 
# ###########################################################
# ## Converting classes from Progression/NoProgression Format to 0-1 Format
# table(Data_train_Mech$usedTrainGroup1)  
# levels(Data_train_Mech$usedTrainGroup1) <- c(0,1) 
# Data_train_Mech$usedTrainGroup1 <- factor(Data_train_Mech$usedTrainGroup1, levels = c(0,1)) # 0=No,1= Yes 
# Train_label <- Data_train_Mech$usedTrainGroup1
# Train_label <- as.vector(Train_label)
# table(Train_label)
# 
# 
# ## The same for validation
# table(Data_val_Mech$usedValGroup)  
# levels(Data_val_Mech$usedValGroup) <- c(0,1) 
# Data_val_Mech$usedValGroup <- factor(Data_val_Mech$usedValGroup, levels = c(0,1)) # 0=No,1= Yes 
# Val_label <- Data_val_Mech$usedValGroup
# Val_label <- as.vector(Val_label)
# table(Val_label)
# 
# ## Combine both the Expression matrix and the phenotype into one matrix
# Testing <- as.data.frame(Testing)
# Data_test_Mech <- cbind(Testing, usedTestGroup)
# 
# ## Converting classes from Progression/NoProgression Format to 0-1 Format
# table(Data_test_Mech$usedTestGroup)  
# levels(Data_test_Mech$usedTestGroup) <- c(0,1)
# Data_test_Mech$usedTestGroup <- factor(Data_test_Mech$usedTestGroup, levels = c(0,1)) #0=No, 1=Yes
# Test_label <- Data_test_Mech$usedTestGroup
# Test_label <- as.vector(Test_label)
# table(Test_label)
# 
# 
# #save(Train_label, Test_label, file = "./Objs/XGB/Labels.rda")
# 
# ## Convert to xgb.DMatrix
# DataTrain_Mech <- xgb.DMatrix(as.matrix(Training), label = Train_label)
# DataVal_Mech <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
# DataTest_Mech <- xgb.DMatrix(as.matrix(Testing), label = Test_label)
# 
# ## Creat a watch list
# watchlist <- list(train  = DataTrain_Mech, test = DataVal_Mech)
# 
# ##########################
# # Scale weight (to compensate for un-balanced class sizes)
# NR <- sum(Train_label == 0)
# R <- sum(Train_label == 1)
# 
# ## Parameters
# set.seed(333)
# 
# parameters <- list(
#   # General Parameters
#   booster            = "gbtree",          # default = "gbtree"
#   silent             = 1,                 # default = 0
#   # Booster Parameters
#   eta                = 0.001,           #0.3  # default = 0.3, range: [0,1]
#   gamma              = 0,             #0   # default = 0,   range: [0,∞]
#   max_depth          = 2,             # 2 
#   min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
#   subsample          = 0.6,             #0.4  # default = 1,   range: (0,1]
#   colsample_bytree   = 1,             # 1   # default = 1,   range: (0,1]
#   colsample_bylevel  = 1,             #1    # default = 1,   range: (0,1]
#   lambda             = 0,             #    # default = 1
#   alpha              = 0,           # 0.1    # default = 0
#   # Task Parameters
#   objective          = "binary:logistic",   # default = "reg:linear"
#   eval_metric        = "auc"
# )
# 
# # The function for bootstraping
# XGBStrap <- function(data, indices) {
#   d <- data[indices, ] # allows boot to select sample
#   Train_label <- as.integer(d$usedTrainGroup1)-1
#   # Scale weight (to compensate for un-balanced class sizes)
#   NR <- sum(Train_label == 0)
#   R <- sum(Train_label == 1)
#   d$usedTrainGroup1 <- NULL
#   d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
#   XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = R/NR)
#   Importance_Mech <- xgb.importance(model = XGB)
#   Importance_Mech <- Importance_Mech[order(Importance_Mech$Gain, decreasing = TRUE), ]
#   Importance_Mech <- Importance_Mech[Importance_Mech$Gain > 0, ]
#   N_ImportanVariables <- length(Importance_Mech$Gain)
#   train_preds <- predict(XGB, d2)
#   test_preds <- predict(XGB, DataTest_Mech)
#   ROCTrainMech <- roc(Train_label, train_preds, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
#   ROCTestMech <- roc(Test_label, test_preds, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
#   return(c(ROCTrainMech$auc, ROCTestMech$auc, N_ImportanVariables))
# }
# 
# 
# set.seed(333)
# bootobjectMech <- boot(data= Data_train_Mech, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15) 
# AUCs_XG_Mech <- bootobjectMech$t
# colnames(AUCs_XG_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

load('./objs/XGB/XGB_bootobjectMech.rda')

#########################################################################3
#########################################################################3

# Agnostic

# 25 pairs

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_25.rda")
load("./Objs/icbData.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_25)

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
Testing <- t(KTSP_STATs_Test_Agnostic_25)

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
NR <- sum(Train_label == 0)
R <- sum(Train_label == 1)


## Parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.3    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 1,             # 2
  min_child_weight   = 2,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.5,             #0.5  # default = 1,   range: (0,1]
  colsample_bytree   = 1,             #1    # default = 1,   range: (0,1]
  colsample_bylevel  = 1,             #1    # default = 1,   range: (0,1]
  lambda             = 0,             # 0  # default = 1
  alpha              = 0,             # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

# The function for bootstraping
XGBStrap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  NR <- sum(Train_label == 0)
  R <- sum(Train_label == 1)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = R/NR)
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
bootobjectAgnostic_25 <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15)

###################################################
#########################################################################3
#########################################################################3

# Agnostic

# 50 pairs

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_50.rda")
load("./Objs/icbData.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_50)

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
Testing <- t(KTSP_STATs_Test_Agnostic_50)

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
NR <- sum(Train_label == 0)
R <- sum(Train_label == 1)


## Parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.3    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 2,             # 2
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.5,             #0.5  # default = 1,   range: (0,1]
  colsample_bytree   = 1,             #1    # default = 1,   range: (0,1]
  colsample_bylevel  = 1,             #1    # default = 1,   range: (0,1]
  lambda             = 0,             # 0  # default = 1
  alpha              = 0,             # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

# The function for bootstraping
XGBStrap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  NR <- sum(Train_label == 0)
  R <- sum(Train_label == 1)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = R/NR)
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
bootobjectAgnostic_50 <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15) 

#########################################################################3
#########################################################################3

# Agnostic

# 100 pairs

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_100.rda")
load("./Objs/icbData.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_100)

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
Testing <- t(KTSP_STATs_Test_Agnostic_100)

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
NR <- sum(Train_label == 0)
R <- sum(Train_label == 1)


## Parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.3    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 2,             # 2
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.5,             #0.5  # default = 1,   range: (0,1]
  colsample_bytree   = 1,             #1    # default = 1,   range: (0,1]
  colsample_bylevel  = 1,             #1    # default = 1,   range: (0,1]
  lambda             = 0,             # 0  # default = 1
  alpha              = 0,             # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

# The function for bootstraping
XGBStrap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  NR <- sum(Train_label == 0)
  R <- sum(Train_label == 1)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = R/NR)
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

# 250 pairs

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_250.rda")
load("./Objs/icbData.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_250)

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
Testing <- t(KTSP_STATs_Test_Agnostic_250)

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
NR <- sum(Train_label == 0)
R <- sum(Train_label == 1)


## Parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.3    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 3,             # 2
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.5,             #0.5  # default = 1,   range: (0,1]
  colsample_bytree   = 1,             #1    # default = 1,   range: (0,1]
  colsample_bylevel  = 1,             #1    # default = 1,   range: (0,1]
  lambda             = 0,             # 0  # default = 1
  alpha              = 0,             # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

# The function for bootstraping
XGBStrap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  Train_label <- as.integer(d$usedTrainGroup1)-1
  # Scale weight (to compensate for un-balanced class sizes)
  NR <- sum(Train_label == 0)
  R <- sum(Train_label == 1)
  d$usedTrainGroup1 <- NULL
  d2 <- xgb.DMatrix(as.matrix(d), label = Train_label)
  XGB <- xgb.train(parameters, data = d2, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = R/NR)
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
bootobjectAgnostic_250 <- boot(data= Data_train_Agnostic, statistic= XGBStrap, R= 1000, parallel = "multicore", ncpus = 15) 

################################################################################################
################################################################################################

### Save boot objects
save(bootobjectMech, bootobjectAgnostic_25, bootobjectAgnostic_50, bootobjectAgnostic_100, bootobjectAgnostic_250, file = "./Objs/XGB/XGBBootObjects_AgnPairs_new.rda")

# Load
load("./Objs/XGB/XGBBootObjects_AgnPairs.rda")

################################################################################################
################################################################################################

## Work with Agnostic bootobject 25 vs mechanistic

AUCs_XGB_Agnostic_25 <- bootobjectAgnostic_25$t
colnames(AUCs_XGB_Agnostic_25) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_XG_Mech <- bootobjectMech$t
colnames(AUCs_XG_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_25 <- AUCs_XGB_Agnostic_25[, "AUC_Train"] - AUCs_XGB_Agnostic_25[, "AUC_Test"]
range(DiffAgnostic_25)
quantile(DiffAgnostic_25, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_XG_Mech[, "AUC_Train"] - AUCs_XG_Mech[, "AUC_Test"]
range(DiffMech)
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUC_Train <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Train"])
AgnosticAUC_Train_25 <- data.frame(AUC = AUCs_XGB_Agnostic_25[, "AUC_Train"])

MechanisticAUC_Train$modelType <- "Mechanistic"
AgnosticAUC_Train_25$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_25 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_25)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUC_Test <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Test"])
AgnosticAUC_Test_25 <- data.frame(AUC = AUCs_XGB_Agnostic_25[, "AUC_Test"])

MechanisticAUC_Test$modelType <- "Mechanistic"
AgnosticAUC_Test_25$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_25 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_25)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_25$data_type <- "Training"
ModelCompareAUC_Test_25$data_type <- "Testing"

ModelCompareAUC_Train_25$NofFeatAgn <- "25_Pairs"
ModelCompareAUC_Test_25$NofFeatAgn <- "25_Pairs"

####################################################################
## For the main figure
ModelCompare_XGB <- rbind(ModelCompareAUC_Train_25, ModelCompareAUC_Test_25)
ModelCompare_XGB$algorithm <- "XGB"
save(ModelCompare_XGB, file = "./Objs/XGB/ModelCompare_XGB_AgnPairs_new.rda")

####################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/XGB/ModelCompareAUC_50_new.rda")

# Combine
ModelCompareAUC_Train_25_50 <- rbind(ModelCompareAUC_Train_25, ModelCompareAUC_Train_50)
ModelCompareAUC_Test_25_50 <- rbind(ModelCompareAUC_Test_25, ModelCompareAUC_Test_50)

################################################################################################
################################################################################################

## Work with Agnostic bootobject 50 vs mechanistic

AUCs_XGB_Agnostic_50 <- bootobjectAgnostic_50$t
colnames(AUCs_XGB_Agnostic_50) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_XG_Mech <- bootobjectMech$t
colnames(AUCs_XG_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_50 <- AUCs_XGB_Agnostic_50[, "AUC_Train"] - AUCs_XGB_Agnostic_50[, "AUC_Test"]
range(DiffAgnostic_50)
quantile(DiffAgnostic_50, c(0.025, 0.975))

# Calculate the difference and the CI of the difference in the testing data
DiffMech <- AUCs_XG_Mech[, "AUC_Train"] - AUCs_XG_Mech[, "AUC_Test"]
range(DiffMech)
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUC_Train <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Train"])
AgnosticAUC_Train_50 <- data.frame(AUC = AUCs_XGB_Agnostic_50[, "AUC_Train"])

MechanisticAUC_Train$modelType <- "Mechanistic"
AgnosticAUC_Train_50$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_50 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_50)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUC_Test <- data.frame(AUC = AUCs_XG_Mech[, "AUC_Test"])
AgnosticAUC_Test_50 <- data.frame(AUC = AUCs_XGB_Agnostic_50[, "AUC_Test"])

MechanisticAUC_Test$modelType <- "Mechanistic"
AgnosticAUC_Test_50$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_50 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_50)

## add info about the data (training / testing) and the number of pairs
ModelCompareAUC_Train_50$data_type <- "Training"
ModelCompareAUC_Test_50$data_type <- "Testing"

ModelCompareAUC_Train_50$NofFeatAgn <- "50_Pairs"
ModelCompareAUC_Test_50$NofFeatAgn <- "50_Pairs"

####################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/XGB/ModelCompareAUC_100_new.rda")

ModelCompareAUC_Train_50_100 <- rbind(ModelCompareAUC_Train_50, ModelCompareAUC_Train_100)
ModelCompareAUC_Test_50_100 <- rbind(ModelCompareAUC_Test_50, ModelCompareAUC_Test_100)

################################################################################################
################################################################################################

## Work with Agnostic bootobject 200 vs mechanistic

AUCs_XGB_Agnostic_100 <- bootobjectAgnostic_100$t
colnames(AUCs_XGB_Agnostic_100) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_100 <- AUCs_XGB_Agnostic_100[, "AUC_Train"] - AUCs_XGB_Agnostic_100[, "AUC_Test"]
range(DiffAgnostic_100)
quantile(DiffAgnostic_100, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_100 <- data.frame(AUC = AUCs_XGB_Agnostic_100[, "AUC_Train"])

AgnosticAUC_Train_100$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_100 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_100)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_100 <- data.frame(AUC = AUCs_XGB_Agnostic_100[, "AUC_Test"])

AgnosticAUC_Test_100$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_100 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_100)

## add info about the data (training / testing) and the number of pairs
ModelCompareAUC_Train_100$data_type <- "Training"
ModelCompareAUC_Test_100$data_type <- "Testing"

ModelCompareAUC_Train_100$NofFeatAgn <- "100_Pairs"
ModelCompareAUC_Test_100$NofFeatAgn <- "100_Pairs"

####################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/XGB/ModelCompareAUC_200_new.rda")

ModelCompareAUC_Train_100_200 <- rbind(ModelCompareAUC_Train_100, ModelCompareAUC_Train_200)
ModelCompareAUC_Test_100_200 <- rbind(ModelCompareAUC_Test_100, ModelCompareAUC_Test_200)

################################################################################################
################################################################################################

## Work with Agnostic bootobject 250 vs mechanistic

AUCs_XGB_Agnostic_250 <- bootobjectAgnostic_250$t
colnames(AUCs_XGB_Agnostic_250) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

# Calculate the difference and the CI of the difference in the training data
DiffAgnostic_250 <- AUCs_XGB_Agnostic_250[, "AUC_Train"] - AUCs_XGB_Agnostic_250[, "AUC_Test"]
range(DiffAgnostic_250)
quantile(DiffAgnostic_250, c(0.025, 0.975))
#colnames(Diff) <- "Diff"


## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_250 <- data.frame(AUC = AUCs_XGB_Agnostic_250[, "AUC_Train"])

AgnosticAUC_Train_250$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_250 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_250)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_250 <- data.frame(AUC = AUCs_XGB_Agnostic_250[, "AUC_Test"])

AgnosticAUC_Test_250$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_250 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_250)

## add info about the data (training / testing) and the number of pairs
ModelCompareAUC_Train_250$data_type <- "Training"
ModelCompareAUC_Test_250$data_type <- "Testing"

ModelCompareAUC_Train_250$NofFeatAgn <- "250_Pairs"
ModelCompareAUC_Test_250$NofFeatAgn <- "250_Pairs"

####################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/XGB/ModelCompareAUC_500_new.rda")

ModelCompareAUC_Train_250_500 <- rbind(ModelCompareAUC_Train_250, ModelCompareAUC_Train_500)
ModelCompareAUC_Test_250_500 <- rbind(ModelCompareAUC_Test_250, ModelCompareAUC_Test_500)

####################################################################################
####################################################################################
####################################################################################
# Combine all together in one dataframe

ModelCompare_XGB_DiffNoFeat <- rbind(ModelCompareAUC_Train_25_50,
                                     ModelCompareAUC_Test_25_50,
                                     ModelCompareAUC_Train_50_100,
                                     ModelCompareAUC_Test_50_100,
                                     ModelCompareAUC_Train_100_200,
                                     ModelCompareAUC_Test_100_200,
                                     ModelCompareAUC_Train_250_500,
                                     ModelCompareAUC_Test_250_500
)

save(ModelCompare_XGB_DiffNoFeat, file = "./Objs/XGB/ModelCompare_XGB_DiffNoFeat_new.rda")

####################################################################################
####################################################################################
####################################################################################

