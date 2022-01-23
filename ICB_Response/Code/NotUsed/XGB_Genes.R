#################################################################################
## Mohamed Omar
## 07/08/2019
## Goal: Creating a predictive model for bladder cancer progression using XGB
### Agnostic (Using all genes)
##############################################################################


rm(list = ls())

setwd("/Volumes/Macintosh/Research/Projects/ICB_Response")

############################################################################

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

## Load the Processed datasets
load("./Objs/ProcessedDatasets.rda")

RegImmEffec <- read.delim("./RegImmEffec.txt")
RegImmEffec <- as.character(RegImmEffec[-1, ])
RegImmEffec <- gsub("-.+", "", RegImmEffec)


BcellAct <- read.delim("./BcellAct.txt")
BcellAct <- as.character(BcellAct[-1,])
BcellAct <- gsub("-.+","",BcellAct)

WNT_B_Cat <- read.delim("./WNT_B_Cat.txt")
WNT_B_Cat <- as.character(WNT_B_Cat[-1,])
WNT_B_Cat <- gsub("-.+", "", WNT_B_Cat)

AllImmune <- c(RegImmEffec, BcellAct, WNT_B_Cat)
AllImmune <- AllImmune[!duplicated(AllImmune)]

keepGns <- intersect(AllImmune, rownames(Expr_NP))

### Quantile normalize
usedTrainMat <- Expr_NP[keepGns, ]
usedTestMat <- Melanoma_Mat[keepGns, ]

### Associated groups
usedTrainGroup <- Pheno_NP$SR
usedTestGroup <- Melanoma_Group

### Transpose usedTrainMat (making samples as rows instead of columns)
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
names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

names(usedValGroup) <- rownames(Validation)
all(rownames(Validation) == names(usedValGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training1)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train <- cbind(Training, usedTrainGroup1)

## The same for validation
Validation <- as.data.frame(Validation)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_val <- cbind(Validation, usedValGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
names_Test <- c(as.vector(rownames(usedTestMat)))
colnames(Testing) <- names_Test
names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train$usedTrainGroup1)  
levels(Data_train$usedTrainGroup1) <- c(0,1) 
Data_train$usedTrainGroup1 <- factor(Data_train$usedTrainGroup1, levels = c(0,1)) # 0=No,1= Yes 
Train_label <- Data_train$usedTrainGroup1
Train_label <- as.vector(Train_label)
table(Train_label)

## The same for validation
table(Data_val$usedValGroup)  
levels(Data_val$usedValGroup) <- c(0,1) 
Data_val$usedValGroup <- factor(Data_val$usedValGroup, levels = c(0,1)) # 0=No,1= Yes 
Val_label <- Data_val$usedValGroup
Val_label <- as.vector(Val_label)
table(Val_label)


## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_test$usedTestGroup)  
levels(Data_test$usedTestGroup) <- c(0,1)
Data_test$usedTestGroup <- factor(Data_test$usedTestGroup, levels = c(0,1)) #0=No, 1=Yes
Test_label <- Data_test$usedTestGroup
Test_label <- as.vector(Test_label)
table(Test_label)



## Convert to xgb.DMatrix
DataTrain <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataVal <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
DataTest <- xgb.DMatrix(as.matrix(Testing), label = Test_label)

## Creat a watch list
watchlist <- list(train  = DataTrain, test = DataVal)


##########################
# Scale weight (to compensate for un-balanced class sizes)
yes <- sum(Train_label == 1)
no <- sum(Train_label == 0)



## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.001,           #0.3   # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 2,             # 3
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,             #1    # default = 1,   range: (0,1]
  colsample_bytree   = 1,             #1    # default = 1,   range: (0,1]
  colsample_bylevel  = 1,               #1  # default = 1,   range: (0,1]
  lambda             = 1,             # 0  # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

## Make the final model  1411  507
xgb.agnostic <- xgb.train(parameters, DataTrain, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = yes/no)

################################################
#################################################
## ROC curve and CM in the training data
XGB_prob_Train_agnostic <- predict(xgb.agnostic, DataTrain)
ROC_Train_agnostic <- roc(Train_label, XGB_prob_Train_agnostic, plot = FALSE, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the training cohort (Mechanistic)")
ROC_Train_agnostic

thr_Train <- coords(roc(Train_label, XGB_prob_Train_agnostic, levels = c("0", "1"), direction = "<"), transpose = TRUE, "best")["threshold"]
thr_Train

## Convert predicted probabilities to binary outcome
prediction_Train_agnostic <- as.numeric(XGB_prob_Train_agnostic > thr_Train)
print(head(prediction_Train_agnostic))

Train_label <- factor(Train_label, levels = c(0,1))
prediction_Train_agnostic <- factor(prediction_Train_agnostic, levels = c(0,1))

# Confusion matrix in the training data
CM_Train <- confusionMatrix(prediction_Train_agnostic, Train_label, positive = "1")
CM_Train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = prediction_Train_agnostic, actuals = Train_label)
MCC_Train

#####################################
######################################

## Predict in the Test data
xgb_prob_test_agnostic <- predict(xgb.agnostic, DataTest)

#########

## Convert predicted probabilities to binary outcome
prediction_Test_agnostic <- as.numeric(xgb_prob_test_agnostic > thr_Train)
print(head(prediction_Test_agnostic))

Test_label <- factor(Test_label, levels = c(0,1))
prediction_Test_agnostic <- factor(prediction_Test_agnostic, levels = c(0,1))

## Confusion matrix
CM_Test <- confusionMatrix(prediction_Test_agnostic, Test_label, positive = "1")
CM_Test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = prediction_Test_agnostic, actuals = Test_label)
MCC_Test

## ROC curve and AUC
roc(Test_label, xgb_prob_test_agnostic, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the testing cohort (Agnostic-Genes)")


##############################################
