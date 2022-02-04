rm(list = ls())

## Load necessary packages
library(xgboost)
library(preprocessCore)
library(limma)
library(pROC)
library(caret)
library(ggplot2)
library(dplyr)
library(mltools)
library(reshape2)
library(RColorBrewer)
library(ggpubr)
library(plotROC)

## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_VanAllenOut.rda")
load("./Objs/icbData_VanAllenOut_Pre.rda")

Training <- t(KTSP_STATs_Train_Mechanistic)


usedTrainGroup <- trainGroup
usedTestGroup <- testGroup


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
Data_train <- cbind(Training, usedTrainGroup1)

## The same for validation
Validation <- as.data.frame(Validation)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_val <- cbind(Validation, usedValGroup)


########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Mechanistic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train$usedTrainGroup1)  
levels(Data_train$usedTrainGroup1) <- c(0,1) 
Data_train$usedTrainGroup <- factor(Data_train$usedTrainGroup1, levels = c(0,1)) # 0=No,1= Yes 
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
R <- sum(Train_label == 1)
NR <- sum(Train_label == 0)

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",        
  silent             = 1,                 
  # Booster Parameters
  eta                = 0.001,           
  gamma              = 0,           
  max_depth          = 2,           
  min_child_weight   = 1,          
  subsample          = 0.5,          
  colsample_bytree   = 1,         
  colsample_bylevel  = 1,    
  lambda             = 1,      
  alpha              = 0,           
  # Task Parameters
  objective          = "binary:logistic", 
  eval_metric        = "auc"
)

## Make the final model
xgb.mechanistic_OnKTSP <- xgb.train(parameters, DataTrain, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = R/NR)

################################################
#################################################
## ROC curve and CM in the training data
XGB_prob_Train_mechanistic_OnKTSP <- predict(xgb.mechanistic_OnKTSP, DataTrain)
ROC_Train_mechanistic_OnKTSP <- roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, plot = FALSE, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the training cohort (Mechanistic)")
ROC_Train_mechanistic_OnKTSP

thr_Train <- coords(roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, levels = c("0", "1"), direction = "<"), transpose = TRUE, "best")["threshold"]
thr_Train

## Convert predicted probabilities to binary outcome
prediction_Train_mechanistic_OnKTSP <- as.numeric(XGB_prob_Train_mechanistic_OnKTSP > thr_Train)
print(head(prediction_Train_mechanistic_OnKTSP))

Train_label <- factor(Train_label, levels = c(0,1))
prediction_Train_mechanistic_OnKTSP <- factor(prediction_Train_mechanistic_OnKTSP, levels = c(0,1))

# Confusion matrix in the training data
CM_Train <- confusionMatrix(prediction_Train_mechanistic_OnKTSP, Train_label, positive = "1")
CM_Train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = prediction_Train_mechanistic_OnKTSP, actuals = Train_label)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROC_Train_mechanistic_OnKTSP$ci, CM_Train$overall["Accuracy"], CM_Train$byClass["Balanced Accuracy"], CM_Train$byClass["Sensitivity"], CM_Train$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#####################################
######################################

## Predict in the Test data
xgb_prob_test_mechanistic_OnKTSP <- predict(xgb.mechanistic_OnKTSP, DataTest)

## Convert predicted probabilities to binary outcome
prediction_Test_mechanistic_OnKTSP <- as.numeric(xgb_prob_test_mechanistic_OnKTSP > thr_Train)
print(head(prediction_Test_mechanistic_OnKTSP))

Test_label <- factor(Test_label, levels = c(0,1))
prediction_Test_mechanistic_OnKTSP <- factor(prediction_Test_mechanistic_OnKTSP, levels = c(0,1))

## Confusion matrix
CM_Test <- confusionMatrix(prediction_Test_mechanistic_OnKTSP, Test_label, positive = "1")
CM_Test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = prediction_Test_mechanistic_OnKTSP, actuals = Test_label)
MCC_Test

## ROC curve and AUC
ROCTest <- roc(Test_label, xgb_prob_test_mechanistic_OnKTSP, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the testing cohort (Mechanistic)")
ROCTest

# Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, CM_Test$overall["Accuracy"], CM_Test$byClass["Balanced Accuracy"], CM_Test$byClass["Sensitivity"], CM_Test$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
VanAllen_Out_XGB_MechPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(VanAllen_Out_XGB_MechPerformance, file = "./Objs/XGB/VanAllen_Out_XGB_MechPerformance.rda")
