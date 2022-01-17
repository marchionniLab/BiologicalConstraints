



rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

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


load("./Objs/KTSP_STATs_Agnostic.rda")
load("./Objs/MetastasisData_New.rda")

usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

Training <- t(KTSP_STATs_Train_Agnostic)




## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))

## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train$usedTrainGroup)  
levels(Data_train$usedTrainGroup) <- c(1,0) 
Data_train$usedTrainGroup <- factor(Data_train$usedTrainGroup, levels = c(0,1)) # 0=No,1= Yes 
Train_label <- Data_train$usedTrainGroup
Train_label <- as.vector(Train_label)
table(Train_label)


## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_test$usedTestGroup)  
levels(Data_test$usedTestGroup) <- c(1,0)
Data_test$usedTestGroup <- factor(Data_test$usedTestGroup, levels = c(0,1)) #0=No, 1=Yes
Test_label <- Data_test$usedTestGroup
Test_label <- as.vector(Test_label)
table(Test_label)



## Convert to xgb.DMatrix
DataTrain <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataTest <- xgb.DMatrix(as.matrix(Testing), label = Test_label)

## Creat a watch list
watchlist <- list(train  = DataTrain, test = DataTest)

##########################
# Scale weight (to compensate for un-balanced class sizes)
Mets <- sum(Train_label == 1)
No_Mets <- sum(Train_label == 0)

## Make a 5-fold CV model to determine the best number of trees/iterations
# hyper_grid <- expand.grid(
#   eta = c(0.1,0.3),
#   max_depth = c(1,3,6),
#   min_child_weight = 1,
#   subsample = c(0.7,0.8,0.9,1),
#   colsample_bytree = c(0.7,0.8,0.9,1),
#   gamma = c(0,1),
#   lambda = c(0,0.5,1),
#   alpha = 0,
#   optimal_trees = 0,               # a place to dump results
#   max_AUC = 0                     # a place to dump results
# )
# 
# 
# ##########################
# # Scale weight (to compensate for un-balanced class sizes)
# Progression <- sum(Train_label == 1)
# NoProgression <- sum(Train_label == 0)
# 
# 
# # grid search
# for(i in 1:nrow(hyper_grid)) {
# 
#   # create parameter list
#   params <- list(
#     eta = hyper_grid$eta[i],
#     max_depth = hyper_grid$max_depth[i],
#     min_child_weight = hyper_grid$min_child_weight[i],
#     subsample = hyper_grid$subsample[i],
#     colsample_bytree = hyper_grid$colsample_bytree[i],
#     gamma = hyper_grid$gamma[i],
#     lambda = hyper_grid$lambda[i],
#     alpha = hyper_grid$alpha[i]
#   )
# 
#   # reproducibility
#   set.seed(333)
# 
#   # train model
#   xgb.tune <- xgb.cv(
#     params = params,
#     data = DataTrain,
#     nrounds = 500,
#     nfold = 5,
#     objective = "binary:logistic",
#     eval_metric        = "auc",# for regression models
#     verbose = 1,               # silent,
#     early_stopping_rounds = 50,
#     scale_pos_weight = NoProgression/Progression
#   )
# 
#   # add max training AUC and trees to grid
#   hyper_grid$optimal_trees[i] <- which.max(xgb.tune$evaluation_log$test_auc_mean)
#   hyper_grid$max_AUC[i] <- max(xgb.tune$evaluation_log$test_auc_mean)
# }
# 
# hyper_grid %>% arrange(max_AUC)
# View(hyper_grid)
# save(hyper_grid, file = "./Objs/XGB/XGBAgnosticPairs_HyperGrid.rda")
# Best: rownumber 203

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.3    # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.5,           #0.5      # default = 1,   range: (0,1]
  colsample_bytree   = 0.7,              #1   # default = 1,   range: (0,1]
  colsample_bylevel  = 1,               #1   # default = 1,   range: (0,1]
  lambda             = 0,             # 0 # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc",
  seed               = 333               # reproducability seed
)

## Make the final model  1411  507
xgb.agnostic_OnKTSP <- xgb.train(parameters, DataTrain, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = No_Mets/Mets)

################################################
#################################################
## ROC curve and CM in the training data
XGB_prob_Train_agnostic_OnKTSP <- predict(xgb.agnostic_OnKTSP, DataTrain)
ROC_Train_agnostic_OnKTSP <- roc(Train_label, XGB_prob_Train_agnostic_OnKTSP, plot = FALSE, print.auc = TRUE, levels = c("0", "1"), col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE, main = "XGB ROC curve in the training cohort (Mechanistic)")

thr_Train <- coords(roc(Train_label, XGB_prob_Train_agnostic_OnKTSP, levels = c("0", "1"), ), transpose = TRUE, "best")["threshold"]
thr_Train

## Convert predicted probabilities to binary outcome
prediction_Train_agnostic_OnKTSP <- as.numeric(XGB_prob_Train_agnostic_OnKTSP > thr_Train)
print(head(prediction_Train_agnostic_OnKTSP))

Train_label <- factor(Train_label, levels = c(0,1))
prediction_Train_agnostic_OnKTSP <- factor(prediction_Train_agnostic_OnKTSP, levels = c(0,1))

# Confusion matrix in the training data
CM_Train <- confusionMatrix(prediction_Train_agnostic_OnKTSP, Train_label, positive = "1")
CM_Train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = prediction_Train_agnostic_OnKTSP, actuals = Train_label)
MCC_Train

#####################################
######################################

## Predict in the Test data
xgb_prob_test_agnostic_OnKTSP <- predict(xgb.agnostic_OnKTSP, DataTest)

### Threshold
thr <- coords(roc(Test_label, xgb_prob_test_agnostic_OnKTSP, levels = c("0", "1"), ), transpose = TRUE, "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(Test_label, xgb_prob_test_agnostic_OnKTSP, levels = c("0", "1"), ), transpose = TRUE ,"local maximas")


#############################################

## Convert predicted probabilities to binary outcome
prediction_Test_agnostic_OnKTSP <- as.numeric(xgb_prob_test_agnostic_OnKTSP > thr)
print(head(prediction_Test_agnostic_OnKTSP))

Test_label <- factor(Test_label, levels = c(0,1))
prediction_Test_agnostic_OnKTSP <- factor(prediction_Test_agnostic_OnKTSP, levels = c(0,1))

## Confusion matrix
CM_Test <- confusionMatrix(prediction_Test_agnostic_OnKTSP, Test_label, positive = "1")
CM_Test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = prediction_Test_agnostic_OnKTSP, actuals = Test_label)
MCC_Test

## ROC curve and AUC
roc(Test_label, xgb_prob_test_agnostic_OnKTSP, plot = F, print.auc = TRUE, levels = c("0", "1"), col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE, main = "XGB ROC curve in the testing cohort (Agnostic)")


##############################################
## Classification error rate
err <- mean(as.numeric(xgb_prob_test_agnostic_OnKTSP > thr) != Test_label)
err


##############

# create importance matrix
importance_matrix <- xgb.importance(model = xgb.agnostic_OnKTSP)
write.csv(importance_matrix, file = "./Objs/XGB/XGB_importanceMatrix_AgnosticPairs.csv")


# variable importance plot
importance_matrix <- importance_matrix[order(importance_matrix$Gain, decreasing = TRUE), ]

png(filename = "./Figs/XGB/XGB_Importance_AgnosticPairs.png", width = 3000, height = 2000, res = 300)
ggplot(data=importance_matrix, aes(x = reorder(Feature, -Gain), y = Gain)) +
  geom_bar(stat="identity", colour = "black", fill = "lightgray") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "XGB Feature Importance (Top 30) (Agnostic Pairs)", x = "Features", y = "Information Gain")
dev.off()

####################
####################

