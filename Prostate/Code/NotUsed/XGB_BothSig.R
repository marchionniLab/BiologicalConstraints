#################################################################################
## Mohamed Omar
## 07/08/2019
## Goal: Creating a predictive model for bladder cancer progression using XGB
### Agnostic (Using all genes)
##############################################################################


rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

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

## Load data
load("./Objs/MetastasisData_New.rda")

## Load TF_MiR genes
#Genes <- load("/Users/mohamedomar/Desktop/Genes/ENCODE_TF2Targets.rda")
load("./Objs/NewSigGenes.rda")
load("./Objs/OncoDxGenes.rda")

BothSig <- c(New_SignatureGenes, OncoDx_Genes)

### Common genes
keepGns <- intersect(BothSig, rownames(trainMat))

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")[keepGns, ]
usedTestMat <- normalizeBetweenArrays(testMat, method = "quantile")[keepGns, ]

#usedTrainMat <- t(scale(t(usedTrainMat), center = TRUE, scale = TRUE))
#usedTestMat <- t(scale(t(usedTestMat), center = TRUE, scale = TRUE))

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
names_train <- c(as.vector(rownames(usedTrainMat)))
colnames(Training) <- names_train



## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))

## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
names_Test <- c(as.vector(rownames(usedTestMat)))
colnames(Testing) <- names_Test
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

######################################################
## For second validation
#Testing2 <- t(usedTestMat2)
#names_Test2 <- c(as.vector(rownames(usedTestMat2)))
#colnames(Testing2) <- names_Test2
#names(usedTestGroup2) <- rownames(Testing2)
#all(rownames(Testing2) == names(usedTestGroup2))

###########################################################
## Combine both the Expression matrix and the phenotype into one matrix
#Testing2 <- as.data.frame(Testing2)
#Data_test2 <- cbind(Testing2, usedTestGroup2)

## Converting classes from Progression/NoProgression Format to 0-1 Format
# table(Data_test2$usedTestGroup2)  
# levels(Data_test2$usedTestGroup2) <- c(1,0)
# Data_test2$usedTestGroup2 <- factor(Data_test2$usedTestGroup2, levels = c(0,1)) #0=No, 1=Yes
# Test_label2 <- Data_test2$usedTestGroup2
# Test_label2 <- as.vector(Test_label2)
# table(Test_label2)
##############################

## Convert to xgb.DMatrix
DataTrain <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataTest <- xgb.DMatrix(as.matrix(Testing), label = Test_label)
#DataTest2 <- xgb.DMatrix(as.matrix(Testing2), label = Test_label2)

## Creat a watch list
watchlist <- list(train  = DataTrain, test = DataTest)



## Make a 5-fold CV model to determine the best number of trees/iterations
# hyper_grid <- expand.grid(
#   eta = 0.1,
#   max_depth = c(3,6),
#   min_child_weight = c(1, 3, 6),
#   subsample = 1,
#   colsample_bytree = 1,
#   gamma = 0,
#   lambda = c(0,0.2,0.5, 1),
#   alpha = c(0, 0.2, 0.5, 1),
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
# hyper_grid %>% arrange(max_AUC) %>% head(100)
# View(hyper_grid)
# write.csv(hyper_grid, file = "./Objs/XGB/hyper_grid_Agnostic.csv")

##########################
# Scale weight (to compensate for un-balanced class sizes)
Metasatsis <- sum(Train_label == 1)
NoMetastasis <- sum(Train_label == 0)



## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.3,           #0.3   # default = 0.3, range: [0,1]
  gamma              = 0,             #0   # default = 0,   range: [0,∞]
  max_depth          = 2,             # 2
  min_child_weight   = 1,             #1    # default = 1,   range: [0,∞]
  subsample          = 0.4,             #0.4    # default = 1,   range: (0,1]
  colsample_bytree   = 1,             # 1     # default = 1,   range: (0,1]
  colsample_bylevel  = 0.6,               #0.6  # default = 1,   range: (0,1]
  lambda             = 2,             # 2  # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc",
  seed               = 333               # reproducability seed
)

## Make the final model  1411  507
xgb.agnostic <- xgb.train(parameters, DataTrain, nrounds = 1000, watchlist,  early_stopping_rounds = 200, scale_pos_weight = NoMetastasis/Metasatsis)

################################################
#################################################
## ROC curve and CM in the training data
XGB_prob_Train_agnostic <- predict(xgb.agnostic, DataTrain)
ROC_Train_agnostic <- roc(Train_label, XGB_prob_Train_agnostic, plot = FALSE, print.auc = TRUE, levels = c("0", "1"), col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE, main = "XGB ROC curve in the training cohort (Mechanistic)")

thr_Train <- coords(roc(Train_label, XGB_prob_Train_agnostic, levels = c("0", "1"), ), transpose = TRUE, "best")["threshold"]
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

### Threshold
thr <- coords(roc(Test_label, xgb_prob_test_agnostic, levels = c("0", "1"), ), transpose = TRUE, "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(Test_label, xgb_prob_test_agnostic, levels = c("0", "1"), ), transpose = TRUE ,"local maximas")


#############################################

## Convert predicted probabilities to binary outcome
prediction_Test_agnostic <- as.numeric(xgb_prob_test_agnostic > thr)
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
png(filename = "./Figs/XGB/AgnosticGenes_ROC_Test.png", width = 2000, height = 2000, res = 300)
roc(Test_label, xgb_prob_test_agnostic, plot = TRUE, print.auc = TRUE, levels = c("0", "1"), col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE, main = "XGB ROC curve in the testing cohort (Agnostic-Genes)")
dev.off()


##############
# create importance matrix
importance_matrix <- xgb.importance(model = xgb.agnostic)
#write.csv(importance_matrix, file = "./Objs/XGB/XGB_importanceMatrix_Agnostic.csv")


# variable importance plot
importance_matrix <- importance_matrix[order(importance_matrix$Gain, decreasing = TRUE), ]
ImportantGns <- importance_matrix$Feature
save(ImportantGns, file = "./ImportantGnsXGB.rda")

png(filename = "./Figs/XGB/XGB_Importance_Agnostic.png", width = 3000, height = 2000, res = 300)
ggplot(data=importance_matrix[1:30,], aes(x = reorder(Feature, -Gain), y = Gain)) +
  geom_bar(stat="identity", colour = "black", fill = "lightgray") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "XGB Feature Importance (Top 30) (Agnostic)", x = "Features", y = "Information Gain")
dev.off()

####################
####################

