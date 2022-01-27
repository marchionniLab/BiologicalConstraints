## Clean the work environment
rm(list = ls())


## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(mltools)
library(doParallel)
#######################################################################

cl <- makePSOCKcluster(12)
registerDoParallel(cl)

#######################################################################

## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_GSE115821Out.rda")
load("./Objs/icbData_GSE115821Out.rda")



### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Mechanistic)


## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Mechanistic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

######
control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)

###########################################################################
###########################################################################

## Model: SVM Poly

## 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=control, tuneLength = 5, metric = "ROC")
fit.svmPoly

## Training using all data (using the best parameters)
Grid <- expand.grid(degree = 1, scale = 0.01, C = 0.25)
set.seed(333)
fit.svmPoly_mechanistic_OnKTSP <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = c("ROC"))
fit.svmPoly_mechanistic_OnKTSP

##########################################################
## Predict in the training data

## ROC stat for the training data
train_pred_prob_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Training, type = "prob")
ROCTrain <- roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2], plot = F, print.auc=TRUE, ci = T, levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Mechanistic)")
ROCTrain

# The best threshold
thr <- coords(roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2], levels = c("NR", "R"), direction = "<"), "best")["threshold"]
thr <- thr$threshold
thr

# Predict classes according to the best threshold
train_pred_classes_svmPoly_mechanistic_OnKTSP <- ifelse(train_pred_prob_svmPoly_mechanistic_OnKTSP[,2] > thr, "R", "NR")
table(train_pred_classes_svmPoly_mechanistic_OnKTSP)
# Convert to factor
train_pred_classes_svmPoly_mechanistic_OnKTSP <- factor(train_pred_classes_svmPoly_mechanistic_OnKTSP, levels = c("NR", "R"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_train_svmPoly_OnKTSP <- confusionMatrix(train_pred_classes_svmPoly_mechanistic_OnKTSP, usedTrainGroup, positive = "R")
Confusion_train_svmPoly_OnKTSP

MCC_Train <- mltools::mcc(pred = train_pred_classes_svmPoly_mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, Confusion_train_svmPoly_OnKTSP$overall["Accuracy"], Confusion_train_svmPoly_OnKTSP$byClass["Balanced Accuracy"], Confusion_train_svmPoly_OnKTSP$byClass["Sensitivity"], Confusion_train_svmPoly_OnKTSP$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#############################################################
## Predict in the testing data 

## ROC/AUC in the Testing set
test_pred_prob_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Testing, type = "prob")
ROCTest <- roc(usedTestGroup, test_pred_prob_svmPoly_mechanistic_OnKTSP[,2], plot = F, print.auc=TRUE, ci = T, levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Mechanistic)")
ROCTest

# Predict classes according to the best threshold
test_pred_classes_svmPoly_mechanistic_OnKTSP <- ifelse(test_pred_prob_svmPoly_mechanistic_OnKTSP[,2] > thr, "R", "NR")
table(test_pred_classes_svmPoly_mechanistic_OnKTSP)
# Convert to factor
test_pred_classes_svmPoly_mechanistic_OnKTSP <- factor(test_pred_classes_svmPoly_mechanistic_OnKTSP, levels = c("NR", "R"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly_OnKTSP <- confusionMatrix(test_pred_classes_svmPoly_mechanistic_OnKTSP, usedTestGroup, positive = "R")
Confusion_test_svmPoly_OnKTSP

MCC_Mechanistic <- mltools::mcc(pred = test_pred_classes_svmPoly_mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Mechanistic

# Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, Confusion_test_svmPoly_OnKTSP$overall["Accuracy"], Confusion_test_svmPoly_OnKTSP$byClass["Balanced Accuracy"], Confusion_test_svmPoly_OnKTSP$byClass["Sensitivity"], Confusion_test_svmPoly_OnKTSP$byClass["Specificity"], MCC_Mechanistic))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
GSE115821_Out_SVM_MechPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(GSE115821_Out_SVM_MechPerformance, file = "./Objs/SVM/GSE115821_Out_SVM_MechPerformance.rda")

stopCluster(cl)
