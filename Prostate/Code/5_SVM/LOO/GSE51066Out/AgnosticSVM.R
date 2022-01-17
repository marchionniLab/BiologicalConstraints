########################################################################### 
## Mohamed Omar
## 27/11/2019
## Goal: SVM (Poly) for predicting prostate cancer metastasis
## Agnostic (Using all gene pairs)
## Cross study validation : GSE51066 out
###########################################################################

## Clean the work environment
rm(list = ls())

## Set the working directory
setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(mltools)

#######################################################################

## Load data
load("./Objs/KTSP/LOO/KTSP_STATs_Agnostic_GSE51066Out.rda")
load("./Objs/LOO/MetastasisData_GSE51066Out.rda")


### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

#################################################################
## Oversampling of the training data set to compensate for the un-balanced classes
# set.seed(333)
# Data_train <- as.data.frame(Data_train)
# Data_train[,332] <- as.factor(Data_train[,332])
# Over_Train <- SMOTE(usedTrainGroup~., data = Data_train, perc.over = 300, perc.under = 134)
# table(Over_Train[,332])

######
control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

###########################################################################
###########################################################################

## Model: SVM Poly

## 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=control, metric = c("ROC"))
fit.svmPoly

## Training using all data (using the best parameters)
Grid <- expand.grid(degree = 3, scale = 0.01, C = 0.25)
set.seed(333)
fit.svmPoly_agnostic_OnKTSP <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
fit.svmPoly_agnostic_OnKTSP

##########################################################
## Predict in the training data

## ROC stat for the training data
train_pred_prob_svmPoly_agnostic_OnKTSP <- predict(fit.svmPoly_agnostic_OnKTSP, Training, type = "prob")
roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP[,2], plot = F, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Agnostic)")

# The best threshold
thr <- coords(roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr <- thr$threshold
thr

# Predict classes according to the best threshold
train_pred_classes_svmPoly_agnostic_OnKTSP <- ifelse(train_pred_prob_svmPoly_agnostic_OnKTSP[,2] > thr, "Mets", "No_Mets")
table(train_pred_classes_svmPoly_agnostic_OnKTSP)
# Convert to factor
train_pred_classes_svmPoly_agnostic_OnKTSP <- factor(train_pred_classes_svmPoly_agnostic_OnKTSP, levels = c("No_Mets", "Mets"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_train_svmPoly_OnKTSP <- confusionMatrix(train_pred_classes_svmPoly_agnostic_OnKTSP, usedTrainGroup, positive = "Mets")
Confusion_train_svmPoly_OnKTSP

MCC_Train <- mltools::mcc(pred = train_pred_classes_svmPoly_agnostic_OnKTSP, actuals = usedTrainGroup)
MCC_Train

#########################################################
## Predict in the testing data

## ROC/AUC in the Testing set
test_pred_prob_svmPoly_agnostic_OnKTSP <- predict(fit.svmPoly_agnostic_OnKTSP, Testing, type = "prob")
roc(usedTestGroup, test_pred_prob_svmPoly_agnostic_OnKTSP[,2], plot = F, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Agnostic)")

# Predict classes according to the best threshold
test_pred_classes_svmPoly_agnostic_OnKTSP <- ifelse(test_pred_prob_svmPoly_agnostic_OnKTSP[,2] > thr, "Mets", "No_Mets")
table(test_pred_classes_svmPoly_agnostic_OnKTSP)
# Convert to factor
test_pred_classes_svmPoly_agnostic_OnKTSP <- factor(test_pred_classes_svmPoly_agnostic_OnKTSP, levels = c("No_Mets", "Mets"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly_OnKTSP <- confusionMatrix(test_pred_classes_svmPoly_agnostic_OnKTSP, usedTestGroup, positive = "Mets")
Confusion_test_svmPoly_OnKTSP

MCC_Test <- mltools::mcc(pred = test_pred_classes_svmPoly_agnostic_OnKTSP, actuals = usedTestGroup)
MCC_Test


## Top predictors
#Importance_SVMPoly <- varImp(fit.svmPoly_agnostic_OnKTSP, scale = TRUE)
# png("./Figs/SVM/SVM_varImp.png", width = 2000, height = 2000, res = 200)
# plot(Importance_SVMPoly, top = 20, main = "Agnostic SVM top 20 genes")
# dev.off()
# 
#Importance_SVMPoly <- Importance_SVMPoly$importance
#Importance_SVMPoly <- Importance_SVMPoly[order(Importance_SVMPoly$Mets, decreasing = TRUE),]
#Importance_SVMPoly <- Importance_SVMPoly[!(Importance_SVMPoly$Mets == 0), ]
# genes_SVMPoly <- rownames(Importance_SVMPoly)[1:20]
# save(Importance_SVMPoly, file = "./Objs/SVM/genes_SVMPloy_Agnostic.rdata")

########################

########################
