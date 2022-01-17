########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: SVM (Poly) for predicting bladder cancer progression
## Agnostic (Using all gene pairs)

###########################################################################

## Clean the work environment
rm(list = ls())

## Set the working directory
#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(mltools)

#######################################################################

## Load data
load("./Objs/KTSP/KTSP_STATs_Agnostic_Combined_500.rda")
load("./Objs/MetastasisDataGood.rda")


### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Agnostic_500)

## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Agnostic_500)

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
# 200: 3 / 0.1 / 0.25
# 500: 2 / 0.01 / 0.25
set.seed(333)
Grid <- expand.grid(degree = 2, scale = 0.01, C = 0.25)
fit.svmPoly_agnostic_OnKTSP_500 <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
fit.svmPoly_agnostic_OnKTSP_500

##########################################################
## Predict in the training data
train_pred_classes_svmPoly_agnostic_OnKTSP_500 <- predict(fit.svmPoly_agnostic_OnKTSP_500, Training, type="raw")
table(train_pred_classes_svmPoly_agnostic_OnKTSP_500)

##
## ROC stat for the training data
train_pred_prob_svmPoly_agnostic_OnKTSP_500 <- predict(fit.svmPoly_agnostic_OnKTSP_500, Training, type = "prob")
roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP_500[,2], plot = F, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC train SVM Poly (Agnostic)")

## Creat a confusion matrix (the predictions against the true classes)
Confusion_train_svmPoly_OnKTSP <- confusionMatrix(train_pred_classes_svmPoly_agnostic_OnKTSP_500, usedTrainGroup, positive = "Mets")
Confusion_train_svmPoly_OnKTSP

MCC_Train <- mltools::mcc(pred = train_pred_classes_svmPoly_agnostic_OnKTSP_500, actuals = usedTrainGroup)
MCC_Train

#########################################################
## Predict in the testing data
test_pred_classes_svmPoly_agnostic_OnKTSP_500 <- predict(fit.svmPoly_agnostic_OnKTSP_500, Testing, type="raw")
table(test_pred_classes_svmPoly_agnostic_OnKTSP_500)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly_OnKTSP <- confusionMatrix(test_pred_classes_svmPoly_agnostic_OnKTSP_500, usedTestGroup, positive = "Mets")
Confusion_test_svmPoly_OnKTSP

MCC_Test <- mltools::mcc(pred = test_pred_classes_svmPoly_agnostic_OnKTSP_500, actuals = usedTestGroup)
MCC_Test

## ROC/AUC in the Testing set
test_pred_prob_svmPoly_agnostic_OnKTSP_500 <- predict(fit.svmPoly_agnostic_OnKTSP_500, Testing, type = "prob")
roc(usedTestGroup, test_pred_prob_svmPoly_agnostic_OnKTSP_500[,2], plot = F, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", ci = T,lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Agnostic)")

## Top 100 predictors
Importance_SVMPoly <- varImp(fit.svmPoly_agnostic_OnKTSP_500, scale = FALSE)
png("./Figs/SVM/SVM_varImp_Agnostic_OnKTSP.png", width = 2000, height = 2000, res = 200)
plot(Importance_SVMPoly, top = 20, main = "Agnostic SVM top 20 features")
dev.off()

Importance_SVMPoly <- Importance_SVMPoly$importance
Importance_SVMPoly <- Importance_SVMPoly[order(Importance_SVMPoly$Progression, decreasing = TRUE),]
genes_SVMPoly <- rownames(Importance_SVMPoly)[1:20]
save(Importance_SVMPoly, file = "./Objs/SVM/genes_SVMPloy_Agnostic.rdata")


########################

####
# Run this one with 200 and 500 then run the agnostic (standard) and the mechanistic then save them

save(train_pred_prob_svmPoly_mechanistic_OnKTSP, test_pred_prob_svmPoly_mechanistic_OnKTSP, 
     train_pred_prob_svmPoly_agnostic_OnKTSP, test_pred_prob_svmPoly_agnostic_OnKTSP,
     train_pred_prob_svmPoly_agnostic_OnKTSP_200, test_pred_prob_svmPoly_agnostic_OnKTSP_200,
     train_pred_prob_svmPoly_agnostic_OnKTSP_500, test_pred_prob_svmPoly_agnostic_OnKTSP_500, 
     file = "./Objs/SVM/SVMAUCStats.rda")


