## Clean the work environment
rm(list = ls())

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(mltools)
library(switchBox)
library(doParallel)

#######################################################################

cl <- makePSOCKcluster(12)
registerDoParallel(cl)

#######################################################################

## Load data
load("./Objs/icbData_GSE91061Out_Pre.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")
usedTestMat <- testMat

######
# remove genes with NA values from the training and testing matrices
sel <- which(apply(usedTrainMat, 1, function(x) all(is.finite(x)) ))
usedTrainMat <- usedTrainMat[sel, ] 

sel2 <- which(apply(usedTestMat, 1, function(x) all(is.finite(x)) ))
usedTestMat <- usedTestMat[sel2, ] 

keep <- intersect(rownames(usedTrainMat), rownames(usedTestMat))
usedTrainMat <- usedTrainMat[keep, ]
usedTestMat <- usedTestMat[keep,]

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#################################################################
### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(usedTrainMat)
#names_train <- c(as.vector(rownames(usedTrainMat)))
#colnames(Training) <- names_train

## Making sure that sample names are identical in both Training and usedTrainGroup
#names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(usedTestMat)
#names_Test <- c(as.vector(rownames(usedTestMat)))
#colnames(Testing) <- names_Test
#names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))


######
control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = F)

###########################################################################
###########################################################################

## Model: SVM Poly

## 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=control, metric = c("ROC"))
fit.svmPoly

## Training using all data (using the best parameters)
Grid <- expand.grid(degree = 1, scale = 0.001, C = 1)
set.seed(333)
fit.svmPoly_agnostic<- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = "ROC")
fit.svmPoly_agnostic

##########################################################
## Predict in the training data

## ROC stat for the training data
train_pred_prob_svmPoly_agnostic <- predict(fit.svmPoly_agnostic, Training, type = "prob")
ROCTrain <- roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic[,2], plot = F, ci = T, print.auc=TRUE, levels = c("NR", "R"), col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Agnostic)")
ROCTrain

# Best threshold
thr <- coords(roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic[,1], levels = c("NR", "R"), direction = "<"), "best")["threshold"]
thr <- thr$threshold
thr

## Convert predicted probabilities to binary outcome
train_pred_classes_svmPoly_agnostic <- ifelse(train_pred_prob_svmPoly_agnostic[,1] >= thr, "R", "NR")
table(train_pred_classes_svmPoly_agnostic)
# Convert to factor
train_pred_classes_svmPoly_agnostic <- factor(train_pred_classes_svmPoly_agnostic, levels = c("NR", "R"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_train_svmPoly <- confusionMatrix(train_pred_classes_svmPoly_agnostic, usedTrainGroup, positive = "R")
Confusion_train_svmPoly

MCC_Train <- mltools::mcc(pred = train_pred_classes_svmPoly_agnostic, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, Confusion_train_svmPoly$overall["Accuracy"], Confusion_train_svmPoly$byClass["Balanced Accuracy"], Confusion_train_svmPoly$byClass["Sensitivity"], Confusion_train_svmPoly$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#########################################################
## Predict in the testing data
## ROC/AUC in the Testing set
test_pred_prob_svmPoly_agnostic <- predict(fit.svmPoly_agnostic, Testing, type = "prob")
ROCTest <- roc(usedTestGroup, test_pred_prob_svmPoly_agnostic[,2], plot = F, print.auc=TRUE, ci = T, levels = c("NR", "R"), col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Agnostic)")
ROCTest

# Predict classes according to the best threshold
test_pred_classes_svmPoly_agnostic <- ifelse(test_pred_prob_svmPoly_agnostic[,1] > thr, "R", "NR")
table(test_pred_classes_svmPoly_agnostic)
# Convert to factor
test_pred_classes_svmPoly_agnostic <- factor(test_pred_classes_svmPoly_agnostic, levels = c("NR", "R"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly <- confusionMatrix(test_pred_classes_svmPoly_agnostic, usedTestGroup, positive = "R")
Confusion_test_svmPoly

MCC_Test <- mltools::mcc(pred = test_pred_classes_svmPoly_agnostic, actuals = usedTestGroup)
MCC_Test

## Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, Confusion_test_svmPoly$overall["Accuracy"], Confusion_test_svmPoly$byClass["Balanced Accuracy"], Confusion_test_svmPoly$byClass["Sensitivity"], Confusion_test_svmPoly$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
GSE91061_Out_SVM_IndvGenes_AgnosticPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(GSE91061_Out_SVM_IndvGenes_AgnosticPerformance, file = "./Objs/SVM/GSE91061_Out_SVM_IndvGenes_AgnosticPerformance.rda")

stopCluster(cl)

########################
