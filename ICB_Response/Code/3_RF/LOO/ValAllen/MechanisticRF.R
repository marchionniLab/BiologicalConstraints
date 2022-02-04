## Clean the work environment
rm(list = ls())

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(randomForest)
library(varSelRF)


## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_VanAllenOut.rda")
load("./Objs/icbData_VanAllenOut_Pre.rda")


usedTrainGroup <- trainGroup
usedTestGroup <- testGroup


predictor_data <- t(KTSP_STATs_Train_Mechanistic)

tmp <- as.vector(table(usedTrainGroup))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

set.seed(333)
tuneRF(x=predictor_data, y=usedTrainGroup, plot = TRUE, improve = 0.01, ntreeTry = 500, proximity = TRUE, sampsize = sampsizes, na.action = na.omit)

set.seed(333)
RF_Mechanistic_OnKTSP <- randomForest(x =predictor_data, y=usedTrainGroup, importance = TRUE, ntree = 50, mtry = 2, proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
RF_Mechanistic_OnKTSP

# Predictions in the training data
train_pred_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data, type="response")
confusionTrain <- confusionMatrix(train_pred_responses_Mechanistic_OnKTSP, usedTrainGroup, positive = "R")
confusionTrain

# Calculate Matthews correlation coefficient
MCC_Mechanistic_Train <- mcc(preds = train_pred_responses_Mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Mechanistic_Train

# ROC curve in training data
train_pred_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, newdata = predictor_data, type = "vote")
ROCTrain <- roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("NR", "R"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTrain

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, confusionTrain$overall["Accuracy"], confusionTrain$byClass["Balanced Accuracy"], confusionTrain$byClass["Sensitivity"], confusionTrain$byClass["Specificity"], MCC_Mechanistic_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_Test_Mechanistic)

## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusionTest <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "R")
confusionTest

# Calculate Matthews correlation coefficient
MCC_Mechanistic_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Mechanistic_Test

ROCTest <- roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], direction = "<", plot = F, print.auc = TRUE, levels = c("NR", "R"), col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

# Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, confusionTest$overall["Accuracy"], confusionTest$byClass["Balanced Accuracy"], confusionTest$byClass["Sensitivity"], confusionTest$byClass["Specificity"], MCC_Mechanistic_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
VanAllen_Out_RF_MechPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(VanAllen_Out_RF_MechPerformance, file = "./Objs/RF/VanAllen_Out_RF_MechPerformance.rda")

