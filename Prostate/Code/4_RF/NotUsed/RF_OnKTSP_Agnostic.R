########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: RF for predicting bladder cancer progression


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
library(randomForest)

## Agnostic
load("./Objs/KTSP/KTSP_STATs_Agnostic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")


usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup



predictor_data <- t(KTSP_STATs_Train_Agnostic)

tmp <- as.vector(table(usedTrainGroup))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

set.seed(333)
tuneRF(x=predictor_data, y=usedTrainGroup, plot = TRUE, improve = 0.01, ntreeTry = 1001, proximity = TRUE, sampsize = sampsizes)

set.seed(333)
RF_Agnostic_OnKTSP <- randomForest(x =predictor_data, y=usedTrainGroup, importance = TRUE, ntree = 1001, mtry = 5 ,proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
RF_Agnostic_OnKTSP
plot(RF_Agnostic_OnKTSP)

#png(filename = "./Figs/RF/AgnosticRF_VarImp.png", width = 2000, height = 2000, res = 300)
#varImpPlot(RF_Agnostic_OnKTSP, type=2, n.var=30, scale=FALSE, main="Agnostic RF variable Importance (Gini)")
#dev.off()


# ROC curve in training data
train_pred_votes_Agnostic_OnKTSP <- predict(RF_Agnostic_OnKTSP, newdata = predictor_data, type = "vote")
roc(usedTrainGroup, train_pred_votes_Agnostic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)

# predictions in the training data
PredResponse_Train <- predict(RF_Agnostic_OnKTSP, predictor_data, type="response")
PredVotes_Train <- predict(RF_Agnostic_OnKTSP, predictor_data, type="vote")

confusion_test <- confusionMatrix(PredResponse_Train, usedTrainGroup, positive = "Mets")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = PredResponse_Train, actuals = usedTrainGroup)
MCC_Train


################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_Test_Agnostic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Agnostic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Agnostic_OnKTSP <- predict(RF_Agnostic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Agnostic_OnKTSP <- predict(RF_Agnostic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Agnostic_OnKTSP, usedTestGroup, positive = "Mets")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Agnostic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Agnostic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)

