########################################################################## 
## Mohamed Omar
## 03/06/2019
## Goal: RF for predicting bladder cancer progression


###########################################################################

## Clean the work environment
rm(list = ls())

## Set the working directory
setwd("/Volumes/Macintosh/Research/Projects/ICB_Response")

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(randomForest)
library(varSelRF)


## Load data
load("./Objs/KTSP_Stats.rda")
load("./Objs/ProcessedDatasets.rda")


usedTrainGroup <- Pheno_NP$SR
usedTestGroup <- Melanoma_Group
levels(usedTestGroup) <- c("no", "yes")

predictor_data <- t(KTSP_STATs_Train_Mechanistic)

tmp <- as.vector(table(usedTrainGroup))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)



set.seed(333)
tuneRF(x=predictor_data, y=usedTrainGroup, plot = TRUE, improve = 0.01, ntreeTry = 500, proximity = TRUE, sampsize = c(10,5), na.action = na.omit)

set.seed(333)
RF_Mechanistic_OnKTSP <- randomForest(x =predictor_data, y=usedTrainGroup, importance = TRUE, ntree = 500, mtry = 1 ,proximity=TRUE, na.action = na.omit, sampsize = c(10,5))
RF_Mechanistic_OnKTSP
#plot(RF_Mechanistic_OnKTSP)


#VAR <- varSelRF(xdata = predictor_data, Class = usedTrainGroup, ntree = 501, ntreeIterat = 20, whole.range = FALSE, keep.forest = TRUE)

#SelectedVars <- VAR$selected.vars

## RandomForest calculates an importance measures for each variable.
#rf_importances <- randomForest::importance(RF_Mechanistic_OnKTSP, scale=FALSE)
#rf_importances <- rf_importances[order(rf_importances[,4], decreasing = TRUE), ]
#MechanisticRF_Pairs <- rownames(rf_importances)[1:50]
#save(MechanisticRF_Pairs, file = "./Objs/RF/MechanisticRF_Pairs.rda")

# Plot importance
#png(filename = "./Figs/RF/MechanisticRF_VarImp.png", width = 2000, height = 2000, res = 300)
#varImpPlot(RF_Mechanistic_OnKTSP, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 pairs")
#dev.off()

# ROC curve in training data
train_pred_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, newdata = predictor_data, type = "vote")
roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

# predictions in the training data
train_pred_response_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data, type="response")

### Predict in the testing data
confusion_train <- confusionMatrix(train_pred_response_Mechanistic_OnKTSP, usedTrainGroup, positive = "yes")
confusion_train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = train_pred_response_Mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Train

################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_Melanoma_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

############################################################################
# HUGO dataset
usedTestGroup <- Pheno_HUGO$Response
levels(usedTestGroup) <- c("no", "yes")
## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_HUGO_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

##################################################################
# TCGA dataset
usedTestGroup <- Pheno_TCGA$Response
levels(usedTestGroup) <- c("no", "yes")
## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_TCGA_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

##################################################################
# RIAZ_Pre_aPD1 dataset
usedTestGroup <- Pheno_RIAZ_Pre_aPD1$Response
levels(usedTestGroup) <- c("no", "yes")
## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_RIAZ_Pre_aPD1_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

##################################################################
# RIAZ_On_aPD1 dataset
usedTestGroup <- Pheno_RIAZ_On_aPD1$Response
levels(usedTestGroup) <- c("no", "yes")
## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_RIAZ_On_aPD1_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

##################################################################3
# MGH_On_aPD1 dataset
usedTestGroup <- Pheno_MGH_On_aPD1$Response
levels(usedTestGroup) <- c("no", "yes")
## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_MGH_On_aPD1_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

##################################################################3
# MGH_On_aCTLA4 dataset
usedTestGroup <- Pheno_MGH_On_aCTLA4$Response
levels(usedTestGroup) <- c("no", "yes")
## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_MGH_On_aCTLA4_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

##################################################################3
# VanAllen dataset
usedTestGroup <- Pheno_VanAllen$Response
levels(usedTestGroup) <- c("no", "yes")
## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(KTSP_STATs_VanAllen_Mechanistic)


## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic_OnKTSP$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="response")
RF_predictions_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data2, type="vote")


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

##################################################################3




