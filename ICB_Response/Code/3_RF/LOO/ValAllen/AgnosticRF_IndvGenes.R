# Clean the work space
rm(list = ls())

## Load necessary libraries
library(randomForest)
require(switchBox)
library(pROC)
library(caret)
library(limma)
library(mltools)

## Load the data
load("./Objs/icbData_VanAllenOut.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")
usedTestMat <- testMat

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))


################################################################################
################################################################################

## transpose the matrix
predictor_data <- t(usedTrainMat)
predictor_names <- c(as.vector(rownames(usedTrainMat))) #gene symbol
colnames(predictor_data) <- predictor_names

## Setting the variable we are trying to predict as our target variable. In this case, it is Progression status.
## train group here is just the column containing the phenptype of interest (Progression vs NoProgression) from the phenotype table

#usedTrainGroup <- ordered(usedTrainGroup, levels=c("NR", "R"))
target <- usedTrainGroup

## Finally we run the RF algorithm. 
## NOTE: use an ODD number for ntree. This is because when the forest is used on test data, ties are broken randomly. Having an odd number of trees avoids this issue.
## Use down-sampling to attempt to compensate for unequal class-sizes (less progression than noProgression).
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

##############

# Tunning RF model (to find the best mtry)
set.seed(333)
tuneRF(x=predictor_data, y=target, plot = TRUE, improve = 0.01, ntreeTry = 1000, proximity = TRUE, sampsize = sampsizes, na.action = na.omit)

set.seed(333)
RF_Agnostic <- randomForest(x =predictor_data, y=target, importance = TRUE, ntree = 1000, mtry = 65 ,proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
print(RF_Agnostic)

# ROC curve in training data
train_pred_votes_Agnostic <- predict(RF_Agnostic, newdata = predictor_data, type = "vote")
ROCTrain <- roc(usedTrainGroup, train_pred_votes_Agnostic[,2], plot = F, print.auc=TRUE, levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE, ci = T)
ROCTrain

### Predict in the training data
train_pred_Response_Agnostic <- predict(RF_Agnostic, newdata = predictor_data, type = "response")

confusionTrain <- confusionMatrix(train_pred_Response_Agnostic, usedTrainGroup, positive = "R")
confusionTrain

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = train_pred_Response_Agnostic, actuals = usedTrainGroup)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, confusionTrain$overall["Accuracy"], confusionTrain$byClass["Balanced Accuracy"], confusionTrain$byClass["Sensitivity"], confusionTrain$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(usedTestMat)
predictor_names2 <- c(as.vector(rownames(usedTestMat))) #gene symbol
colnames(predictor_data2) <- predictor_names2

## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Agnostic$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Agnostic <- predict(RF_Agnostic, predictor_data2, type="response")
RF_predictions_votes_Agnostic <- predict(RF_Agnostic, predictor_data2, type="vote")


### Predict in the testing data
confusionTest <- confusionMatrix(RF_predictions_responses_Agnostic, usedTestGroup, positive = "R")
confusionTest

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Agnostic, actuals = usedTestGroup)
MCC_Test

## ROC curve and AUC
ROCTest <- roc(usedTestGroup, RF_predictions_votes_Agnostic[,2], plot = F, print.auc = TRUE, levels = c("NR", "R"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
ROCTest

#######
# Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, confusionTest$overall["Accuracy"], confusionTest$byClass["Balanced Accuracy"], confusionTest$byClass["Sensitivity"], confusionTest$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
VanAllen_Out_RF_IndvGenes_AgnosticPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(VanAllen_Out_RF_IndvGenes_AgnosticPerformance, file = "./Objs/RF/VanAllen_Out_RF_IndvGenes_AgnosticPerformance.rda")

#########
