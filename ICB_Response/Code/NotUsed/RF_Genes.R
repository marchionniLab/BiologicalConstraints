#################################################################################
### Mohamed Omar
### 18/4/2019
### GOAL: Creating a rondom forest classifier for bladder cancer progression (Non-muscle invasive << Muscle-invasive)
### Using: TF_MiR genes (Mechanistic)
#################################################################################

# Clean the work space
rm(list = ls())

## settng the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Bladder")

## Load necessary libraries
library(randomForest)
library(pROC)
library(caret)
library(limma)
library(mltools)

## Load the data
load("./Objs/ProcessedDatasets.rda")


RegImmEffec <- read.delim("./RegImmEffec.txt")
RegImmEffec <- as.character(RegImmEffec[-1, ])
RegImmEffec <- gsub("-.+", "", RegImmEffec)


BcellAct <- read.delim("./BcellAct.txt")
BcellAct <- as.character(BcellAct[-1,])
BcellAct <- gsub("-.+","",BcellAct)

WNT_B_Cat <- read.delim("./WNT_B_Cat.txt")
WNT_B_Cat <- as.character(WNT_B_Cat[-1,])
WNT_B_Cat <- gsub("-.+", "", WNT_B_Cat)

# CytoProd <- read.delim("./CytoProd.txt")
# CytoProd <- as.character(CytoProd[-1,])
# CytoProd <- gsub("-.+","",CytoProd)

#Genes <- read.delim("./geneset1.txt")
#Genes <- as.character(Genes[-1,])
#Genes <- gsub("-.+", "", Genes)

AllImmune <- c(RegImmEffec, BcellAct, WNT_B_Cat)
AllImmune <- AllImmune[!duplicated(AllImmune)]

keepGns <- intersect(AllImmune, rownames(Expr_NP))


### Normalization
usedTrainMat <- Expr_NP[keepGns, ]
usedTestMat <- Melanoma_Mat[keepGns, ]


### Associated groups
usedTrainGroup <- Pheno_NP$SR
usedTestGroup <- Melanoma_Group

names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

levels(usedTestGroup) <- c("no", "yes")

################################################################################
################################################################################
################################################################################


## transpose the matrix
predictor_data <- t(usedTrainMat)
predictor_names <- c(as.vector(rownames(usedTrainMat))) #gene symbol
colnames(predictor_data) <- predictor_names

## Setting the variable we are trying to predict as our target variable. In this case, it is Progression status.
## train group here is just the column containing the phenptype of interest (Progression vs NoProgression) from the phenotype table

#usedTrainGroup <- ordered(usedTrainGroup, levels=c("NoProgression", "Progression"))
target <- usedTrainGroup

## Finally we run the RF algorithm. 
## NOTE: use an ODD number for ntree. This is because when the forest is used on test data, ties are broken randomly. Having an odd number of trees avoids this issue.
## Use down-sampling to attempt to compensate for unequal class-sizes (less progression than noProgression).
tmp <- as.vector(table(target))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)

######
# # Find the important variables
# set.seed(333)
# Boruta <- Boruta(predictor_data, usedTrainGroup, maxRuns = 500, pValue = 0.01, getImp = getImpRfGini)
# Boruta
# plot(Boruta)
# Variables <- getSelectedAttributes(Boruta, withTentative = F)



# Tunning RF model (to find the best mtry)
set.seed(333)
tuneRF(x=predictor_data, y=target, plot = TRUE, improve = 0.01, ntreeTry = 501, proximity = TRUE, sampsize = c(2,2))

set.seed(333)
RF_Mechanistic <- randomForest(x =predictor_data, y=target,importance = TRUE, ntree = 501, mtry =40 ,proximity=TRUE, na.action = na.omit, sampsize = c(2,2))

print(RF_Mechanistic)
plot(RF_Mechanistic)


## RandomForest calculates an importance measures for each variable.
#rf_importances <- randomForest::importance(RF_Mechanistic, scale=FALSE)
#rf_importances <- rf_importances[order(rf_importances[,4], decreasing = TRUE), ]

## Predict in the training data
Prediction_Train_Mechanistic <- predict(RF_Mechanistic, predictor_data)
confusion_train <- confusionMatrix(Prediction_Train_Mechanistic, usedTrainGroup, positive = "yes")
confusion_train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = Prediction_Train_Mechanistic, actuals = usedTrainGroup)
MCC_Train

## Create a representation of the top 50 variables categorized by importance.
#png("./Figs/RF/RF_varImp_Mechanistic.png", width = 2000, height = 2000, res = 300)
#varImpPlot(RF_Mechanistic, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors (Mechanistic)")
#dev.off()

## An MDS plot provides a sense of the separation of classes.
#png("./Figs/RF/MDS_Mechanistic.png", width = 2000, height = 2000, res = 300)
#target_labels=as.vector(usedTrainGroup)
#MDSplot(RF_Mechanistic, usedTrainGroup, k=2, xlab="", ylab="", pch=target_labels, palette=c("red", "blue"), main="MDS plot")
#dev.off()



# ROC curve in training data
train_pred_votes_Mechanistic <- predict(RF_Mechanistic, newdata = predictor_data, type = "vote")
roc(target, train_pred_votes_Mechanistic[,2], plot = F, print.auc=TRUE, levels = c("no", "yes"), direction = "<", col="blue", lwd=2, grid=TRUE)


################################################################################ 
################################################################################ 
## Testing the classifier using the testing set

## transpose the matrix. RandomForests expects the predictor variables (genes) to be represented as columns instead of rows. Finally, assign gene symbol as the predictor name.
predictor_data2 <- t(usedTestMat)
predictor_names2 <- c(as.vector(rownames(usedTestMat))) #gene symbol
colnames(predictor_data2) <- predictor_names2

## Extract predictor (gene) names from RF model built above and subset the test data to just these predictors
RF_predictor_names <- rownames(RF_Mechanistic$importance)
predictor_data2 <- predictor_data2[,RF_predictor_names]

## Run the test data through forest!
RF_predictions_responses_Mechanistic <- predict(RF_Mechanistic, predictor_data2, type="response")
RF_predictions_votes_Mechanistic <- predict(RF_Mechanistic, predictor_data2, type="vote")

### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic, usedTestGroup, positive = "yes")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic, actuals = usedTestGroup)
MCC_Test

## ROC curve and AUC
roc(usedTestGroup, RF_predictions_votes_Mechanistic[,2], plot = F, print.auc = TRUE, levels = c("no", "yes"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE)

#########
## Save RF classifier
#save(RF_Mechanistic, file = "./Objs/RF/RF_Classifier_Mechanistic")

## Save RF importances
#save(rf_importances, file = "./Objs/RF/RF_Importance_Mechanistic")

