############################################################################
## Moahmed Omar
## 06/07/2019
## Goal: Creating a Deep Neural Network to predict the metastasis of prostate cancer 
## Using Keras/Tensor Flow
#############################################################################

## Clean the Working Space
rm(list = ls())

## Setting the Working Directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

## Load Necessary Libraries
library(limma)
library(caret)
library(pROC)
library(DMwR)
library(keras)
#install_keras()
####################################

## Load data
load("./Objs/MetastasisData.rda")

## Load the RF important genes (ordered according to MeanDecreaseGini)
load("./Objs/RandomForest/RandomForest_Importance.rda")

## Keep only the top 500 genes
keepGns <- rownames(rf_importances)[1:500]

## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(mixTrainMat)[keepGns, ]
boxplot(usedTrainMat, outline = FALSE)

usedTestMat <- normalizeBetweenArrays(mixTestMat)[keepGns, ]
boxplot(usedTestMat, outline = FALSE)

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

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
Data_train$usedTrainGroup <- as.integer(Data_train$usedTrainGroup)

Train_label <- Data_train$usedTrainGroup-1
table(Train_label)

########################################################

## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
Data_test$usedTestGroup <- as.integer(Data_test$usedTestGroup)
table(Data_test$usedTestGroup)
Test_label <- Data_test$usedTestGroup-1
table(Test_label)

################################################################
## Before making the NN we need to put all Data in a Matrix Format

Data_train <- as.matrix(Data_train)
Data_train[,1:500] <- normalize(Data_train[,1:500])
Data_train[,501] <- as.numeric(Data_train[,501])-1
table(Data_train[,501])  ## 0: Mets, 1: No_Mets

## The same For the Test Data
Data_test <- as.matrix(Data_test)
Data_test[,1:500] <- normalize(Data_test[,1:500])
Data_test[,501] <- as.numeric(Data_test[,501])-1
table(Data_test[,501])

Train_target <- Data_train[,501]
Test_target <- Data_test[,501]

## One Hot Encoding
TrainLabels <- to_categorical(Train_target) 
TestLabels <- to_categorical(Test_target)

## Performing Oversampling using SMOTE (to compensate for the unbalanced classes)
set.seed(333)
Data_train <- as.data.frame(Data_train)
Data_train[,501] <- as.factor(Data_train[,501])
Data_train <- SMOTE(usedTrainGroup~., Data_train, perc.over = 650, perc.under = 120)
table(Data_train$usedTrainGroup)

Data_train <- as.matrix(Data_train)
Data_train[, 501] <- as.numeric(Data_train[,501])
table(Data_train[,501])

Train_target <- Data_train[,501]
TrainLabels <- to_categorical(Train_target) 

########################################################################
## Creating the Network

set.seed(333)
Model <- keras_model_sequential()

Model %>% 
  layer_dense(units = 500, activation = "relu", input_shape = c(500)) %>%
  layer_dropout(0.4) %>%
  layer_dense(units = 200, activation = "relu") %>%
  layer_dropout(0.2, seed = 333) %>%
  layer_dense(units = 50, activation = "relu") %>%
  layer_dense(units = 2, activation = "softmax")
summary(Model)

Model %>% 
  compile(loss="binary_crossentropy", optimizer= "adam", metrics="accuracy")

History <- Model %>% fit(Data_train[,1:500], TrainLabels, epoch=300, batch_size=32, validation_split=0.3)

## Plot the model (epoch/iterations VS accuracy)
png("./Figs/Model_accuracy_PrognGns.png", width = 3000, height = 1500, res = 400)
plot(History)
dev.off()

## Evaluation of the NN Model on the Testing set
Model %>% evaluate(Data_test[,1:500], TestLabels, batch_size=64) 

## Predict Probabilities in the Test Data
prob <- Model %>% predict_proba(Data_test[,1:500], batch_size = 64)

## Predict classes in the Test Data
pred <- Model %>% predict_classes(Data_test[,1:500], batch_size = 64)
table(pred)

## Create a confusion table 
table1 <- table(Predicted = pred, Actual=Test_target)
table1

## ROC and AUC in the Testing Data
png("./Figs/ROC_Test_PrognGns.png", width = 2000, height = 2000, res = 300)
roc(Test_target, pred, plot = TRUE, print.auc=TRUE, levels = c("1", "0"), direction = ">", col="blue", lwd=2, grid=TRUE, main="Keras ROC test using Prognosis genes (631 genes)")
dev.off()
