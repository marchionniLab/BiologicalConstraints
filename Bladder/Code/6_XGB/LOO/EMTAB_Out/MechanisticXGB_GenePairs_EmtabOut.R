# Mechanistic XGB
# Gene pairs
# Cross-Study validation (EMTAB out)


rm(list = ls())

setwd("/Volumes/Macintosh/Research/Projects/Bladder")

## Load necessary packages
library(xgboost)
library(pdp)
library(lime)
library(preprocessCore)
library(limma)
library(pROC)
library(caret)
library(DiagrammeR)
library(ggplot2)
library(xgboostExplainer)
library(dplyr)
library(Ckmeans.1d.dp)
library(mltools)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
library(plotROC)

## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_EMTABOut.rda")
load("./Objs/ProgressionDataGood_EMTABOut.rda")

Training <- t(KTSP_STATs_Train_Mechanistic)


usedTrainGroup <- trainGroup
usedTestGroup <- testGroup


#####################
## Here, WE divide the training data into "actual training" and another validation
# This validation data will be used in the "WatchList" parameter. It is independent of the testing data.
set.seed(333)
inds <- createDataPartition(usedTrainGroup, p = 0.7, list = F)
Training1 <- Training[inds, ]
Validation <- Training[-inds, ]

usedTrainGroup1 <- usedTrainGroup[inds]
usedValGroup <- usedTrainGroup[-inds]

table(usedTrainGroup1)
table(usedValGroup)


## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup1) <- rownames(Training1)
all(rownames(Training1) == names(usedTrainGroup1))

all(rownames(Validation) == names(usedValGroup))

## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training1)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_train <- cbind(Training, usedTrainGroup1)

## The same for validation
Validation <- as.data.frame(Validation)
#usedTrainGroup <- as.data.frame(usedTrainGroup)
Data_val <- cbind(Validation, usedValGroup)


########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Mechanistic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

###########################################################
## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_train$usedTrainGroup1)  
levels(Data_train$usedTrainGroup1) <- c(0,1) 
Data_train$usedTrainGroup <- factor(Data_train$usedTrainGroup1, levels = c(0,1)) # 0=No,1= Yes 
Train_label <- Data_train$usedTrainGroup1
Train_label <- as.vector(Train_label)
table(Train_label)


## The same for validation
table(Data_val$usedValGroup)  
levels(Data_val$usedValGroup) <- c(0,1) 
Data_val$usedValGroup <- factor(Data_val$usedValGroup, levels = c(0,1)) # 0=No,1= Yes 
Val_label <- Data_val$usedValGroup
Val_label <- as.vector(Val_label)
table(Val_label)

## Combine both the Expression matrix and the phenotype into one matrix
Testing <- as.data.frame(Testing)
Data_test <- cbind(Testing, usedTestGroup)

## Converting classes from Progression/NoProgression Format to 0-1 Format
table(Data_test$usedTestGroup)  
levels(Data_test$usedTestGroup) <- c(0,1)
Data_test$usedTestGroup <- factor(Data_test$usedTestGroup, levels = c(0,1)) #0=No, 1=Yes
Test_label <- Data_test$usedTestGroup
Test_label <- as.vector(Test_label)
table(Test_label)



## Convert to xgb.DMatrix
DataTrain <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataVal <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
DataTest <- xgb.DMatrix(as.matrix(Testing), label = Test_label)

## Creat a watch list
watchlist <- list(train  = DataTrain, test = DataVal)

##########################
# Scale weight (to compensate for un-balanced class sizes)
Progression <- sum(Train_label == 1)
NoProgression <- sum(Train_label == 0)

## Make a 5-fold CV model to determine the best number of trees/iterations
# hyper_grid <- expand.grid(
#   eta = c(0.1, 0.3),
#   max_depth = c(1,3,6),
#   min_child_weight = 1,
#   subsample = c(0.7,0.8,0.9,1),
#   colsample_bytree = c(0.7,0.8,0.9,1),
#   colsample_bylevel = c(0.7,0.8,0.9,1),
#   gamma = c(0,1),
#   lambda = c(0,0.5,1),
#   alpha = 0,
#   optimal_trees = 0,               # a place to dump results
#   max_AUC = 0                     # a place to dump results
# )
# 
# 
# ##########################
# # grid search
# set.seed(333)
# 
# for(i in 1:nrow(hyper_grid)) {
# 
#   # create parameter list
#   params <- list(
#     eta = hyper_grid$eta[i],
#     max_depth = hyper_grid$max_depth[i],
#     min_child_weight = hyper_grid$min_child_weight[i],
#     subsample = hyper_grid$subsample[i],
#     colsample_bytree = hyper_grid$colsample_bytree[i],
#     colsample_bylevel = hyper_grid$colsample_bylevel[i],
#     gamma = hyper_grid$gamma[i],
#     lambda = hyper_grid$lambda[i],
#     alpha = hyper_grid$alpha[i]
#   )
# 
#   # reproducibility
#   set.seed(333)
# 
#   # train model
#   xgb.tune <- xgb.cv(
#     params = params,
#     data = DataTrain,
#     nrounds = 500,
#     nfold = 5,
#     objective = "binary:logistic",
#     eval_metric        = "auc",# for regression models
#     verbose = 1,               # silent,
#     early_stopping_rounds = 50,
#     scale_pos_weight = NoProgression/Progression
#   )
# 
#   # add max training AUC and trees to grid
#   hyper_grid$optimal_trees[i] <- which.max(xgb.tune$evaluation_log$test_auc_mean)
#   hyper_grid$max_AUC[i] <- max(xgb.tune$evaluation_log$test_auc_mean)
# }
# 
# hyper_grid %>% arrange(max_AUC)
# View(hyper_grid)
# save(hyper_grid, file = "./Objs/XGB/CV_XGBMechanisticPairs_HyperGrid.rda")

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1    # default = 0.3, range: [0,1]
  gamma              = 0,             #0  # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #0    # default = 1,   range: [0,∞]
  subsample          = 0.5,           #0.5      # default = 1,   range: (0,1]
  colsample_bytree   = 1,           #1     # default = 1,   range: (0,1]
  colsample_bylevel  = 1,    #1  # default = 1,   range: (0,1]
  lambda             = 1,             #0.1 # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
  )

## Make the final model
xgb.mechanistic_OnKTSP <- xgb.train(parameters, DataTrain, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = NoProgression/Progression)

################################################
#################################################
## ROC curve and CM in the training data
XGB_prob_Train_mechanistic_OnKTSP <- predict(xgb.mechanistic_OnKTSP, DataTrain)
ROC_Train_mechanistic_OnKTSP <- roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, plot = FALSE, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the training cohort (Mechanistic)")
ROC_Train_mechanistic_OnKTSP

thr_Train <- coords(roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, levels = c("0", "1"), direction = "<"), transpose = TRUE, "best")["threshold"]
thr_Train

## Convert predicted probabilities to binary outcome
prediction_Train_mechanistic_OnKTSP <- as.numeric(XGB_prob_Train_mechanistic_OnKTSP > thr_Train)
print(head(prediction_Train_mechanistic_OnKTSP))

Train_label <- factor(Train_label, levels = c(0,1))
prediction_Train_mechanistic_OnKTSP <- factor(prediction_Train_mechanistic_OnKTSP, levels = c(0,1))

# Confusion matrix in the training data
CM_Train <- confusionMatrix(prediction_Train_mechanistic_OnKTSP, Train_label, positive = "1")
CM_Train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = prediction_Train_mechanistic_OnKTSP, actuals = Train_label)
MCC_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROC_Train_mechanistic_OnKTSP$ci, CM_Train$overall["Accuracy"], CM_Train$byClass["Balanced Accuracy"], CM_Train$byClass["Sensitivity"], CM_Train$byClass["Specificity"], MCC_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#####################################
######################################

## Predict in the Test data
xgb_prob_test_mechanistic_OnKTSP <- predict(xgb.mechanistic_OnKTSP, DataTest)

## Convert predicted probabilities to binary outcome
prediction_Test_mechanistic_OnKTSP <- as.numeric(xgb_prob_test_mechanistic_OnKTSP > thr_Train)
print(head(prediction_Test_mechanistic_OnKTSP))

Test_label <- factor(Test_label, levels = c(0,1))
prediction_Test_mechanistic_OnKTSP <- factor(prediction_Test_mechanistic_OnKTSP, levels = c(0,1))

## Confusion matrix
CM_Test <- confusionMatrix(prediction_Test_mechanistic_OnKTSP, Test_label, positive = "1")
CM_Test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = prediction_Test_mechanistic_OnKTSP, actuals = Test_label)
MCC_Test

## ROC curve and AUC
ROCTest <- roc(Test_label, xgb_prob_test_mechanistic_OnKTSP, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, main = "XGB ROC curve in the testing cohort (Mechanistic)")
ROCTest

# Put the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, CM_Test$overall["Accuracy"], CM_Test$byClass["Balanced Accuracy"], CM_Test$byClass["Sensitivity"], CM_Test$byClass["Specificity"], MCC_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
EMTAB_Out_XGB_MechPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(EMTAB_Out_XGB_MechPerformance, file = "./Objs/XGB/EMTAB_Out_XGB_MechPerformance.rda")


##############################################
#############

# create importance matrix
importance_matrix <- xgb.importance(model = xgb.mechanistic_OnKTSP)
#write.csv(importance_matrix, file = "./Objs/BEST_XGB_importanceMatrix_Mechanistic.csv")
# 
# 
# # variable importance plot
# importance_matrix <- importance_matrix[order(importance_matrix$Gain, decreasing = TRUE), ]
# 
# png(filename = "./Figs/BEST_XGB_Importance_Mechanistic_GSE57813Out.png", width = 3000, height = 2000, res = 300)
# ggplot(data=importance_matrix[1:30,], aes(x = reorder(Feature, -Gain), y = Gain)) +
#   geom_bar(stat="identity", colour = "black", fill = "lightgray") +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
#   labs(title = "XGB Feature Importance (Top 30) (Mechanistic)", x = "Features", y = "Information Gain")
# dev.off()
# 
# ####################
############################################################################
## Plot ggplot2 ROC curves for both the agnostic and mechanistic XGBoost 

### Prepare the legend
forLegend_XGB <- apply(rbind(
  ci(roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP)),
  ci(roc(Test_label, xgb_prob_test_mechanistic_OnKTSP)),
  ci(roc(Train_label, XGB_prob_Train_agnostic_OnKTSP)),
  ci(roc(Test_label, xgb_prob_test_agnostic_OnKTSP))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


####################
### ROC curves Using ggplot2

### Training
datTrn_XGB <- melt(data.frame(
  ## Training Group
  Training=factor(usedTrainGroup1, levels=c("NoProgression", "Progression")),
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Training.Pairs=XGB_prob_Train_agnostic_OnKTSP,
  Mechanistic.Training.Pairs=XGB_prob_Train_mechanistic_OnKTSP
))

### Change Colnames
colnames(datTrn_XGB) <- c("Progression", "XGB_type", "XGB_sum")


### Testing
datTst_XGB <- melt(data.frame(
  ## Testing group
  Testing=factor(usedTestGroup, levels=c("NoProgression", "Progression")),
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Testing.Pairs=xgb_prob_test_agnostic_OnKTSP,
  Mechanistic.Testing.Pairs=xgb_prob_test_mechanistic_OnKTSP
))
### Change Colnames
colnames(datTst_XGB) <- c("Progression", "XGB_type", "XGB_sum")

### Combine
dat_XGB <- rbind(datTrn_XGB, datTst_XGB)
dat_XGB$Progression <- as.numeric(dat_XGB$Progression)-1

### Replace levels
levels(dat_XGB$XGB_type) <- gsub("\\.", "-", levels(dat_XGB$XGB_type))
levels(dat_XGB$XGB_type) <- paste(levels(dat_XGB$XGB_type), forLegend_XGB[c(3,1,4,2)])

###################
### Plot Curve
png("./Figs/XGB/CompareAUCggplot_EmtabOut.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic vs agnostic XGB (EMTAB-4321 for validation)"
legendTitle <- "Mechanistic VS Agnostic"
### Plot
basicplot_XGB_EMTABOut <- ggplot(dat_XGB, aes(d=Progression, m=XGB_sum, color=XGB_type,
                                                 linetype = XGB_type)) +
  geom_roc(cutoffs.at = seq(0.1,1,0.1)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_XGB_EMTABOut
### Close device
dev.off()

save(basicplot_XGB_EMTABOut, file = "./Objs/XGB/BasicPlot_XGB_EMTABOut.rda")


####################

