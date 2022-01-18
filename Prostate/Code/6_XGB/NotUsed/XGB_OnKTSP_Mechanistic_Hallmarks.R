rm(list = ls())

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

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
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")

Training <- t(KTSP_STATs_Train_Mechanistic)


usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

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

##
#save(Train_label, Test_label, file = "./Objs/XGB/Labels.rda")

## Convert to xgb.DMatrix
DataTrain <- xgb.DMatrix(as.matrix(Training), label = Train_label)
DataVal <- xgb.DMatrix(as.matrix(Validation), label = Val_label)
DataTest <- xgb.DMatrix(as.matrix(Testing), label = Test_label)

## Creat a watch list
watchlist <- list(train  = DataTrain, test = DataVal)

##########################
# Scale weight (to compensate for un-balanced class sizes)
Mets <- sum(Train_label == 1)
No_Mets <- sum(Train_label == 0)

## Make a 5-fold CV model to determine the best number of trees/iterations
# hyper_grid <- expand.grid(
#   eta = c(0.1,0.3),
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
# 
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
# save(hyper_grid, file = "./Objs/XGB/XGBMechanisticPairs_HyperGrid.rda")

# best nrow:99

## Make a list of model parameters
set.seed(333)

parameters <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.01,           #0.01 /0.01 /0.1/0.01   # default = 0.3, range: [0,1]
  gamma              = 0,             #0  # default = 0,   range: [0,∞]
  max_depth          = 3,             # 3 /6 /3
  min_child_weight   = 1,             #1   # default = 1,   range: [0,∞]
  subsample          = 0.3,           #0.3   # default = 1,   range: (0,1]
  colsample_bytree   = 1,           #1      # default = 1,   range: (0,1]
  colsample_bylevel  = 1,    #1  # default = 1,   range: (0,1]
  lambda             = 0,             #0   # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

## Make the final model
xgb.mechanistic_OnKTSP <- xgb.train(parameters, DataTrain, nrounds = 500, watchlist,  early_stopping_rounds = 50, scale_pos_weight = No_Mets/Mets)

#save(xgb.mechanistic_OnKTSP, file = "./Objs/XGB/MechanisticXGBoostModel.rda")
################################################
#################################################
## ROC curve and CM in the training data
XGB_prob_Train_mechanistic_OnKTSP <- predict(xgb.mechanistic_OnKTSP, DataTrain)
ROC_Train_mechanistic_OnKTSP <- roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, plot = FALSE, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE, main = "XGB ROC curve in the training cohort (Mechanistic)")

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
roc(Test_label, xgb_prob_test_mechanistic_OnKTSP, plot = F, print.auc = TRUE, levels = c("0", "1"), direction = "<", col = "black", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE, percent = TRUE, main = "XGB ROC curve in the testing cohort (Mechanistic)")


##############################################
##############

# create importance matrix
importance_matrix <- xgb.importance(model = xgb.mechanistic_OnKTSP)
write.csv(importance_matrix, file = "./Objs/XGB/XGB_importanceMatrix_MechanisticPairs.csv")


# variable importance plot
importance_matrix <- importance_matrix[order(importance_matrix$Gain, decreasing = TRUE), ]
MechanisticXGB_Pairs <- importance_matrix$Feature[1:50]
save(MechanisticXGB_Pairs, file = "./Objs/XGB/MechanisticXGB_Pairs.rda")

png(filename = "./Figs/XGB/XGB_Importance_MechanisticPairs.png", width = 3000, height = 2000, res = 300)
ggplot(data=importance_matrix[1:30,], aes(x = reorder(Feature, -Gain), y = Gain)) +
  geom_bar(stat="identity", colour = "black", fill = "lightgray") + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "XGB Feature Importance (Top 30) (Mechanistic Pairs)", x = "Features", y = "Information Gain")
dev.off()

####################
############################################################################
## Plot ggplot2 ROC curves for both the agnostic and mechanistic XGBoost 

### Prepare the legend
forLegend_XGB <- apply(rbind(
  ci(roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_mechanistic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Train_label, XGB_prob_Train_agnostic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_agnostic_OnKTSP, levels = c(0,1), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


####################
### ROC curves Using ggplot2

### Training
datTrn_XGB <- melt(data.frame(
  ## Training Group
  Training= usedTrainGroup1,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Training.Pairs=XGB_prob_Train_agnostic_OnKTSP,
  Mechanistic.Training.Pairs=XGB_prob_Train_mechanistic_OnKTSP
))

### Change Colnames
colnames(datTrn_XGB) <- c("Status", "XGB_type", "XGB_sum")


### Testing
datTst_XGB <- melt(data.frame(
  ## Testing group
  Testing= usedTestGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Testing.Pairs=xgb_prob_test_agnostic_OnKTSP,
  Mechanistic.Testing.Pairs=xgb_prob_test_mechanistic_OnKTSP
))
### Change Colnames
colnames(datTst_XGB) <- c("Status", "XGB_type", "XGB_sum")

### Combine
dat_XGB <- rbind(datTrn_XGB, datTst_XGB)
dat_XGB$Status <- as.numeric(dat_XGB$Status)-1

### Replace levels
levels(dat_XGB$XGB_type) <- gsub("\\.", "-", levels(dat_XGB$XGB_type))
levels(dat_XGB$XGB_type) <- paste(levels(dat_XGB$XGB_type), forLegend_XGB[c(3,1,4,2)])

###################
### Plot Curve
png("./Figs/XGB/CompareAUCggplot_Hallmarks.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic vs agnostic XGB (Pairs)"
legendTitle <- "Mechanistic VS Agnostic"
### Plot
basicplot_XGB_Pairs <- ggplot(dat_XGB, aes(d=Status, m=XGB_sum, color=XGB_type,
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
basicplot_XGB_Pairs
### Close device
dev.off()

save(basicplot_XGB_Pairs, file = "./Objs/XGB/BasicPlot_XGB_Pairs_Hallmarks.rda")


####################
# Load the pre-saved ROC stats for all the RF models (Mechanistic and agnostic with different N of pairs)
load("./Objs/XGB/XGBAUCStats.rda")

# Load the original data just for the phenotype variable
load("./Objs/XGB/Labels.rda")

## For the forest plot
AUCs_XGB_StrSampling <-  rbind(
  ci(roc(Train_label, XGB_prob_Train_mechanistic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_mechanistic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Train_label, XGB_prob_Train_agnostic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_agnostic_OnKTSP, levels = c(0,1), direction = "<")),
  ci(roc(Train_label, XGB_prob_Train_agnostic_OnKTSP_200, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_agnostic_OnKTSP_200, levels = c(0,1), direction = "<")),
  ci(roc(Train_label, XGB_prob_Train_agnostic_OnKTSP_500, levels = c(0,1), direction = "<")),
  ci(roc(Test_label, xgb_prob_test_agnostic_OnKTSP_500, levels = c(0,1), direction = "<"))
)


AUCs_XGB_StrSampling[,c(1,2,3)] <- AUCs_XGB_StrSampling[,c(2,1,3)]
colnames(AUCs_XGB_StrSampling) <- c("AUC", "CI_low", "CI_High")



AUCs_XGB_StrSampling <- as.data.frame(AUCs_XGB_StrSampling)

AUCs_XGB_StrSampling$model_type <- c("Mechanistic", "Mechanistic", "Agnostic", "Agnostic", "Agnostic200", "Agnostic200", "Agnostic500", "Agnostic500")
AUCs_XGB_StrSampling$data_type <- c("Training", "Testing", "Training", "Testing", "Training", "Testing", "Training", "Testing")
AUCs_XGB_StrSampling$algorithm <- rep("XGB", 8)
AUCs_XGB_StrSampling$approach <- rep("Str.Sampling", 8)  

save(AUCs_XGB_StrSampling, file = "./Objs/XGB/AUCs_XGB_StrSampling.rda")

