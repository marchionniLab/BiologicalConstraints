########################################################################### 
## Mohamed Omar
## 27/11/2019
## Goal: SVM (Poly) for predicting prostate cancer metastasis
## Mechanistic (Using TF-targets genes)
## Cross study validation: JHU Out
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
library(mltools)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(plotROC)
#######################################################################

## Load data
load("./Objs/KTSP/LOO/KTSP_STATs_Mechanistic_JHUOut.rda")
load("./Objs/LOO/MetastasisData_JHUOut.rda")



### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

### Transpose usedTrainMat (making samples as rows instead of columns)
Training <- t(KTSP_STATs_Train_Mechanistic)


## Making sure that sample names are identical in both Training and usedTrainGroup
names(usedTrainGroup) <- rownames(Training)
all(rownames(Training) == names(usedTrainGroup))


## Combining the expression matrix and the phenotype in one data frame
Training <- as.data.frame(Training)
Data_train <- cbind(Training, usedTrainGroup)

########################################################
# Transpose usedTestMat and make the sample names identical 
Testing <- t(KTSP_STATs_Test_Mechanistic)

names(usedTestGroup) <- rownames(Testing)
all(rownames(Testing) == names(usedTestGroup))

#################################################################
## Oversampling of the training data set to compensate for the un-balanced classes
# set.seed(333)
# Data_train <- as.data.frame(Data_train)
# Data_train[,629] <- as.factor(Data_train[,629])
# Over_Train <- SMOTE(usedTrainGroup~., data = Data_train, perc.over = 300, perc.under = 134)
# table(Over_Train[,629])

######
control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary)

###########################################################################
###########################################################################

## Model: SVM Poly

## 5-fold cross validation repeated 5 times (to find the best parameters)
set.seed(333)
fit.svmPoly <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=control, metric = "ROC")
fit.svmPoly

## Training using all data (using the best parameters)
Grid <- expand.grid(degree = 2, scale = 0.001, C = 0.25)
set.seed(333)
fit.svmPoly_mechanistic_OnKTSP <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = c("ROC"))
fit.svmPoly_mechanistic_OnKTSP

##########################################################
## Predict in the training data

## ROC stat for the training data
train_pred_prob_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Training, type = "prob")
roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2], plot = F, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Mechanistic)")

# The best threshold
thr <- coords(roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr <- thr$threshold
thr

# Predict classes according to the best threshold
train_pred_classes_svmPoly_mechanistic_OnKTSP <- ifelse(train_pred_prob_svmPoly_mechanistic_OnKTSP[,2] > thr, "Mets", "No_Mets")
table(train_pred_classes_svmPoly_mechanistic_OnKTSP)
# Convert to factor
train_pred_classes_svmPoly_mechanistic_OnKTSP <- factor(train_pred_classes_svmPoly_mechanistic_OnKTSP, levels = c("No_Mets", "Mets"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_train_svmPoly_OnKTSP <- confusionMatrix(train_pred_classes_svmPoly_mechanistic_OnKTSP, usedTrainGroup, positive = "Mets")
Confusion_train_svmPoly_OnKTSP

MCC_Train <- mltools::mcc(pred = train_pred_classes_svmPoly_mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Train

#############################################################
## Predict in the testing data 

## ROC/AUC in the Testing set
test_pred_prob_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Testing, type = "prob")
roc(usedTestGroup, test_pred_prob_svmPoly_mechanistic_OnKTSP[,2], plot = F, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Mechanistic)")

# Predict classes according to the best threshold
test_pred_classes_svmPoly_mechanistic_OnKTSP <- ifelse(test_pred_prob_svmPoly_mechanistic_OnKTSP[,2] > thr, "Mets", "No_Mets")
table(test_pred_classes_svmPoly_mechanistic_OnKTSP)
# Convert to factor
test_pred_classes_svmPoly_mechanistic_OnKTSP <- factor(test_pred_classes_svmPoly_mechanistic_OnKTSP, levels = c("No_Mets", "Mets"))

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly_OnKTSP <- confusionMatrix(test_pred_classes_svmPoly_mechanistic_OnKTSP, usedTestGroup, positive = "Mets")
Confusion_test_svmPoly_OnKTSP

MCC_Mechanistic <- mltools::mcc(pred = test_pred_classes_svmPoly_mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Mechanistic


## Top predictors
#Importance_SVMPoly_Mechanistic <- varImp(fit.svmPoly_mechanistic_OnKTSP, scale = TRUE)
# png("./Figs/SVM/SVM_varImp_Mechanistic_OnKTSP.png", width = 2000, height = 2000, res = 200)
# plot(Importance_SVMPoly_Mechanistic, top = 20, main = "Mechanistic SVM top 20 genes")
# dev.off()
# 
#Importance_SVMPoly_Mechanistic <- Importance_SVMPoly_Mechanistic$importance
#Importance_SVMPoly_Mechanistic <- Importance_SVMPoly_Mechanistic[order(Importance_SVMPoly_Mechanistic$Mets, decreasing = TRUE),]
#Importance_SVMPoly_Mechanistic <- Importance_SVMPoly_Mechanistic[!(Importance_SVMPoly_Mechanistic$Mets == 0), ]

# genes_SVMPoly <- rownames(Importance_SVMPoly)[1:50]
# save(Importance_SVMPoly, file = "./Objs/SVM/genes_SVMPloy_Mechanistic_OnKTSP.rdata")

########################
############################################################################
## Plot ggplot2 ROC curves for both the agnostic and mechanistic XGBoost 

### Prepare the legend
forLegend_SVM <- apply(rbind(
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2])),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_mechanistic_OnKTSP[,2])),
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP[,2])),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_agnostic_OnKTSP[,2]))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


####################
### ROC curves Using ggplot2

### Training
datTrn_SVM <- melt(data.frame(
  ## Training Group
  Training= usedTrainGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  # Agnostic pairs
  Agnostic.Training.Pairs = train_pred_prob_svmPoly_agnostic_OnKTSP[,2],
  # Mechanistic pairs
  Mechanistic.Training.Pairs = train_pred_prob_svmPoly_mechanistic_OnKTSP[,2]))
### Change Colnames
colnames(datTrn_SVM) <- c("Status", "SVM_type", "SVM_sum")


### Testing
datTst_SVM <- melt(data.frame(
  ## Testing group
  Testing= usedTestGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  # Agnostic.testing.Pairs
  Agnostic.Testing.Pairs = test_pred_prob_svmPoly_agnostic_OnKTSP[,2],
  # Mechanistic.Testing.Pairs
  Mechanistic.Testing.Pairs = test_pred_prob_svmPoly_mechanistic_OnKTSP[,2]))
### Change Colnames
colnames(datTst_SVM) <- c("Status", "SVM_type", "SVM_sum")

### Combine
dat_SVM <- rbind(datTrn_SVM, datTst_SVM)
dat_SVM$Status <- as.numeric(dat_SVM$Status)-1

### Replace levels
levels(dat_SVM$SVM_type) <- gsub("\\.", "-", levels(dat_SVM$SVM_type)) 
levels(dat_SVM$SVM_type) <- paste(levels(dat_SVM$SVM_type),forLegend_SVM[c(3,1,4,2)])

###################
### Plot Curve
png("./Figs/SVM/LOO/CompareAUCggplot_SVM_JHUOut.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic vs agnostic SVM (JHU for validation)"
legendTitle <- "Mechanistic VS Agnostic"
### Plot
basicplot_SVM_JHUOut <- ggplot(dat_SVM, aes(d=Status, m=SVM_sum, color=SVM_type,
                                                 linetype = SVM_type)) +
  geom_roc(cutoffs.at = seq(0.1,1,0.1)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=13)) +
  scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_SVM_JHUOut
### Close device
dev.off()

########################

save(basicplot_SVM_JHUOut, file = "./Objs/SVM/LOO/BasicPlot_SVM_JHUOut.rda")

