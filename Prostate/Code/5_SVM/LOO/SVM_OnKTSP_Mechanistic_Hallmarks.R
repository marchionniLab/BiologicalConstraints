########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: SVM (Poly) for predicting bladder cancer progression
## Mechanistic (Using TF-MiR genes)

###########################################################################

## Clean the work environment
rm(list = ls())

## Set the working directory
#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(mltools)

#######################################################################

## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Hallmarks.rda")
load("./Objs/MetastasisDataGood.rda")



### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

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
Grid <- expand.grid(degree = 3, scale = 0.1, C = 0.25)
set.seed(333)
fit.svmPoly_mechanistic_OnKTSP <- train(usedTrainGroup~., data=Data_train, method="svmPoly", trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), tuneGrid = Grid, metric = c("ROC"))
fit.svmPoly_mechanistic_OnKTSP

##########################################################
## Predict in the training data
train_pred_classes_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Training, type="raw")
table(train_pred_classes_svmPoly_mechanistic_OnKTSP)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_train_svmPoly_OnKTSP <- confusionMatrix(train_pred_classes_svmPoly_mechanistic_OnKTSP, usedTrainGroup, positive = "Mets")
Confusion_train_svmPoly_OnKTSP

MCC_Train <- mltools::mcc(pred = train_pred_classes_svmPoly_mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Train

#############################################################
## Predict in the testing data 

test_pred_classes_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Testing, type="raw")
table(test_pred_classes_svmPoly_mechanistic_OnKTSP)

## Creat a confusion matrix (the predictions against the true classes)
Confusion_test_svmPoly_OnKTSP <- confusionMatrix(test_pred_classes_svmPoly_mechanistic_OnKTSP, usedTestGroup, positive = "Mets")
Confusion_test_svmPoly_OnKTSP

MCC_Mechanistic <- mltools::mcc(pred = test_pred_classes_svmPoly_mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Mechanistic

## ROC/AUC in the Testing set
test_pred_prob_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Testing, type = "prob")
roc(usedTestGroup, test_pred_prob_svmPoly_mechanistic_OnKTSP[,2], plot = F, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Mechanistic)")

## Top 100 predictors
Importance_SVMPoly_Mechanistic <- varImp(fit.svmPoly_mechanistic_OnKTSP, scale = FALSE)
png("./Figs/SVM/SVM_varImp_Mechanistic_OnKTSP.png", width = 2000, height = 2000, res = 200)
plot(Importance_SVMPoly_Mechanistic, top = 20, main = "Mechanistic SVM top features")
dev.off()

Importance_SVMPoly_Mechanistic <- Importance_SVMPoly_Mechanistic$importance
Importance_SVMPoly_Mechanistic <- Importance_SVMPoly_Mechanistic[order(Importance_SVMPoly_Mechanistic$Progression, decreasing = TRUE),]
MechanisticSVM_Pairs <- rownames(Importance_SVMPoly_Mechanistic)[1:50]
save(MechanisticSVM_Pairs, file = "./Objs/SVM/MechanisticSVM_Pairs.rda")

########################
## ROC stat for training data
train_pred_prob_svmPoly_mechanistic_OnKTSP <- predict(fit.svmPoly_mechanistic_OnKTSP, Training, type = "prob")
roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2], plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main = "ROC Test SVM Poly (Mechanistic)")

############################################################################
## Plot ggplot2 ROC curves for both the agnostic and mechanistic XGBoost 

### Prepare the legend
forLegend_SVM <- apply(rbind(
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2],levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


####################
### ROC curves Using ggplot2

### Training
datTrn_SVM <- melt(data.frame(
  ## Training Group
  Training=factor(usedTrainGroup, levels=c("No_Mets", "Mets")),
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
  Testing=factor(usedTestGroup, levels=c("No_Mets", "Mets")),
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
png("./Figs/SVM/CompareAUCggplot_SVM_Hallmarks.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic vs agnostic SVM"
legendTitle <- "Mechanistic VS Agnostic"
### Plot
basicplot_SVM <- ggplot(dat_SVM, aes(d=Status, m=SVM_sum, color=SVM_type,
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
basicplot_SVM
### Close device
dev.off()

########
save(basicplot_SVM, file = "./Objs/SVM/BasicPlot_SVM_Hallmarks.rda")

###########################
# Load the pre-saved ROC stats for all the RF models (Mechanistic and agnostic with different N of pairs)
load("./Objs/SVM/SVMAUCStats.rda")

# Load the original data just for the phenotype variable
load("./Objs/PhenoStrSampling.rda")

## For the forest plot
AUCs_SVM_StrSampling <-  rbind(
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP_200[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_agnostic_OnKTSP_200[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_prob_svmPoly_agnostic_OnKTSP_500[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, test_pred_prob_svmPoly_agnostic_OnKTSP_500[,2], levels = c("No_Mets", "Mets"), direction = "<"))
)


AUCs_SVM_StrSampling[,c(1,2,3)] <- AUCs_SVM_StrSampling[,c(2,1,3)]
colnames(AUCs_SVM_StrSampling) <- c("AUC", "CI_low", "CI_High")



AUCs_SVM_StrSampling <- as.data.frame(AUCs_SVM_StrSampling)

AUCs_SVM_StrSampling$model_type <- c("Mechanistic", "Mechanistic", "Agnostic", "Agnostic","Agnostic200", "Agnostic200", "Agnostic500", "Agnostic500")
AUCs_SVM_StrSampling$data_type <- c("Training", "Testing", "Training", "Testing", "Training", "Testing", "Training", "Testing")
AUCs_SVM_StrSampling$algorithm <- rep("SVM", 8)
AUCs_SVM_StrSampling$approach <- rep("Str.Sampling", 8)  

save(AUCs_SVM_StrSampling, file = "./Objs/SVM/AUCs_SVM_StrSampling.rda")



