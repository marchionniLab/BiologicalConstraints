########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: RF for predicting Prostate cancer Metastasis
# Mechanistic: gene pairs
# Corss validation: JHU Out

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
library(varSelRF)


## Load data
load("./Objs/KTSP/LOO/KTSP_STATs_Mechanistic_JHUOut.rda")
load("./Objs/LOO/MetastasisData_JHUOut.rda")


usedTrainGroup <- trainGroup
usedTestGroup <- testGroup


predictor_data <- t(KTSP_STATs_Train_Mechanistic)

tmp <- as.vector(table(usedTrainGroup))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)



set.seed(333)
tuneRF(x=predictor_data, y=usedTrainGroup, plot = TRUE, improve = 0.01, ntreeTry = 501, proximity = TRUE, sampsize = sampsizes, na.action = na.omit)

set.seed(333)
RF_Mechanistic_OnKTSP <- randomForest(x =predictor_data, y=usedTrainGroup, importance = TRUE, ntree = 501, mtry = 20 ,proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
RF_Mechanistic_OnKTSP
plot(RF_Mechanistic_OnKTSP)


#VAR <- varSelRF(xdata = predictor_data, Class = usedTrainGroup, ntree = 501, ntreeIterat = 20, whole.range = FALSE, keep.forest = TRUE)

#SelectedVars <- VAR$selected.vars

## RandomForest calculates an importance measures for each variable.
# rf_importances <- randomForest::importance(RF_Mechanistic_OnKTSP, scale=FALSE)
# rf_importances <- rf_importances[order(rf_importances[,4], decreasing = TRUE), ]
# ImportantVariables <- rownames(rf_importances)[1:50]
# 
# varImpPlot(RF_Mechanistic_OnKTSP, type=2, n.var=30, scale=FALSE, main="Variable Importance (Gini) for top 30 predictors")

# ROC curve in training data
train_pred_votes_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, newdata = predictor_data, type = "vote")
roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)

# Confusion Matrix
train_pred_responses_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data, type="response")
confusion_train <- confusionMatrix(train_pred_responses_Mechanistic_OnKTSP, usedTrainGroup, positive = "Mets")
confusion_train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = train_pred_responses_Mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Train


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

# ROC curve in the testing data
roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)

### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "Mets")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test




############################################################################
## Plot ggplot2 ROC curves for both the agnostic and mechanistic XGBoost 

### Prepare the legend
forLegend_RF <- apply(rbind(
  ci(roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2])),
  ci(roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2])),
  ci(roc(usedTrainGroup, train_pred_votes_Agnostic_OnKTSP[,2])),
  ci(roc(usedTestGroup, RF_predictions_votes_Agnostic_OnKTSP[,2]))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


####################
### ROC curves Using ggplot2

### Training
datTrn_RF <- melt(data.frame(
  ## Training Group
  Training= usedTrainGroup,
  Agnostic.Training.Pairs = train_pred_votes_Agnostic_OnKTSP[,2],
  Mechanistic.Training.Pairs = train_pred_votes_Mechanistic_OnKTSP[,2]
))

### Change Colnames
colnames(datTrn_RF) <- c("Status", "RF_type", "RF_sum")


### Testing
datTst_RF <- melt(data.frame(
  ## Testing group
  Testing= usedTestGroup, 
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Testing.Pairs = RF_predictions_votes_Agnostic_OnKTSP[,2],
  Mechanistic.Testing.Pairs = RF_predictions_votes_Mechanistic_OnKTSP[,2]
))
### Change Colnames
colnames(datTst_RF) <- c("Status", "RF_type", "RF_sum")

### Combine
dat_RF <- rbind(datTrn_RF, datTst_RF)
dat_RF$Status <- as.numeric(dat_RF$Status)-1

### Replace levels
levels(dat_RF$RF_type) <- gsub("\\.", "-", levels(dat_RF$RF_type))
levels(dat_RF$RF_type) <- paste(levels(dat_RF$RF_type), forLegend_RF[c(3,1,4,2)])

###################
### Plot Curve
png("./Figs/RF/LOO/CompareAUCggplot_JHUOut.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic vs agnostic RF (JHU for validation)"
legendTitle <- "Mechanistic VS Agnostic"
### Plot
basicplot_RF_JHUOut <- ggplot(dat_RF, aes(d=Status, m=RF_sum, color=RF_type,
                                               linetype = RF_type)) +
  geom_roc(cutoffs.at = seq(0.1,1,0.1)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="bold", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_RF_JHUOut
### Close device
dev.off()


save(basicplot_RF_JHUOut, file = "./Objs/RF/LOO/BasicPlot_RF_JHUOut.rda")

