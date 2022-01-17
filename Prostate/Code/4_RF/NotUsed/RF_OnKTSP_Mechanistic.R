
########################################################################### 
## Mohamed Omar
## 03/06/2019
## Goal: RF for predicting bladder cancer progression


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
library(randomForest)
library(varSelRF)


## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Combined.rda")
load("./Objs/MetastasisDataGood.rda")


usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#save(usedTrainGroup, usedTestGroup, file = "./Objs/PhenoStrSampling.rda")

predictor_data <- t(KTSP_STATs_Train_Mechanistic)

tmp <- as.vector(table(usedTrainGroup))
num_classes <- length(tmp)
min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
sampsizes <- rep(min_size,num_classes)



set.seed(333)
tuneRF(x=predictor_data, y=usedTrainGroup, plot = TRUE, improve = 0.01, ntreeTry = 1001, proximity = TRUE, sampsize = sampsizes, na.action = na.omit)

set.seed(333)
RF_Mechanistic_OnKTSP <- randomForest(x =predictor_data, y=usedTrainGroup, importance = TRUE, ntree = 1001, mtry = 10, proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
RF_Mechanistic_OnKTSP
plot(RF_Mechanistic_OnKTSP)


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
roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)

# predictions in the training data
train_pred_response_Mechanistic_OnKTSP <- predict(RF_Mechanistic_OnKTSP, predictor_data, type="response")

### Predict in the testing data
confusion_train <- confusionMatrix(train_pred_response_Mechanistic_OnKTSP, usedTrainGroup, positive = "Mets")
confusion_train

# Calculate Matthews correlation coefficient
MCC_Train <- mcc(preds = train_pred_response_Mechanistic_OnKTSP, actuals = usedTrainGroup)
MCC_Train

################################################################################ 
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


### Predict in the testing data
confusion_test <- confusionMatrix(RF_predictions_responses_Mechanistic_OnKTSP, usedTestGroup, positive = "Mets")
confusion_test

# Calculate Matthews correlation coefficient
MCC_Test <- mcc(preds = RF_predictions_responses_Mechanistic_OnKTSP, actuals = usedTestGroup)
MCC_Test

roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)



############################################################################
## Plot ggplot2 ROC curves for both the agnostic and mechanistic XGBoost 

### Prepare the legend
forLegend_RF <- apply(rbind(
  ci(roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_votes_Agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


####################
### ROC curves Using ggplot2

### Training
datTrn_RF <- melt(data.frame(
  ## Training Group
  Training=factor(usedTrainGroup, levels=c("No_Mets", "Mets")),
  Agnostic.Training.Pairs = train_pred_votes_Agnostic_OnKTSP[,2],
  Mechanistic.Training.Pairs = train_pred_votes_Mechanistic_OnKTSP[,2]
  ))

### Change Colnames
colnames(datTrn_RF) <- c("Status", "RF_type", "RF_sum")


### Testing
datTst_RF <- melt(data.frame(
  ## Testing group
  Testing=factor(usedTestGroup, levels=c("No_Mets", "Mets")),
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
png("./Figs/RF/CompareAUCggplot.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic vs agnostic RF (Pairs))"
legendTitle <- "Mechanistic VS Agnostic"
### Plot
basicplot_RF_Pairs <- ggplot(dat_RF, aes(d=Status, m=RF_sum, color=RF_type,
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
basicplot_RF_Pairs
### Close device
dev.off()


save(basicplot_RF_Pairs, file = "./Objs/RF/BasicPlot_RF_Pairs.rda")


####################
# Load the pre-saved ROC stats for all the RF models (Mechanistic and agnostic with different N of pairs)
load("./Objs/RF/RFAUCStats.rda")

# Load the original data just for the phenotype variable
load("./Objs/PhenoStrSampling.rda")

## For the forest plot
AUCs_RF_StrSampling <-  rbind(
  ci(roc(usedTrainGroup, train_pred_votes_Mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Mechanistic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_votes_Agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Agnostic_OnKTSP[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_votes_Agnostic_OnKTSP_200[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Agnostic_OnKTSP_200[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, train_pred_votes_Agnostic_OnKTSP_500[,2], levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, RF_predictions_votes_Agnostic_OnKTSP_500[,2], levels = c("No_Mets", "Mets"), direction = "<"))
)


AUCs_RF_StrSampling[,c(1,2,3)] <- AUCs_RF_StrSampling[,c(2,1,3)]
colnames(AUCs_RF_StrSampling) <- c("AUC", "CI_low", "CI_High")



AUCs_RF_StrSampling <- as.data.frame(AUCs_RF_StrSampling)

AUCs_RF_StrSampling$model_type <- c("Mechanistic", "Mechanistic", "Agnostic", "Agnostic", "Agnostic200", "Agnostic200", "Agnostic500", "Agnostic500")
AUCs_RF_StrSampling$data_type <- c("Training", "Testing", "Training", "Testing", "Training", "Testing", "Training", "Testing")
AUCs_RF_StrSampling$algorithm <- rep("RF", 8)
AUCs_RF_StrSampling$approach <- rep("Str.Sampling", 8)  

save(AUCs_RF_StrSampling, file = "./Objs/RF/AUCs_RF_StrSampling.rda")

