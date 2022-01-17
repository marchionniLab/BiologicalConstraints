###################################################################################
### Mohamed Omar
### 5/5/2019
### ### Goal: Creating the restricted ktsp classifier.
### Combined Adhesion and activation genes
#################################################################################

###### 
# Clean Work space
rm(list = ls())

# Set work directory
#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(enrichR)
library(mltools)
library(xtable)

###########################################################################
### Load expression and phenotype data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

############################################################################
# load mechanistic pairs
load("./Objs/StromalPairs_C12C1.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Common genes
myTSPs <- as.matrix(StromalPairs_C12C1)
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))
#keepGns_TF_MiR <- keepGns
#save(keepGns_TF_MiR, file = "./Objs/KTSP/KeepGns_TF_MiR.rda")

#usedTrainMat <- usedTrainMat
#usedTestMat <- usedTestMat

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

#print(xtable(myTSPs, type = "latex"), file = "./Objs/KTSP/Restricted_Pairs.tex")
#write.csv(myTSPs, file = "./Objs/KTSP/Restricted_Pairs.csv")

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #18
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange = ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

ktspPredictorRes

###########################################################################
### Check consistency with biology
keep <- ktspPredictorRes$TSPs[,1] %in% myTSPs[,"c12"] & ktspPredictorRes$TSPs[,2] %in% myTSPs[,"c1"]
table(keep)

#Subset
ktspPredictorRes$name <- paste(sum(keep), "TSPs", sep="")
ktspPredictorRes$TSPs <- ktspPredictorRes$TSPs[ keep, ]
ktspPredictorRes$score <- ktspPredictorRes$score[ keep ]
ktspPredictorRes$tieVote <- droplevels(ktspPredictorRes$tieVote[keep])

# Visualize the classifier
ktspPredictorRes


############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

x = as.data.frame(ktspStatsTrainRes$comparisons)
x$label = usedTrainGroup

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, ci = T, print.thres=thr, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")
ROCTrain

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
ConfusionTrain <- confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Mets", mode = "everything")
ConfusionTrain

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = usedTrainGroup)
MCC_Mechanistic_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, ConfusionTrain$overall["Accuracy"], ConfusionTrain$byClass["Balanced Accuracy"], ConfusionTrain$byClass["Sensitivity"], ConfusionTrain$byClass["Specificity"], MCC_Mechanistic_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
ROCTest <- roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
ROCTest

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
ConfusionTest <- confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Mets", mode = "everything")
ConfusionTest

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Mechanistic_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, ConfusionTest$overall["Accuracy"], ConfusionTest$byClass["Balanced Accuracy"], ConfusionTest$byClass["Sensitivity"], ConfusionTest$byClass["Specificity"], MCC_Mechanistic_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
MechanisticKTSP_Perf <- cbind(TrainPerf, TestPerf)

# Save
save(MechanisticKTSP_Perf, file = "./Objs/KTSP/MechanisticKTSP_Perf.rda")


## GGPLOT to compare the two classifiers
### Prepare the legend
forLegend_KTSP <- apply(rbind(
  ci(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


#################################################################
### ROC curves Using ggplot2

### Training
datTrn_KTSP <- melt(data.frame(
  ## Training Group
  Training= usedTrainGroup,
  ## Mechanistic KTSP SUM training
  Training_data= ktspStatsTrainRes$statistics))
### Change Colnames
colnames(datTrn_KTSP) <- c("Status", "data_type", "KTSP_sum")


### Testing
datTst_KTSP <- melt(data.frame(
  ## Testing group
  Testing= usedTestGroup,
  ## Mechanistic KTSP SUM training
  Testing_data=ktspStatsTestRes$statistics))
### Change Colnames
colnames(datTst_KTSP) <- c("Status", "data_type", "KTSP_sum")

### Combine
dat_KTSP <- rbind(datTrn_KTSP, datTst_KTSP)
dat_KTSP$Status <- as.numeric(dat_KTSP$Status)-1
####
### Replace levels
levels(dat_KTSP$data_type) <- gsub("\\.", "-", levels(dat_KTSP$data_type))
levels(dat_KTSP$data_type) <- paste(levels(dat_KTSP$data_type), forLegend_KTSP[c(1,2)])

#################################################################
### Plot Curve
png("./Figs/KTSP/AUCggplot_SingleCellStroma_C12C1.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "Performance of the 7-TSPs model at predicting prostate cancer metastasis"
legendTitle <- paste("Mechanistic (", nrow(ktspPredictorRes$TSPs), " pairs)", "C1-C12")
### Plot
basicplot_KTSP_SingleCellStroma <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=data_type,
                                                        linetype = data_type)) +
  geom_roc(cutoffs.at = seq(1,20,2)) +
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
basicplot_KTSP_SingleCellStroma
### Close device
dev.off()

save(basicplot_KTSP_SingleCellStroma, file = "./Objs/KTSP/basicplot_KTSP_SingleCellStroma.rda")




