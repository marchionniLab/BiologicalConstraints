###################################################################################
### Mohamed Omar
### 25/11/2019
### ### Goal: Creating the restricted ktsp classifier.
### TF_targets genes
## Leaving 1 dataset out (GSE41408)
#################################################################################

###### 
# Clean Work space
rm(list = ls())
# Set work directory
setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

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
load("./Objs/LOO/MetastasisData_GSE41408Out.rda")
load("./Objs/Correlation/RGenes.rda")

############################################################################
## Load the selected genes
## Load the selected genes
Genes1 <- read.delim("./GO_Adhesion.txt")
Genes1 <- as.matrix(Genes1)
Genes1 <- Genes1[-1,]

Genes2 <- read.delim("./GO_Activation.txt")
Genes2 <- as.matrix(Genes2)
Genes2 <- Genes2[-1,]

Genes3 <- read.delim("./GO_O2Response.txt")
Genes3 <- as.matrix(Genes3)
Genes3 <- Genes3[-1,]

Genes <- c(Genes1,Genes2, Genes3)
Genes <- Genes[!duplicated(Genes)]

myTSPs <- t(combn(Genes,2))

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")[RGenes,]
usedTestMat <- testMat[RGenes, ]

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))

usedTrainMat <- usedTrainMat[keepGns, ]
usedTestMat <- usedTestMat[keepGns, ]

### Associated groups
usedTrainGroup <- trainGroup
usedTestGroup <- testGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #7
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

ktspPredictorRes

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<" ), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
ROCTrain <- roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, ci = T, print.thres.adj=c(0.01,1.25), print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")
ROCTrain

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionTrain <- confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Mets", mode = "everything")
confusionTrain

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = usedTrainGroup)
MCC_Mechanistic_Train

# Put the performance metrics together
TrainPerf <- data.frame("Training" = c(ROCTrain$ci, confusionTrain$overall["Accuracy"], confusionTrain$byClass["Balanced Accuracy"], confusionTrain$byClass["Sensitivity"], confusionTrain$byClass["Specificity"], MCC_Mechanistic_Train))
TrainPerf[1:3, ] <- TrainPerf[c(2,1,3), ]
rownames(TrainPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
ROCTest <- roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, ci = T, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_targets Gns")
ROCTest

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionTest <- confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Mets", mode = "everything")
confusionTest

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Mechanistic_Test

## Group the performance metrics together
TestPerf <- data.frame("Testing" = c(ROCTest$ci, confusionTest$overall["Accuracy"], confusionTest$byClass["Balanced Accuracy"], confusionTest$byClass["Sensitivity"], confusionTest$byClass["Specificity"], MCC_Mechanistic_Test))
TestPerf[1:3, ] <- TestPerf[c(2,1,3), ]
rownames(TestPerf) <- c("AUC", "AUC_CI_low", "AUC_CI_high", "Accuracy", "Bal.Accuracy", "Sensitivity", "Specificity", "MCC")

## Group the performance metrics of the classifier in one data frame
GSE41408_Out_MechPerformance <- cbind(TrainPerf, TestPerf)

# Save
save(GSE41408_Out_MechPerformance, file = "./Objs/KTSP/GSE41408_Out_MechPerformance.rda")

#########################################################################
######################################################

## GGPLOT to compare the two classifiers
### Prepare the legend
forLegend_KTSP <- apply(rbind(
  ci(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"))
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
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Training = ktspStatsTrainUnRes$statistics,
  ## Mechanistic KTSP SUM training
  Mechanistic.Training= ktspStatsTrainRes$statistics))
### Change Colnames
colnames(datTrn_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")


### Testing
datTst_KTSP <- melt(data.frame(
  ## Testing group
  Testing= usedTestGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Testing=ktspStatsTestUnRes$statistics,
  ## Mechanistic KTSP SUM training
  Mechanistic.Testing=ktspStatsTestRes$statistics))
### Change Colnames
colnames(datTst_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")

### Combine
dat_KTSP <- rbind(datTrn_KTSP, datTst_KTSP)
dat_KTSP$Status <- as.numeric(dat_KTSP$Status)-1
####
### Replace levels
levels(dat_KTSP$KTSP_type) <- gsub("\\.", "-", levels(dat_KTSP$KTSP_type))
levels(dat_KTSP$KTSP_type) <- paste(levels(dat_KTSP$KTSP_type), forLegend_KTSP[c(3,1,4,2)])

#################################################################
### Plot Curve
png("./Figs/KTSP/LOO/CompareAUCggplot_GSE41408Out.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic (TFs_Targets) vs agnostic K-TSP"
legendTitle <- paste("Mechanistic (", nrow(ktspPredictorRes$TSPs), " pairs)",
                     " Agnostic (", nrow(ktspPredictorUnRes$TSPs), " pairs)",  sep="")
### Plot
basicplot_KTSP_GSE41408Out <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=KTSP_type,
                                                   linetype = KTSP_type)) +
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
basicplot_KTSP_GSE41408Out
### Close device
dev.off()

save(basicplot_KTSP_GSE41408Out, file = "./Objs/KTSP/LOO/BasicPlot_KTSP_GSE41408Out.rda")

##############################################################
