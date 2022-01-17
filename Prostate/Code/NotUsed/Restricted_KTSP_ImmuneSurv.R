###############################################################################
### Mohamed Omar
### 01/07/2019
### Goal : Creating the restricted K-TSP classifier for prostate cancer metastasis
### Using : TF genes
#################################################################################

rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

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

###################################################################
## Load Metastasis Data
load("./Objs/MetastasisData.rda")
load("./mytspsCombined.rdata")

## Load TF_MiR genes
immuneSurv <- load("/Users/mohamedomar/Desktop/ImmuneSurvGenes.rda")

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(mixTrainMat))

## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[keepGns, ]
boxplot(usedTrainMat)

usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[keepGns, ]
boxplot(usedTestMat)

#
usedTrainGroup <- mixTrainGroup
names(usedTrainGroup) <- colnames(usedTrainMat)
usedTestGroup <- mixTestGroup
names(usedTestGroup) <- colnames(usedTestMat)


### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

### Set Feature number and max k
featN <- nrow(usedTrainMat) #nrow(usedTrainMat)
ktsp <- 9 #22 #16

### Train a classifier using default filtering function based on Wilcoxon
ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featN, RestrictedPairs = myTSPs)
ktspPredictorRes

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = ">",), transpose = TRUE, "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = ">", ), transpose = TRUE ,"local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = TRUE, print.thres=thr, print.thres.adj=c(0.01,1.25), print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = ">", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Mets")

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Threshold
thr_test <- coords(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No_Mets", "Mets"), direction = ">",),"best", transpose = FALSE)["threshold"]
thr_test

## Print ROC curve local maximas
coords(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No_Mets", "Mets"), direction = ">",), transpose = TRUE, "local maximas")

## Plot curve
png("./Figs/ImmuneSurv/ROC_Test_TF_MiR.png", width = 2000, height = 2000, res = 300)
roc(usedTestGroup, ktspStatsTestRes$statistics, plot = TRUE, print.thres=thr_test, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = ">", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")
dev.off()

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr_test)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Mets")

## OR (calculate accuracy, sens, spec and AUC)
SWAP.GetKTSP.PredictionStats(usedTestPredictionRes, usedTestGroup, decision_values = ktspStatsTestRes$statistics)

############################################################################
###########################################################################
### Plot genes in the training set

## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ktspPredictorRes$TSPs)

## Assemble
dfTspTrain <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=usedTrainMat, g=usedTrainGroup)

## Reduce
datTrain <- Reduce("rbind", dfTspTrain)

## Rename columns
colnames(datTrain)[colnames(datTrain) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
png("./Figs/TF/mechanistic.trainKTSPexprs.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTrain), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(group=Group), position = position_jitterdodge(1.2), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

########################################################################
#######################################################################
### Plot genes in the testing set

## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ktspPredictorRes$TSPs)

## Assemble
dfTspTest <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), melt(out))
}, x=usedTestMat, g=usedTestGroup)

## Reduce
datTest <- Reduce("rbind", dfTspTest)

## Rename columns
colnames(datTest)[colnames(datTest) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
png("./Figs/TF/mechanistic.testKTSPexprs.png", width = 3000, height = 1500, res = 200)
bxplt <- ggplot(na.omit(datTest), aes(x=Gene, y=Expression, fill=Group)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(aes(group=Group), position = position_jitterdodge(1.2), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
  facet_wrap(~pair, scales = "free", nrow = 2) +
  theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
bxplt
dev.off()

#######################################################################
#### SAVE
save(list = ls(pattern = "^ktsp"), file = "./Objs/TF/mechanistic.ktspPredictor.rda")


#################################################################################
##################################################################################
##################################################################################

