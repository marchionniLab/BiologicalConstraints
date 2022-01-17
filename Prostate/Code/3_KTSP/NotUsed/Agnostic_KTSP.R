################################################################################
### Mohamed Omar
### 10/4/2019
### Goal: Creating unrestricted K-TSP classifier
###############################################################################

rm(list = ls()) 

#setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

####################### 
##Load required packages
library(switchBox)
library(Biobase)
library(limma)
library(pROC)
library(caret)
library(RColorBrewer)
library(ggplot2)
library(reshape)
library(plotROC)
library(xtable)
library(mltools)

#####################################################################
### Load data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

##########################################################
## Quantile normalize the datasets
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

###########################################################################
### TRAINING using all expressed genes

## Set Feature number and max K
featN <- nrow(usedTrainMat) # the same as in the mechanistic classifier 
ktsp <- c(3:25)  # the same as in the mechanistic classifier


###########
### Train a classifier using the default filter function
set.seed(333)

ktspPredictorUnRes <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, 
                                      krange = ktsp, 
                                      FilterFunc = SWAP.Filter.Wilcoxon, featureNo= featN)
ktspPredictorUnRes

#Agnostic_KTSP <- cbind(ktspPredictorUnRes$TSPs, ktspPredictorUnRes$score)
#colnames(Agnostic_KTSP) <- c("gene1", "gene2", "score")

#print(xtable(Agnostic_KTSP, type = "latex"), file = "./Objs/KTSP/Agnostic.tex")

#####################################################################
## Compute the sum and find the best threshold: in the training samples
ktspStatsTrainUnRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTrainUnRes$statistics)

## Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

## Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: reorder the levels ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, plot = TRUE, ci = TRUE, print.auc = TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, main = "Agnostic KTSP ROC Train")

##############################################################
## Get prediction based on the best threshold from ROC curve
usedTrainPredictionUnRes <- SWAP.KTSP.Classify(inputMat = usedTrainMat, classifier = ktspPredictorUnRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionUnRes, usedTrainGroup, positive = "Mets")

MCC_Agnostic_Train <- mltools::mcc(pred = usedTrainPredictionUnRes, actuals = usedTrainGroup)
MCC_Agnostic_Train

#################################################################################
###############################################################################

#### Testing

# Compute the sum and find the best threshold
ktspStatsTestUnRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorUnRes, CombineFunc = sum)
summary(ktspStatsTestUnRes$statistics)

# plot the curve
roc(usedTestGroup, ktspStatsTestUnRes$statistics, plot = TRUE, print.thres = thr, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col= "blue", lwd=2, grid = TRUE, main = "Agnostic KTSP ROC Test")

### Get prediction based on best threshold from ROC curve
usedTestPredictionUnRes <- SWAP.KTSP.Classify(inputMat = usedTestMat, classifier = ktspPredictorUnRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionUnRes, usedTestGroup, positive = "Mets", mode = "everything")

MCC_Agnostic_Test <- mltools::mcc(pred = usedTestPredictionUnRes, actuals = usedTestGroup)
MCC_Agnostic_Test

### basic plots of various metrics
#library(precrec)
#precrec_obj <- evalmod(scores = ktspStatsTestUnRes$statistics, labels = usedTestGroup, mode="basic")
#autoplot(precrec_obj)   

#################################################################################
##################
### Plot genes in the training set
## Which TSPs
i <- 1:nrow(ktspPredictorUnRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ktspPredictorUnRes$TSPs)

## Assemble
dfTspTrain <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), out)
}, x=usedTrainMat, g=usedTrainGroup)

names(dfTspTrain) <- rownames(ktspPredictorUnRes$TSPs)

# Change the names of elements inside each element in dfTspTrain (For Plotting)  
for(i in seq_along(dfTspTrain)) names(dfTspTrain[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")


## Reduce
datTrain <- Reduce("rbind", dfTspTrain)

## Rename columns
#colnames(datTrain)[colnames(datTrain) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
# png("./Figs/KTSP/mechanistic.trainKTSPexprs.png", width = 3000, height = 1500, res = 200)
# bxplt <- ggplot(na.omit(datTrain), aes(x=Gene, y=Expression, fill=Group)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
#   facet_wrap(~pair, scales = "free", nrow = 2) +
#   theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
# bxplt
# dev.off()

#####################
## Scatter plots
png(filename = "./Figs/KTSP/AgnosticKTSP_Train_ScatterPlot.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datTrain), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 12),
        axis.title=element_text(face="bold", size = 12),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=17.5),
        legend.title = element_text(face="bold", size=17.5),
        strip.text.x = element_text(face="bold", size=11))
sctplt
dev.off()


########################################################################
#######################################################################
### Plot genes in the testing set

## Which TSPs
i <- 1:nrow(ktspPredictorUnRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ktspPredictorUnRes$TSPs)

## Assemble
dfTspTest <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), out)
}, x=usedTestMat, g=usedTestGroup)


names(dfTspTest) <- rownames(ktspPredictorUnRes$TSPs)

# Change the names of elements inside each element in dfTspTrain 
for(i in seq_along(dfTspTest)) names(dfTspTest[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")


## Reduce
datTest <- Reduce("rbind", dfTspTest)

## Rename columns
#colnames(datTest)[colnames(datTest) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
# png("./Figs/KTSP/mechanistic.testKTSPexprs.png", width = 3000, height = 1500, res = 200)
# bxplt <- ggplot(na.omit(datTest), aes(x=Gene, y=Expression, fill=Group)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
#   facet_wrap(~pair, scales = "free", nrow = 2) +
#   theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
# bxplt
# dev.off()

#####################
## Scatter plots
png(filename = "./Figs/KTSP/AgnosticKTSP_Test_ScatterPlot.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datTest), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 12),
        axis.title=element_text(face="bold", size = 12),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=17.5),
        legend.title = element_text(face="bold", size=17.5),
        strip.text.x = element_text(face="bold", size=11))
sctplt
dev.off()

