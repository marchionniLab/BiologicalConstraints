###### 
# Clean Work space
rm(list = ls())


setwd("/Volumes/Macintosh/Research/Projects/ICB_Response")

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

######################################################
## Load the Processed datasets
load("./Objs/ProcessedDatasets.rda")


#Get the Immune checkpoint genes
# CellActivImm <- read.delim("./CellActivImm.txt")
# CellActivImm <- as.character(CellActivImm[-1,])
# CellActivImm <- gsub("-.+", "", CellActivImm)


RegImmEffec <- read.delim("./RegImmEffec.txt")
RegImmEffec <- as.character(RegImmEffec[-1, ])
RegImmEffec <- gsub("-.+", "", RegImmEffec)


BcellAct <- read.delim("./BcellAct.txt")
BcellAct <- as.character(BcellAct[-1,])
BcellAct <- gsub("-.+","",BcellAct)

WNT_B_Cat <- read.delim("./WNT_B_Cat.txt")
WNT_B_Cat <- as.character(WNT_B_Cat[-1,])
WNT_B_Cat <- gsub("-.+", "", WNT_B_Cat)
 
# CytoProd <- read.delim("./CytoProd.txt")
# CytoProd <- as.character(CytoProd[-1,])
# CytoProd <- gsub("-.+","",CytoProd)

CaSignaling <- read.delim("./CaSignaling.txt")
CaSignaling <- as.character(CaSignaling[-1,])
CaSignaling <- gsub("-.+", "", CaSignaling)

# Genes <- read.delim("./geneset1.txt")
# Genes <- as.character(Genes[-1,])
# Genes <- gsub("-.+", "", Genes)

AllImmune <- c(RegImmEffec, BcellAct, WNT_B_Cat,CaSignaling)
AllImmune <- AllImmune[!duplicated(AllImmune)]

#myTSPs <- expand.grid(RegImmEffec, BcellAct)
#myTSPs <- as.matrix(myTSPs)
myTSPs <- t(combn(AllImmune, 2))

keepGns <- intersect(as.vector(myTSPs), rownames(Expr_NP))

myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

#keepGns2 <- c("PDCD1", "CD274", "CTLA4", "CD28", "CD80", "CD86")

#myTSPs <- myTSPs[myTSPs[,1] %in% keepGns2 | myTSPs[,2] %in% keepGns2 , ]

#########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) #13
featNo <- nrow(Expr_NP)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333) 

ktspPredictorRes <- SWAP.Train.KTSP(
  Expr_NP, Pheno_NP$SR, krange=20,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, 
  RestrictedPairs = myTSPs, disjoint = T)

ktspPredictorRes


#ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = Expr_NP, classifier = ktspPredictorRes, CombineFunc = sum)

#### make it numeric
#kMat <- t(1*ktspStatsTrainRes$comparisons)
#rownames(kMat) <- paste(rownames(kMat),  1:nrow(kMat), sep="_")

### Compute mutual info
#kMatMI <- mutualInfo(kMat)

### Plot
#plot(hclust(kMatMI))

# remove pairs N 3,5,6,7,10
X <- ktspPredictorRes$TSPs[-c(2,8,9,12,13,14,16,19), ]

ktspPredictorRes <- SWAP.Train.KTSP(
   Expr_NP, Pheno_NP$SR, krange=12,
   FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo,
   RestrictedPairs = X, disjoint = T)

ktspPredictorRes
 
############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = Expr_NP, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

KTSP_STATs_Train_Mechanistic <- t(ktspStatsTrainRes$comparisons)
KTSP_STATs_Train_Mechanistic[KTSP_STATs_Train_Mechanistic == FALSE] <- 0

### Threshold
thr <- coords(roc(Pheno_NP$SR, ktspStatsTrainRes$statistics, levels = c("no", "yes"), direction = "<"), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(Pheno_NP$SR, ktspStatsTrainRes$statistics, levels = c("no", "yes"), direction = "<" ), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(Pheno_NP$SR, ktspStatsTrainRes$statistics, plot = T, print.thres="all", print.auc=TRUE, print.auc.col="black", levels = c("no", "yes"), direction = "<", col="black", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(Expr_NP, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, Pheno_NP$SR, positive = "yes", mode = "everything")

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = Pheno_NP$SR)
MCC_Mechanistic_Train


###################################
###########################################################################
### Plot genes in the training set
## Which TSPs
# i <- 1:nrow(ktspPredictorRes$TSPs)
# 
# ## Assemble in a data frame
# tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ktspPredictorRes$TSPs)
# 
# ## Assemble
# dfTspTrain <- lapply(tsp, function(i,x,g){
#   out <- data.frame(t(x[i, ]), Group=g)
#   out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), out)
# }, x=Expr_NP, g=Pheno_NP$SR)
# 
# names(dfTspTrain) <- rownames(ktspPredictorRes$TSPs)
# 
# # Change the names of elements inside each element in dfTspTrain (For Plotting)  
# for(i in seq_along(dfTspTrain)) names(dfTspTrain[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")
# 
# 
# ## Reduce
# datTrain <- Reduce("rbind", dfTspTrain)

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
# png(filename = "./Figs/MechanisticKTSP_Train_ScatterPlot.png", width=3000, height=1500, res=300)
# ### Prepare ggplot
# sctplt <- ggplot(na.omit(datTrain), aes(x=Gene1, y=Gene2, color=Group)) +
#   geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
#   facet_wrap(~ pair, scales="free", nrow=2) +
#   theme(axis.text=element_text(face="bold", size = 12),
#         axis.title=element_text(face="bold", size = 12),
#         legend.position="bottom",
#         legend.text=element_text(face="bold", size=17.5),
#         legend.title = element_text(face="bold", size=17.5),
#         strip.text.x = element_text(face="bold", size=11))
# sctplt
# dev.off()


#########################################################################
#########################################################################
### Testing

# Test on Chen_Pre_aCTLA4

## Compute the sum and find the best threshold
ktspPredictorRes$labels <- c("NR", "R")

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_Chen_Pre_aCTLA4, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_Chen_Pre_aCTLA4$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_Chen_Pre_aCTLA4, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_Chen_Pre_aCTLA4$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_Chen_Pre_aCTLA4$Response)
MCC_Mechanistic_Test

###########################################3
# Test on Chen_Pre_aPD1

## Compute the sum and find the best threshold

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_Chen_Pre_aPD1, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_Chen_Pre_aPD1$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_Chen_Pre_aPD1, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_Chen_Pre_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_Chen_Pre_aPD1$Response)
MCC_Mechanistic_Test


###########################################3
# Test on Chen_On_aPD1

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_Chen_On_aPD1, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_Chen_On_aPD1$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_Chen_On_aPD1, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_Chen_On_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_Chen_On_aPD1$Response)
MCC_Mechanistic_Test

###########################################3
# Test on PRAT

ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = Expr_PRAT, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(Pheno_PRAT$Response, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Expr_PRAT, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_PRAT$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_PRAT$Response)
MCC_Mechanistic_Test


###########################################3
# Test on HUGO

ktspStatsHUGORes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_HUGO), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsHUGORes$statistics)

KTSP_STATs_HUGO_Mechanistic <- t(ktspStatsHUGORes$comparisons)
KTSP_STATs_HUGO_Mechanistic[KTSP_STATs_HUGO_Mechanistic == FALSE] <- 0

## Plot curve
RocHugo <- roc(Pheno_HUGO$Response, ktspStatsHUGORes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)
RocHugo
### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_HUGO), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_HUGO$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_HUGO$Response)
MCC_Mechanistic_Test

###########################################3
# Test on TCGA

ktspStatsTCGARes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_TCGA), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTCGARes$statistics)

KTSP_STATs_TCGA_Mechanistic <- t(ktspStatsTCGARes$comparisons)
KTSP_STATs_TCGA_Mechanistic[KTSP_STATs_TCGA_Mechanistic == FALSE] <- 0

## Plot curve
RocTcga <- roc(Pheno_TCGA$Response, ktspStatsTCGARes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)
RocTcga

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_TCGA), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_TCGA$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_TCGA$Response)
MCC_Mechanistic_Test

####################################################
# Test on RIAZ_Pre_aPD1

ktspStatsRIAZ_Pre_aPD1_Res <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_RIAZ_Pre_aPD1), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsRIAZ_Pre_aPD1_Res$statistics)

KTSP_STATs_RIAZ_Pre_aPD1_Mechanistic <- t(ktspStatsRIAZ_Pre_aPD1_Res$comparisons)
KTSP_STATs_RIAZ_Pre_aPD1_Mechanistic[KTSP_STATs_RIAZ_Pre_aPD1_Mechanistic == FALSE] <- 0

## Plot curve
RocRiaz_Pre_aPD1 <- roc(Pheno_RIAZ_Pre_aPD1$Response, ktspStatsRIAZ_Pre_aPD1_Res$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)
RocRiaz_Pre_aPD1

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_RIAZ_Pre_aPD1), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_RIAZ_Pre_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_RIAZ_Pre_aPD1$Response)
MCC_Mechanistic_Test

####################################################
# Test on RIAZ_On_aPD1

ktspStatsRIAZ_On_aPD1_Res <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_RIAZ_On_aPD1), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsRIAZ_On_aPD1_Res$statistics)

KTSP_STATs_RIAZ_On_aPD1_Mechanistic <- t(ktspStatsRIAZ_On_aPD1_Res$comparisons)
KTSP_STATs_RIAZ_On_aPD1_Mechanistic[KTSP_STATs_RIAZ_On_aPD1_Mechanistic == FALSE] <- 0

## Plot curve
RocRiaz_On_aPD1 <- roc(Pheno_RIAZ_On_aPD1$Response, ktspStatsRIAZ_On_aPD1_Res$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)
RocRiaz_On_aPD1

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_RIAZ_On_aPD1), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_RIAZ_On_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_RIAZ_On_aPD1$Response)
MCC_Mechanistic_Test

####################################################
# Test on MGH_On_aPD1

ktspStatsMGH_On_aPD1_Res <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_MGH_On_aPD1), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsMGH_On_aPD1_Res$statistics)

KTSP_STATs_MGH_On_aPD1_Mechanistic <- t(ktspStatsMGH_On_aPD1_Res$comparisons)
KTSP_STATs_MGH_On_aPD1_Mechanistic[KTSP_STATs_MGH_On_aPD1_Mechanistic == FALSE] <- 0

## Plot curve
RocMgh_On_aPD1 <- roc(Pheno_MGH_On_aPD1$Response, ktspStatsMGH_On_aPD1_Res$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)
RocMgh_On_aPD1

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_MGH_On_aPD1), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_MGH_On_aPD1$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_MGH_On_aPD1$Response)
MCC_Mechanistic_Test

####################################################
# Test on MGH_On_aCTLA4

ktspStatsMGH_On_aCTLA4_Res <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_MGH_On_aCTLA4), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsMGH_On_aCTLA4_Res$statistics)

KTSP_STATs_MGH_On_aCTLA4_Mechanistic <- t(ktspStatsMGH_On_aCTLA4_Res$comparisons)
KTSP_STATs_MGH_On_aCTLA4_Mechanistic[KTSP_STATs_MGH_On_aCTLA4_Mechanistic == FALSE] <- 0

## Plot curve
RocMgh_On_aCTLA4 <- roc(Pheno_MGH_On_aCTLA4$Response, ktspStatsMGH_On_aCTLA4_Res$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)
RocMgh_On_aCTLA4

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_MGH_On_aCTLA4), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_MGH_On_aCTLA4$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_MGH_On_aCTLA4$Response)
MCC_Mechanistic_Test

####################################################
# Test on VanAllen

ktspStatsVanAllenRes <- SWAP.KTSP.Statistics(inputMat = as.matrix(Expr_VanAllen), classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsVanAllenRes$statistics)

KTSP_STATs_VanAllen_Mechanistic <- t(ktspStatsVanAllenRes$comparisons)
KTSP_STATs_VanAllen_Mechanistic[KTSP_STATs_VanAllen_Mechanistic == FALSE] <- 0

## Plot curve
roc(Pheno_VanAllen$Response, ktspStatsVanAllenRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(as.matrix(Expr_VanAllen), ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Pheno_VanAllen$Response, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Pheno_VanAllen$Response)
MCC_Mechanistic_Test

####################################################
# Test on all Melanoma

ktspStatsMelanoma_Res <- SWAP.KTSP.Statistics(inputMat = Melanoma_Mat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsMelanoma_Res$statistics)

KTSP_STATs_Melanoma_Mechanistic <- t(ktspStatsMelanoma_Res$comparisons)
KTSP_STATs_Melanoma_Mechanistic[KTSP_STATs_Melanoma_Mechanistic == FALSE] <- 0

## Plot curve
RocMelanoma <- roc(Melanoma_Group, ktspStatsMelanoma_Res$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("NR", "R"), direction = "<", col="blue", lwd=2, grid=TRUE)
RocMelanoma

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(Melanoma_Mat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, Melanoma_Group, positive = "R", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = Melanoma_Group)
MCC_Mechanistic_Test

ktspPredictorRes
################3
### Plot genes in the testing set

## Which TSPs
# i <- 1:nrow(ktspPredictorRes$TSPs)
# 
# ## Assemble in a data frame
# tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ktspPredictorRes$TSPs)
# 
# ## Assemble
# dfTspTest <- lapply(tsp, function(i, x, g){
#   out <- data.frame(t(x[i, ]), Group=g)
#   out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), out)
# }, x=Melanoma_Mat, g=Melanoma_Group)
# 
# 
# names(dfTspTest) <- rownames(ktspPredictorRes$TSPs)
# 
# # Change the names of elements inside each element in dfTspTrain 
# for(i in seq_along(dfTspTest)) names(dfTspTest[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")
# 
# 
# ## Reduce
# datTest <- Reduce("rbind", dfTspTest)
# 
# #####################
# ## Scatter plots
# png(filename = "./Figs/MechanisticKTSP_Test_ScatterPlot.png", width=3000, height=1500, res=300)
# ### Prepare ggplot
# sctplt <- ggplot(na.omit(datTest), aes(x=Gene1, y=Gene2, color=Group)) +
#   geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
#   facet_wrap(~ pair, scales="free", nrow=2) +
#   theme(axis.text=element_text(face="bold", size = 12),
#         axis.title=element_text(face="bold", size = 12),
#         legend.position="bottom",
#         legend.text=element_text(face="bold", size=17.5),
#         legend.title = element_text(face="bold", size=17.5),
#         strip.text.x = element_text(face="bold", size=11))
# sctplt
# dev.off()

##################################
## Save statistics for each dataset

# save(KTSP_STATs_Train_Mechanistic, KTSP_STATs_HUGO_Mechanistic,
#      KTSP_STATs_TCGA_Mechanistic, KTSP_STATs_RIAZ_Pre_aPD1_Mechanistic,
#      KTSP_STATs_RIAZ_On_aPD1_Mechanistic, KTSP_STATs_MGH_On_aPD1_Mechanistic,
#      KTSP_STATs_VanAllen_Mechanistic,
#      KTSP_STATs_MGH_On_aCTLA4_Mechanistic, KTSP_STATs_Melanoma_Mechanistic,
#      file = "./Objs/KTSP_Stats.rda")






