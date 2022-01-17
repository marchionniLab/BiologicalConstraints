###################################################################################
### Mohamed Omar
### 29/6/2019
### ### Goal: Creating a classifier that can classify prostate cancer metastasis into Bone // Visceral metastasis
### Using: K-TSPs
### TF_MIR genes (as restricted pairs)
#################################################################################

###### 
# Clean Work space
rm(list = ls())

# Set work directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

############################################################################
### Load library
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(GEOquery)

############################################################################
## Download the data sets

Dataset1 <- getGEO("GSE101607", GSEMatrix = TRUE, AnnotGPL = TRUE)
Dataset1 <- Dataset1$GSE101607_series_matrix.txt.gz

## part of Dataset3
#Dataset2 <- getGEO("GSE74367", GSEMatrix = TRUE, AnnotGPL = TRUE)
#Dataset2 <- Dataset2$GSE74367_series_matrix.txt.gz

Dataset3 <- getGEO("GSE74685", GSEMatrix = TRUE, AnnotGPL = TRUE)
Dataset3 <- Dataset3$GSE74685_series_matrix.txt.gz


############################################################

## Get the expression matrices

expr1 <- exprs(Dataset1)
expr3 <- exprs(Dataset3)

## Get feature data
FeatData1 <- fData(Dataset1)
FeatData3 <- fData(Dataset3)

## Get phenotype
pheno1 <- pData(Dataset1)
pheno3 <- pData(Dataset3)

#############################################################

## Annotate exprs1
# Log2 transform
expr1 <- log2(expr1)
rownames(expr1) <- FeatData1$`Gene symbol`
summary(is.na(rownames(expr1)))
rownames(expr1) <- gsub("-","", rownames(expr1))
rownames(expr1) <- gsub("_","",rownames(expr1))
sel <- which(apply(expr1, 1, function(x) all(is.finite(x)) ))
expr1 <- expr1[sel, ]
expr1 <- expr1[!(rownames(expr1) == ""), ]

## Annotate expr2
rownames(expr2) <- FeatData2$GENE_SYMBOL

## Annotate expr3
rownames(expr3) <- FeatData3$GENE_SYMBOL
summary(is.na(rownames(expr3)))
rownames(expr3) <- gsub("-","", rownames(expr3))
rownames(expr3) <- gsub("_","",rownames(expr3))
sel <- which(apply(expr3, 1, function(x) all(is.finite(x)) ))
expr3 <- expr3[sel, ]
expr3 <- expr3[!(rownames(expr3) == ""), ]
dim(expr3)


#rownames(expr4) <- FeatData4$GENE_SYMBOL

#########################################

## Modify pheno1 

# Keep only prostate samples
pheno1 <- pheno1[grep("prostate", pheno1$`disease state:ch1`), ]
# Remove untreated prostate cancer (keep only castration-resistant)
pheno1 <- pheno1[grep("Castration", pheno1$`disease state:ch1`), ]

pheno1$Mets <- as.factor(pheno1$`tissue:ch1`)
levels(pheno1$Mets) <- c("Bone_Mets")

## Modify expr1
expr1 <- expr1[,colnames(expr1) %in% rownames(pheno1)]


## Replace both expr1 and pheno1 in Dataset1

exprs(Dataset1) <- expr1

############################################################
## Modify pheno2

# Remove xenograft samples
pheno2 <- pheno2[-c(1:20), ]

# Remove Primary tumors
pheno2 <- pheno2[-c(46:56), ]

pheno2$Mets <- pheno2$title

pheno2$Mets <- gsub(".+_", "", pheno2$Mets)

pheno2$Mets[pheno2$Mets == "LIVER CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "LN CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "RETROPERITONEAL CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "LUNG CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "KIDNEY CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "APPENDIX CRPC metastasis"] <- "Visceral_Mets"
pheno2$Mets[pheno2$Mets == "PERITONEUM CRPC metastasis"] <- "Visceral_Mets"

pheno2$Mets[pheno2$Mets == "BONE CRPC metastasis"] <- "Bone_Mets"

pheno2$Mets <- as.factor(pheno2$Mets)
levels(pheno2$Mets)
table(pheno2$Mets)

## Modify expr2
expr2 <- expr2[, colnames(expr2) %in% rownames(pheno2)]


###################################################################
## Modify pheno3
pheno3$Mets <- pheno3$title
pheno3$Mets <- gsub(".+_", "", pheno3$Mets)

pheno3$Mets[pheno3$Mets == "LIVER CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "LN CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "LUNG CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "RETROPERITONEAL CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "RENAL CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "ADRENAL CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "KIDNEY CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "APPENDIX CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "PERITONEUM CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "PERITONEAL CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "SCROTUM CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "SKIN CRPC metastasis"] <- "Visceral_Mets"
pheno3$Mets[pheno3$Mets == "SPLEEN CRPC metastasis"] <- "Visceral_Mets"

pheno3$Mets[pheno3$Mets == "BONE CRPC metastasis"] <- "Bone_Mets"

pheno3$Mets <- as.factor(pheno3$Mets)

#############################################################
#############################################################

## Combine the expression and phenotype

AllExprs <- list(expr1, expr3)
names(AllExprs) <- c("GSE101607", "GSE74685")

AllPheno <- list(pheno1, pheno3)
names(AllPheno) <- c("GSE101607", "GSE74685")

GroupMets <- c(pheno1$BoneMets, pheno3$Mets)
GroupMets <- as.factor(GroupMets)
levels(GroupMets) <- c("Bone_Mets", "Visceral_Mets")
table(GroupMets)

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(AllExprs, rownames))

### Filter expression for the common genes
exprsMetastasis <- mapply(x=AllExprs, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

##########
## Assemble in one data frame
AllMat <- do.call("cbind", exprsMetastasis)
## Check if sample names are identical
all(colnames(AllMat) == names(GroupMets))

## Divide data into training and testing data
ind <- createDataPartition(y = GroupMets, p = 0.7, list = FALSE)

TrainGroup <- GroupMets[ind]
TestGroup <- GroupMets[-ind]

trainMat <- AllMat[, ind]
testMat <- AllMat[, -ind]

###################################################################
## Load TF_MiR genes
TF_MiR <- load("/Users/mohamedomar/Desktop/TF_MiR.rda")

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(trainMat))

## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(trainMat, method = "quantile")
boxplot(usedTrainMat)

usedTestMat <- normalizeBetweenArrays(testMat, method = "quantile")
boxplot(usedTestMat)

#
usedTrainGroup <- TrainGroup
usedTestGroup <- TestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

### Set Feature number and max k
featN <- nrow(usedTrainMat) #nrow(usedTrainMat)
ktsp <- c(1:50)

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
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("Bone_Mets", "Visceral_Mets"), direction = "<",), transpose = TRUE, "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("Bone_Mets", "Visceral_Mets"), direction = "<", ), transpose = TRUE ,"local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = TRUE, print.thres=thr, print.thres.adj=c(0.01,1.25), print.auc=TRUE, print.auc.col="black", levels = c("Bone_Mets", "Visceral_Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Bone_Mets")

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = expr2, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Threshold
thr_test <- coords(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("Bone_Mets", "Visceral_Mets"), direction = "<",),"best", transpose = FALSE)["threshold"]
thr_test

## Print ROC curve local maximas
coords(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("Bone_Mets", "Visceral_Mets"), direction = "<",), transpose = TRUE, "local maximas")

## Plot curve
roc(usedTestGroup, ktspStatsTestRes$statistics, plot = TRUE, print.thres=thr_test, print.auc=TRUE, print.auc.col="black", levels = c("Bone_Mets", "Visceral_Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr_test)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Bone_Mets")

## OR (calculate accuracy, sens, spec and AUC)
SWAP.GetKTSP.PredictionStats(usedTestPredictionRes, usedTestGroup, decision_values = ktspStatsTestRes$statistics)
