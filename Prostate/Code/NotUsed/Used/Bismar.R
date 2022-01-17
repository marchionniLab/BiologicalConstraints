###############################################################################################
## Mohamed Omar
## 21/07/2019
## Goal: Discovery and validation of a small gene signature that can predict the metastasis potential in primary prostate cancer
################################################################################################

## Clean work space
rm(list = ls())

## Set the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")


## Load necessary packages
library(MetaIntegrator)
library(GEOquery)
library(pROC)
library(caret)
library(genefilter)
library(mltools)

########################################################
## Load Training data sets (4)
load("./Data/ProstateData.rda")

## Load Testing data set
## Load Testing data set
load("./Data/Dataset1.rda")
load("./Data/Dataset6.rda")
load("./Objs/Dataset7.rda")
load("./Objs/Dataset8.rda")
load("./Objs/Dataset9.rda")
load("./Objs/Dataset10.rda")


#######################################################

## Getting the phenotype data for each data set
#pheno4 <- pData(Dataset4)
#ProstateData$originalData$GSE46691$pheno <- pheno4
pheno2 <- Dataset2$pheno
pheno3 <- Dataset3$pheno
pheno4 <- Dataset4$pheno
Pheno5 <- Dataset5$pheno
Pheno6 <- Dataset6$pheno
Pheno7 <- Dataset7$pheno
pheno8 <- Dataset8$pheno
pheno9 <- Dataset9$pheno
pheno10 <- Dataset10$pheno

################################

## Getting the expression data for each data set
## load expr4
#load("/Users/mohamedomar/Documents/Research/Projects/Prostate/Data/expr4.rda")
#ProstateData$originalData$GSE46691$expr <- expr4
expr1 <- Dataset1$expr
expr2 <- Dataset2$expr
expr3 <- Dataset3$expr
expr4 <- Dataset4$expr
expr5 <- Dataset5$expr
expr6 <- Dataset6$expr
expr7 <- Dataset7$expr
expr8 <- Dataset8$expr
expr9 <- Dataset9$expr
expr10 <- Dataset10$expr


## Checking if the expression data are normalized and log2 transformed
# boxplot(expr2[,1:15], outline= FALSE)
# boxplot(expr3[,1:15], outline = FALSE)
# boxplot(expr4[,1:15], outline= FALSE)

#################################################################
## Create a list containing training data sets
#AllDataSets <- list(Dataset2, Dataset3, Dataset4, Dataset5)
#names(AllDataSets) <- c("GSE55935", "GSE51066", "GSE46691", "GSE41408")

###################################################################

## Annotate expression

## Expr1
# head(rownames(expr1))
# rownames(expr1) <- Dataset1$keys
# expr1 <- expr1[!is.na(rownames(expr1)), ]
# 
# rownames(expr1) <- gsub("\\,.+", "", rownames(expr1))
# rownames(expr1) <- gsub("\\-.+", "", rownames(expr1))
# 
# expr1 <- aggregate(expr1[,], list(Gene = rownames(expr1)), FUN = mean)
# rownames(expr1) <- expr1$Gene
# expr1$Gene <- NULL
# 
# expr1 <- as.matrix(expr1)
# dim(expr1)
# 
#####################
## expr2
#head(rownames(expr2))
#rownames(expr2) <- Dataset2$keys
expr2 <- expr2[!is.na(rownames(expr2)), ]
dim(expr2)

# X2 <- expr2
# ffun <- filterfun(pOverA(p = 0.8, A = 100))
# filt2 <- genefilter(2^X2, ffun)
# expr2 <- expr2[filt2, ]
#####################
## expr3
#rownames(expr3) <- Dataset3$keys
#expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)
# X3 <- expr3
# filt3 <- genefilter(2^X3, ffun)
# expr3 <- expr3[filt3, ]
# #######################
# dim(expr4)
# X4 <- expr4
# filt4 <- genefilter(2^X4, ffun)
# expr4 <- expr4[filt4, ]

######################
## expr5
head(rownames(expr5))

rownames(expr5) <- Dataset5$keys
expr5 <- expr5[!is.na(rownames(expr5)), ]
dim(expr5)

#####################
## expr6
# processed

#####################
## expr7 
# Processed
expr7 <- expr7+2

#######################
## Expr8 
# Processed

#######################
## Expr9 (Processed)
rownames(expr9) <- gsub("\\,.+", "", rownames(expr9))
rownames(expr9) <- gsub("\\-.+", "", rownames(expr9))

expr9 <- aggregate(expr9[,], list(Gene = rownames(expr9)), FUN = mean)
rownames(expr9) <- expr9$Gene
expr9$Gene <- NULL
expr9 <- as.matrix(expr9)

expr9 <- expr9+4

########################
########################
## Expr10 (Processed)
# rownames(expr10) <- Dataset10$keys
# summary(is.na(rownames(expr10)))
# 
# expr10 <- expr10[!is.na(rownames(expr10)), ]
dim(expr10)
sel = which(apply(expr10, 1, function(x) all(is.finite(x)) ))
expr10 <- expr10[sel,]


# ############################################################
####################################################################
pheno1$Metastasis <- pheno1$`met event (1=yes, 0=no):ch1`
pheno1$Metastasis[pheno1$Metastasis == 0] <- "No_Mets"
pheno1$Metastasis[pheno1$Metastasis == 1] <- "Mets"
pheno1$Metastasis <- as.factor(pheno1$Metastasis)
table(pheno1$Metastasis)
all(rownames(pheno1) == colnames(expr1))

Dataset1$pheno <- pheno1
Dataset1$expr <- expr1
Dataset1$keys <- rownames(expr1)
################

## Modify pheno2
# Remove cell lines
pheno2 <- pheno2[-c(47:54), ]
pheno2$Metastasis <- pheno2$`lymph node metastasis status:ch1`
pheno2$Metastasis[pheno2$Metastasis == 0] <- "No_Mets"
pheno2$Metastasis[pheno2$Metastasis == 1] <- "Mets"
pheno2 <- pheno2[!(pheno2$Metastasis == "NA"), ]

pheno2$Metastasis <- as.factor(pheno2$Metastasis)
table(pheno2$Metastasis)
# Modify expr2
expr2 <- expr2[,colnames(expr2) %in% rownames(pheno2)]
all(rownames(pheno2) == colnames(expr2))
## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset2$expr <- expr2
Dataset2$pheno <- pheno2

####### 
## Modify pheno4
pheno4$Metastasis <- pheno4$`metastatic event:ch1`
pheno4$Metastasis[pheno4$Metastasis == 0] <- "No_Mets"
pheno4$Metastasis[pheno4$Metastasis == 1] <- "Mets"
pheno4 <- pheno4[!(pheno4$Metastasis == "NA"), ]

pheno4$Metastasis <- as.factor(pheno4$Metastasis)
table(pheno4$Metastasis)
levels(pheno4$Metastasis)

## Modify sample names to match sample names of pheno4
head(colnames(expr4))
colnames(expr4) <- gsub(".CEL", "", colnames(expr4))
colnames(expr4) <- gsub(".+\\.", "", colnames(expr4))
expr4 <- expr4[,order(colnames(expr4))]

rownames(pheno4) <- pheno4$description
head(rownames(pheno4))
rownames(pheno4) <- gsub(".+\\.", "", rownames(pheno4))
pheno4 <- pheno4[order(rownames(pheno4)), ]

all(rownames(pheno4) == colnames(expr4))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset4$expr <- expr4
Dataset4$pheno <- pheno4

######
# Modify pheno3
pheno3$Metastasis <- pheno3$`metastatic event:ch1`
pheno3$Metastasis[pheno3$Metastasis == 0] <- "No_Mets"
pheno3$Metastasis[pheno3$Metastasis == 1] <- "Mets"

pheno3$Metastasis <- as.factor(pheno3$Metastasis)
table(pheno3$Metastasis)
levels(pheno3$Metastasis)

all(rownames(pheno3) == colnames(expr3))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset3$pheno <- pheno3
Dataset3$expr <- expr3

#####
## Modify Pheno5
table(Pheno5$`development of metastasis:ch1`)

Pheno5$Metastasis <- Pheno5$`development of metastasis:ch1`

Pheno5$Metastasis[Pheno5$Metastasis == "no"] <- "No_Mets"
Pheno5$Metastasis[Pheno5$Metastasis == "yes"] <- "Mets"

Pheno5$Metastasis <- as.factor(Pheno5$Metastasis)
table(Pheno5$Metastasis)
all(rownames(Pheno5) == colnames(expr5))

Dataset5$expr <- expr5
Dataset5$pheno <- Pheno5

######
# Pheno6
table(Pheno6$`extreme:ch1`)
Pheno6$Metastasis <- Pheno6$`extreme:ch1`

Pheno6$Metastasis[Pheno6$Metastasis == "Indolent"] <- "No_Mets"
Pheno6$Metastasis[Pheno6$Metastasis == "Lethal"] <- "Mets"

Pheno6$Metastasis <- as.factor(Pheno6$Metastasis)
table(Pheno6$Metastasis)
all(rownames(Pheno6) == colnames(expr6))

Dataset6$expr <- expr6
Dataset6$pheno <- Pheno6

########
# Pheno7
# Processed
Pheno7$Metastasis <- as.factor(Pheno7$Metastasis) 

Dataset7$expr <- expr7
Dataset7$pheno <- Pheno7
Dataset7$keys <- rownames(expr7)

#########
## Pheno8 
# processed

###########
## Pheno9 
# Processed
#pheno9$Metastasis <- as.factor(pheno9$`recurrence status:ch1`)
#table(pheno9$Metastasis)
#levels(pheno9$Metastasis) <- c("No_Mets", "Mets")

#all(rownames(pheno9) == colnames(expr9))

#Dataset9$expr <- expr9
#Dataset9$pheno <- pheno9
#Dataset9$keys <- rownames(expr9)

Dataset9$expr <- expr9
Dataset9$pheno <- pheno9
Dataset9$keys <- rownames(expr9)

###########
## Pheno10
# Processed
pheno10$Metastasis <- pheno10$`pathology stage:ch1`
pheno10$Metastasis <- gsub("T.+M", "M", pheno10$Metastasis)
pheno10$Metastasis <- gsub("N.+", "", pheno10$Metastasis)

table(pheno10$Metastasis)
pheno10 <- pheno10[!(pheno10$Metastasis == "Mx"), ]
pheno10 <- pheno10[!(pheno10$Metastasis == "pMx"), ]
pheno10 <- pheno10[!(pheno10$Metastasis == "U"), ]
pheno10$Metastasis <- as.factor(pheno10$Metastasis)
levels(pheno10$Metastasis) <- c("No_Mets", "No_Mets", "Mets", "Mets")
table(pheno10$Metastasis)

expr10 <- expr10[, colnames(expr10) %in% rownames(pheno10)]
all(rownames(pheno10) == colnames(expr10))

Dataset10$expr <- expr10
Dataset10$pheno <- pheno10
Dataset10$keys <- rownames(expr10)

#########################################################################
############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘disease’ or 0 for ‘control’)
Dataset1 <- classFunction(Dataset1, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset2 <- classFunction(Dataset2, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset3 <- classFunction(Dataset3, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset4 <- classFunction(Dataset4, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset5 <- classFunction(Dataset5, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset9 <- classFunction(Dataset9, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset10 <- classFunction(Dataset10, column = "Metastasis", diseaseTerms = c("Mets"))


############################################################################
############################################################################
#########################################################################
## The metaanalysis

## Creating the meta object
AllDataSets <- list(Dataset2, Dataset3, Dataset4, Dataset5, Dataset9, Dataset10)
names(AllDataSets) <- c(Dataset2$formattedName, Dataset3$formattedName, Dataset4$formattedName, Dataset5$formattedName, Dataset9$formattedName, Dataset10$formattedName)

Prostate_meta <- list()
Prostate_meta$originalData <- AllDataSets

## Replace keys within each data set
#Prostate_meta$originalData$GSE116918$keys <- rownames(expr1)
Prostate_meta$originalData$GSE55935$keys <- rownames(expr2)
Prostate_meta$originalData$GSE51066$keys <- rownames(expr3)
Prostate_meta$originalData$GSE46691$keys <- rownames(expr4)
Prostate_meta$originalData$GSE41408$keys <- rownames(expr5)
#Prostate_meta$originalData$GSE16560$keys <- rownames(expr6)
Prostate_meta$originalData$GSE25136$keys <- rownames(expr9)
Prostate_meta$originalData$GSE70769$keys <- rownames(expr10)


## Check the meta object before the metaanalysis
checkDataObject(Prostate_meta, "Meta", "Pre-Analysis") ## If true, Proceed to the meta analysis


## Run the meta analysis
Prostate_metaanalysis <- runMetaAnalysis(Prostate_meta, runLeaveOneOutAnalysis = F, maxCores = 3)

## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
Prostate_metaanalysis <- filterGenes(Prostate_metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.1, numberStudiesThresh = 6)


## Stats of all genes
AllGenes <- Prostate_metaanalysis$metaAnalysis$pooledResults

## The OncoDx genes (15 genes)
Bismar_Sig <- c("FLNA", "AMACR", "CDK7", "ITGA5", "JAG1", "SLC4A1AP", "MIB1", "MTA1", "MUC1", "TP63", "KLK3", "TPD52")

#OncoDx_Genes <- c("BGN", "COL1A1", "SFRP4", "FLNC", "GSN", "TPM2", "GSTM2", "TPX2", "LAMB3", "FAM13C", "KLK2", "AZGP1", "SRD5A2", "DUSP1", "FOS")
save(Bismar_Sig, file = "./Objs/Bismar_Sig_Genes.rda")



summary(Bismar_Sig %in% rownames(AllGenes))


# Keep only the stats of OncoDx genes (2 genes are absent)
Bismar_Sig <- AllGenes[Bismar_Sig, ]
#OncoDx <- na.omit(OncoDx)

# Replace the gene names with the OncoDx genes
Prostate_metaanalysis$filterResults$FDR0.1_es0_nStudies6_looaFALSE_hetero0$posGeneNames <- c("JAG1", "SLC4A1AP", "MIB1", "MUC1", "TPD52")
Prostate_metaanalysis$filterResults$FDR0.1_es0_nStudies6_looaFALSE_hetero0$negGeneNames <- c("FLNA", "AMACR", "CDK7", "ITGA5", "MTA1", "TP63", "KLK3")

Bismar_Filter <- Prostate_metaanalysis$filterResults[[1]]
Bismar_Filter_Summary <- summarizeFilterResults(metaObject = Prostate_metaanalysis, getMostRecentFilter(Prostate_metaanalysis))


## Create a summary ROC curve (Training data sets)
set.seed(333)
pooledROCPlot(metaObject = Prostate_metaanalysis, filterObject = Bismar_Filter)
###########################################################################################
#############################################################################

## The next step is validation on an indepndent data set (GSE116918)

# Get the validation data set
#valDataSet <- getGEOData("GSE116918")
valDataSet1 <- Dataset1
valDataSet2 <- Dataset6
valDataSet3 <- Dataset7
valDataSet4 <- Dataset8
#valDataSet5 <- Dataset9
#valDataSet6 <- Dataset10

# Acess the phenotype and expression data
val_pheno1 <- Dataset1$pheno
val_expr1 <- Dataset1$expr

val_pheno2 <- Pheno6
val_expr2 <- expr6

val_pheno3 <- Pheno7
val_expr3 <- expr7

val_pheno4 <- pheno8
val_expr4 <- expr8

#val_pheno5 <- pheno9
#val_expr5 <- expr9

#val_pheno6 <- pheno10
#val_expr6 <- expr10


## Modify val_expr1
rownames(val_expr1) <- valDataSet1$keys
summary(is.na(rownames(val_expr1)))
val_expr1 <- val_expr1[!is.na(rownames(val_expr1)), ]
dim(val_expr1)


## Modify val_pheno1
val_pheno1$Metastasis <- val_pheno1$`met event (1=yes, 0=no):ch1`
val_pheno1$Metastasis[val_pheno1$Metastasis == 0] <- "No_Mets"
val_pheno1$Metastasis[val_pheno1$Metastasis == 1] <- "Mets"
val_pheno1$Metastasis <- as.factor(val_pheno1$Metastasis)
table(val_pheno1$Metastasis)
all(rownames(val_pheno1) == colnames(val_expr1))

valDataSet1$pheno <- val_pheno1
valDataSet1$expr <- val_expr1
valDataSet1$keys <- rownames(val_expr1)

## Label the samples
valDataSet1 <- classFunction(valDataSet1, column = "Metastasis", diseaseTerms = c("Mets"))
valDataSet2 <- classFunction(valDataSet2, column = "Metastasis", diseaseTerms = c("Mets"))
valDataSet3 <- classFunction(valDataSet3, column = "Metastasis", diseaseTerms = c("Mets"))
valDataSet4 <- classFunction(valDataSet4, column = "Metastasis", diseaseTerms = c("Mets"))
#valDataSet5 <- classFunction(valDataSet5, column = "Metastasis", diseaseTerms = c("Mets"))
#valDataSet6 <- classFunction(valDataSet6, column = "Metastasis", diseaseTerms = c("Mets"))

##############################################
## Examine the performance of OncoDx filter

## ROC plot

# Validation dataset1
rocPlot(datasetObject = valDataSet1, filterObject = Bismar_Filter)

# Validation dataset2
valDataSet2$keys <- rownames(val_expr2)
rocPlot(datasetObject = valDataSet2, filterObject = Bismar_Filter)

# Validation dataset3
valDataSet3$keys <- rownames(val_expr3)
rocPlot(datasetObject = valDataSet3, filterObject = Bismar_Filter)

# Validation dataset4
valDataSet4$keys <- rownames(val_expr4)
rocPlot(datasetObject = valDataSet4, filterObject = Bismar_Filter)

# Validation dataset5
#valDataSet5$keys <- rownames(val_expr5)
#rocPlot(datasetObject = valDataSet5, filterObject = OncoDx_Filter)

# Validation dataset6
#valDataSet6$keys <- rownames(val_expr6)
#rocPlot(datasetObject = valDataSet6, filterObject = OncoDx_Filter)

###########################
# PRC plot
prcPlot(datasetObject = valDataSet1, filterObject = Bismar_Filter)

# val dataset2
prcPlot(datasetObject = valDataSet2, filterObject = Bismar_Filter)

# val dataset3
prcPlot(datasetObject = valDataSet3, filterObject = Bismar_Filter)

# val dataset4
prcPlot(datasetObject = valDataSet4, filterObject = Bismar_Filter)

# val dataset5
#prcPlot(datasetObject = valDataSet5, filterObject =OncoDx_Filter)

# val dataset6
#prcPlot(datasetObject = valDataSet6, filterObject =OncoDx_Filter)

##########################
# Violin plot
violinPlot(filterObject = Bismar_Filter, datasetObject = valDataSet1, labelColumn = "Metastasis")

violinPlot(filterObject = Bismar_Filter, datasetObject = valDataSet2, labelColumn = "Metastasis")

violinPlot(filterObject = Bismar_Filter, datasetObject = valDataSet3, labelColumn = "Metastasis")

violinPlot(filterObject = Bismar_Filter, datasetObject = valDataSet4, labelColumn = "Metastasis")

#violinPlot(filterObject = OncoDx_Filter, datasetObject = valDataSet5, labelColumn = "Metastasis")

#violinPlot(filterObject = OncoDx_Filter, datasetObject = valDataSet6, labelColumn = "Metastasis")

##############
## Calculate a signature score (Z score) and add it to the phenotype table
val_pheno1$score <- calculateScore(filterObject = Bismar_Filter, datasetObject = valDataSet1)

# Val Dataset2
val_pheno2$score <- calculateScore(filterObject = Bismar_Filter, datasetObject = valDataSet2)

# Val Dataset3
val_pheno3$score <- calculateScore(filterObject = Bismar_Filter, datasetObject = valDataSet3)

# Val Dataset4
val_pheno4$score <- calculateScore(filterObject = Bismar_Filter, datasetObject = valDataSet4)

# Val Dataset5
#val_pheno5$score <- calculateScore(filterObject = OncoDx_Filter, datasetObject = valDataSet5)

# Val Dataset6
#val_pheno6$score <- calculateScore(filterObject = OncoDx_Filter, datasetObject = valDataSet6)


## Find the best threshold (for further use for the classification in the Test data)
thr_test1 <- coords(roc(val_pheno1$Metastasis, val_pheno1$score, levels = c("No_Mets", "Mets")),"best", transpose = TRUE)["threshold"]
thr_test1

# Val Dataset2
thr_test2 <- coords(roc(val_pheno2$Metastasis, val_pheno2$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
thr_test2 

# Val Dataset3
thr_test3 <- coords(roc(val_pheno3$Metastasis, val_pheno3$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
thr_test3 

# Val Dataset4
thr_test4 <- coords(roc(val_pheno4$Metastasis, val_pheno4$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
thr_test4 

# Val Dataset5
#thr_test5 <- coords(roc(val_pheno5$Metastasis, val_pheno5$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
#thr_test5 

# Val Dataset6
#thr_test6 <- coords(roc(val_pheno6$Metastasis, val_pheno6$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
#thr_test6 

## Find the optimal trade off between the sensitivity and specificity (We want good sensitivity with preservation of a decent accuracy)
coords(roc(val_pheno1$Metastasis, val_pheno1$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")
# Val Dataset2
coords(roc(val_pheno2$Metastasis, val_pheno2$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")
# Val Dataset3
coords(roc(val_pheno3$Metastasis, val_pheno3$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")
# Val Dataset4
coords(roc(val_pheno4$Metastasis, val_pheno4$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")
# Val Dataset5
#coords(roc(val_pheno5$Metastasis, val_pheno5$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")
# Val Dataset6
#coords(roc(val_pheno6$Metastasis, val_pheno6$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")


## Predictions
# Val Dataset1
Test_Predictions1 <- ifelse(val_pheno1$score >= thr_test1, "Mets", "No_Mets")
confusionMatrix(as.factor(Test_Predictions1), val_pheno1$Metastasis, positive = "Mets")

mcc(preds = as.factor(Test_Predictions1), actuals = val_pheno1$Metastasis)

# Val Dataset2
Test_Predictions2 <- ifelse(val_pheno2$score >= thr_test2, "Mets", "No_Mets")
confusionMatrix(as.factor(Test_Predictions2), val_pheno2$Metastasis, positive = "Mets")

mcc(preds = as.factor(Test_Predictions2), actuals = val_pheno2$Metastasis)

# Val Dataset3
Test_Predictions3 <- ifelse(val_pheno3$score >= thr_test3, "Mets", "No_Mets")
confusionMatrix(as.factor(Test_Predictions3), val_pheno3$Metastasis, positive = "Mets")

mcc(preds = as.factor(Test_Predictions3), actuals = val_pheno3$Metastasis)

# Val Dataset4
Test_Predictions4 <- ifelse(val_pheno4$score >= thr_test4, "Mets", "No_Mets")
confusionMatrix(as.factor(Test_Predictions4), val_pheno4$Metastasis, positive = "Mets")

mcc(preds = as.factor(Test_Predictions4), actuals = val_pheno4$Metastasis)

# Val Dataset5
#Test_Predictions5 <- ifelse(val_pheno5$score >= thr_test5, "Mets", "No_Mets")
#confusionMatrix(as.factor(Test_Predictions5), val_pheno5$Metastasis, positive = "Mets")

#mcc(preds = as.factor(Test_Predictions5), actuals = val_pheno5$Metastasis)

# Val Dataset6
#Test_Predictions6 <- ifelse(val_pheno6$score >= thr_test6, "Mets", "No_Mets")
#confusionMatrix(as.factor(Test_Predictions6), val_pheno6$Metastasis, positive = "Mets")

#mcc(preds = as.factor(Test_Predictions6), actuals = val_pheno6$Metastasis)

#########################
## Meta -ROC Test
## Meta-Test
Meta_Test <- list() 
Meta_Test$originalData$GSE116918 <- valDataSet1
#Meta_Test$originalData$GSE16560 <- valDataSet2
Meta_Test$originalData$JHU <- valDataSet3
Meta_Test$originalData$HPFS <- valDataSet4
#Meta_Test$originalData$GSE25136 <- valDataSet5
#Meta_Test$originalData$GSE70769 <- valDataSet6

#Meta_Test$originalData$GSE16560$keys <- rownames(Meta_Test$originalData$GSE16560$expr)
Meta_Test$originalData$JHU$keys <- rownames(Meta_Test$originalData$JHU$expr)
Meta_Test$originalData$HPFS$keys <- rownames(Meta_Test$originalData$HPFS$expr)
#Meta_Test$originalData$GSE25136$keys <- rownames(Meta_Test$originalData$GSE25136$expr)
#Meta_Test$originalData$GSE70769$keys <- rownames(Meta_Test$originalData$GSE70769$expr)


Meta_Test$originalData$GSE116918 <- classFunction(Meta_Test$originalData$GSE116918, column = "Metastasis", diseaseTerms = c("Mets"))
#Meta_Test$originalData$GSE16560 <- classFunction(Meta_Test$originalData$GSE16560, column = "Metastasis", diseaseTerms = c("Mets"))
Meta_Test$originalData$JHU <- classFunction(Meta_Test$originalData$JHU, column = "Metastasis", diseaseTerms = c("Mets"))
Meta_Test$originalData$HPFS <- classFunction(Meta_Test$originalData$HPFS, column = "Metastasis", diseaseTerms = c("Mets"))
#Meta_Test$originalData$GSE25136 <- classFunction(Meta_Test$originalData$GSE25136, column = "Metastasis", diseaseTerms = c("Mets"))
#Meta_Test$originalData$GSE70769 <- classFunction(Meta_Test$originalData$GSE70769, column = "Metastasis", diseaseTerms = c("Mets"))


set.seed(333)
summaryROCPlot(metaObject = Meta_Test, filterObject = Bismar_Filter)

set.seed(333)
png(filename = "./Objs/Meta/Bismar_pooledROC_Testing.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = Meta_Test, filterObject = Bismar_Filter)
dev.off()

############################################


# Survival (Meta-score)

## Load necessary packages
library(survival)
library(survminer)

## Turn the score column into numeric
val_pheno1[,"score"] <- as.numeric(val_pheno1[,"score"]) 
val_pheno2[,"score"] <- as.numeric(val_pheno2[, "score"])
val_pheno3[,"score"] <- as.numeric(val_pheno3[, "score"])

## Divide the dataset into quantiles based on the Meta score risk.
#quantiles <- quantile(val_pheno1$score, probs=c(0.33333,0.66667))
# val_Pheno1
val_pheno1$Meta_Score <- val_pheno1$score
val_pheno1[which(val_pheno1$score <= thr_test1), "Meta_Score"] = "low"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
val_pheno1[which(val_pheno1$score > thr_test1),"Meta_Score"] = "high"

# val_pheno2
val_pheno2$Meta_Score <- val_pheno2$score
val_pheno2[which(val_pheno2$score <= thr_test2), "Meta_Score"] = "low"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
val_pheno2[which(val_pheno2$score > thr_test2),"Meta_Score"] = "high"

# val_pheno3
val_pheno3$Meta_Score <- val_pheno3$score
val_pheno3[which(val_pheno3$score <= thr_test3), "Meta_Score"] = "low"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
val_pheno3[which(val_pheno3$score > thr_test3),"Meta_Score"] = "high"


val_pheno1$Time <- as.numeric(val_pheno1$`follow-up time (met, months):ch1`)
val_pheno2$Time <- as.numeric(val_pheno2$`fup.month:ch1`)
val_pheno3$Time <- as.numeric(val_pheno3$met_time)

val_pheno1$Event <- as.numeric(val_pheno1$`met event (1=yes, 0=no):ch1`)

val_pheno2$Event <- as.numeric(val_pheno2$Metastasis)
val_pheno2$Event[val_pheno2$Event == 2] <- 0 

val_pheno3$Event <- as.numeric(val_pheno3$Metastasis)
val_pheno3$Event[val_pheno3$Event == 2] <- 0 
table(val_pheno3$Event)

## Keep only relevant information from the phenotype table
surv_data1 <- val_pheno1[,c("Time","Event","score","Meta_Score","gleason grade:ch1", "patient age (years):ch1", "psa (ng/ml):ch1", "t-stage:ch1")]
surv_data1$Meta_Score <- factor(surv_data1$Meta_Score, levels = c("low", "high"))

surv_data2 <- val_pheno2[, c("Time", "Event", "score", "Meta_Score", "gleason:ch1", "age:ch1", "cancer.percent:ch1")]
surv_data2$Meta_Score <- factor(surv_data2$Meta_Score, levels = c("low", "high"))

surv_data3 <- val_pheno3[, c("Time", "Event", "score", "Meta_Score", "Re.Gleasonsum", "age.y", "ERG", "re.PositiveMargin", "re.Stage", "PTEN_LOSS")]
surv_data3$Meta_Score <- factor(surv_data3$Meta_Score, levels = c("low", "high"))

# Create a survival object
surv_data.surv1 <-  with(surv_data1, Surv(Time, Event == 1))
surv_data.surv2 <-  with(surv_data2, Surv(Time, Event == 1))
surv_data.surv3 <-  with(surv_data3, Surv(Time, Event == 1))

#Calculate p-value
survdifftest1 <- survdiff(surv_data.surv1 ~ Meta_Score, data = surv_data1)
survpvalue1 <- 1 - pchisq(survdifftest1$chisq, length(survdifftest1$n) - 1)

survdifftest2 <- survdiff(surv_data.surv2 ~ Meta_Score, data = surv_data2)
survpvalue2 <- 1 - pchisq(survdifftest2$chisq, length(survdifftest2$n) - 1)

survdifftest3 <- survdiff(surv_data.surv3 ~ Meta_Score, data = surv_data3)
survpvalue3 <- 1 - pchisq(survdifftest3$chisq, length(survdifftest3$n) - 1)

## Create a linear test p-value
# surv_data_lin <- val_pheno1[,c("Time","Event","Meta_Score")]
# surv_data_lin[,"Meta_Score"] = as.vector(surv_data_lin[,"Meta_Score"])
# surv_data_lin[which(surv_data_lin[,"Meta_Score"]=="low"),"Meta_Score"] = 1
# surv_data_lin[which(surv_data_lin[,"Meta_Score"]=="int"),"Meta_Score"] = 2
# surv_data_lin[which(surv_data_lin[,"Meta_Score"]=="high"),"Meta_Score"] = 3
# surv_data_lin[ , "Meta_Score"] <- as.numeric(surv_data_lin[ , "Meta_Score"])
# survpvalue_linear <- summary(coxph(Surv(Time, as.numeric(Event))~Meta_Score, data=surv_data_lin))$sctest[3]
# survpvalue_linear <-format(as.numeric(survpvalue_linear), digits=3)


## Kaplan-Meier curve
png(filename = "./Figs/Meta/Bismar_KM_Survival_GSE116918.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data1),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the Bismar signature Meta score in GSE116918", xlab = "Time in months", ylab = "Survival")

dev.off()

###############

png(filename = "./Figs/Meta/KM_Survival2.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data2),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the gene signature Meta score in the Test data", xlab = "Time in months", ylab = "Survival")

dev.off()

###################
png(filename = "./Figs/Meta/Bismar_KM_Survival_JHU.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data3),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the Bismar signature Meta score in the JHU cohort", xlab = "Time in months", ylab = "Survival")

dev.off()


#############################################################################
###### COX 

## First validation dataset
Cox_data1 <- surv_data1
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data1$Age <- ifelse(Cox_data1$`patient age (years):ch1` >= 70, ">=70 years", "<70 years")
Cox_data1$Age <- factor(Cox_data1$Age, levels = c("<70 years", ">=70 years"))

Cox_data1$PSA <- ifelse(Cox_data1$`psa (ng/ml):ch1` > 20, "High PSA (> 20)", "Low PSA (<= 20)")
Cox_data1$PSA <- factor(Cox_data1$PSA, levels = c("Low PSA (<= 20)", "High PSA (> 20)"))

Cox_data1$Gleason <- ifelse(Cox_data1$`gleason grade:ch1` == 6 , "Gleason = 6", "Gleason > 6")
Cox_data1$Gleason <- factor(Cox_data1$Gleason, levels = c("Gleason = 6", "Gleason > 6"))

Cox_data1$T_Stage <- ifelse(Cox_data1$`t-stage:ch1` == c("T1", "T2"), "T1/T2", "T3/T4/Unknown")
Cox_data1$T_Stage <- factor(Cox_data1$T_Stage, levels = c("T1/T2", "T3/T4/Unknown"))

## Univariate Cox model
# PFS
covariates1 <- c("Meta_Score", "Age", "PSA", "Gleason", "T_Stage")
univ_formulas1 <- sapply(covariates1,
                         function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models1 <- lapply( univ_formulas1, function(x){coxph(x, data = Cox_data1)})
# Extract data 
univ_results1 <- lapply(univ_models1,
                        function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                        })

UniCOX1 <- t(as.data.frame(univ_results1, check.names = F))
as.data.frame(UniCOX1)
########################################

## Multivariate

CoxModel1 <- coxph(Surv(Time, Event)~ Meta_Score+Age+PSA+Gleason+T_Stage, data=Cox_data1)

png(filename = "./Figs/Meta/Bismar_CoxModel_GSE116918.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel1, data = Cox_data1, main = "Cox proportional hazards model using Bismar signature in GSE116918")
dev.off()

############################################################################
############################################################################
## Second validation dataset
Cox_data2 <- surv_data2
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data2$Age <- ifelse(Cox_data2$`age:ch1` >= 70, ">=70 years", "<70 years")
Cox_data2$Age <- factor(Cox_data2$Age, levels = c("<70 years", ">=70 years"))

#Cox_data1$PSA <- ifelse(Cox_data1$`psa (ng/ml):ch1` > 20, "High PSA (> 20)", "Low PSA (<= 20)")
#Cox_data1$PSA <- factor(Cox_data1$PSA, levels = c("Low PSA (<= 20)", "High PSA (> 20)"))

Cox_data2$Gleason <- ifelse(Cox_data2$`gleason:ch1` == 6 , "Gleason = 6", "Gleason > 6")
Cox_data2$Gleason <- factor(Cox_data2$Gleason, levels = c("Gleason = 6", "Gleason > 6"))

Cox_data2$`cancer.percent:ch1` <- as.integer(Cox_data2$`cancer.percent:ch1`)
Cox_data2$CancerPercent <- ifelse(Cox_data2$`cancer.percent:ch1` >= 50, ">= 50", "< 50")
Cox_data2$CancerPercent <- factor(Cox_data2$CancerPercent, levels = c("< 50", ">= 50"))

## Univariate Cox model
# PFS
covariates2 <- c("Meta_Score", "Age", "Gleason", "CancerPercent")
univ_formulas2 <- sapply(covariates2,
                         function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models2 <- lapply( univ_formulas2, function(x){coxph(x, data = Cox_data2)})
# Extract data 
univ_results2 <- lapply(univ_models2,
                        function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                        })

UniCOX2 <- t(as.data.frame(univ_results2, check.names = F))
as.data.frame(UniCOX2)
########################################

## Multivariate

CoxModel2 <- coxph(Surv(Time, Event)~ Meta_Score+Age+Gleason+CancerPercent, data=Cox_data2)

png(filename = "./Figs/Meta/CoxModel2.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel2, data = Cox_data2, main = "Cox proportional hazards model")
dev.off()

#######################################################
#######################################################
## Third validation dataset
Cox_data3 <- surv_data3
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data3$Age <- ifelse(Cox_data3$age.y >= 65, ">=65 years", "<65 years")
Cox_data3$Age <- factor(Cox_data3$Age, levels = c("<65 years", ">=65 years"))

#Cox_data1$PSA <- ifelse(Cox_data1$`psa (ng/ml):ch1` > 20, "High PSA (> 20)", "Low PSA (<= 20)")
#Cox_data1$PSA <- factor(Cox_data1$PSA, levels = c("Low PSA (<= 20)", "High PSA (> 20)"))

Cox_data3$Gleason <- ifelse(Cox_data3$Re.Gleasonsum <= 7 , "Gleason <= 7", "Gleason > 7")
Cox_data3$Gleason <- factor(Cox_data3$Gleason, levels = c("Gleason <= 7", "Gleason > 7"))

Cox_data3$ERG <- as.factor(Cox_data3$ERG)
table(Cox_data3$ERG)
levels(Cox_data3$ERG)

Cox_data3$Margin <- as.factor(Cox_data3$re.PositiveMargin)
table(Cox_data3$Margin)
levels(Cox_data3$Margin) <- c("Negative", "Negative", "Positive", "Positive")


Cox_data3$Stage <- as.factor(Cox_data3$re.Stage)
table(Cox_data3$Stage)
levels(Cox_data3$Stage) <- c("< T3", "< T3", "< T3", "< T3", "< T3", "< T3", "T3", "T3")

Cox_data3$PTEN <- as.factor(Cox_data3$PTEN_LOSS)
table(Cox_data3$PTEN)
levels(Cox_data3$PTEN)

## Univariate Cox model
# PFS
covariates3 <- c("Meta_Score", "Age", "Gleason", "ERG", "Margin", "Stage", "PTEN")
univ_formulas3 <- sapply(covariates3,
                         function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models3 <- lapply( univ_formulas3, function(x){coxph(x, data = Cox_data3)})
# Extract data 
univ_results3 <- lapply(univ_models3,
                        function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                        })

UniCOX3 <- t(as.data.frame(univ_results3, check.names = F))
as.data.frame(UniCOX3)
########################################

## Multivariate

CoxModel3 <- coxph(Surv(Time, Event)~ Meta_Score+Age+Gleason+ERG+Margin+Stage+PTEN, data=Cox_data3)

png(filename = "./Figs/Meta/Bismar_CoxModel_JHU.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel3, data = Cox_data3, main = "Cox proportional hazards model using Bismar signature in JHU cohort")
dev.off()

