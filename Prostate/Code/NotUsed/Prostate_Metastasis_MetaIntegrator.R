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

## Load data

#load("./Data/Datasets.rda")

#ProstateData <- getGEOData(c("GSE55935", "GSE46691", "GSE41408"))
#save(ProstateData, file = "./Data/ProstateData.rda")

## Load Training data sets (4)
load("./Data/ProstateData.rda")

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
pheno1 <- Dataset1$pheno
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
# expr1 <- as.matrix(expr1)
# dim(expr1)
# X1 <- expr1
# ffun <- filterfun(pOverA(p = 0.5, A = 50))
# filt1 <- genefilter(2^X1, ffun)
# expr1 <- expr1[filt1, ]
# dim(expr1)

#expr1 <- t(scale(t(expr1), center = T, scale = T))
#####################
## expr2
#head(rownames(expr2))
#rownames(expr2) <- Dataset2$keys
expr2 <- expr2[!is.na(rownames(expr2)), ]
dim(expr2)

X2 <- expr2
ffun <- filterfun(cv(a = 0.1, b = 10))
filt2 <- genefilter(2^X2, ffun)
table(filt2)
expr2 <- expr2[filt2, ]
dim(expr2)

#expr2 <- t(scale(t(expr2), center = T, scale = T))

#####################
## expr3
#rownames(expr3) <- Dataset3$keys
#expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)

X3 <- expr3
#ffun <- filterfun(cv(a = 0.1, b = 1))
filt3 <- genefilter(2^X3, ffun)
table(filt3)
expr3 <- expr3[filt3, ]

#expr3 <- t(scale(t(expr3), center = T, scale = T))

# #######################
dim(expr4)
X4 <- expr4
filt4 <- genefilter(2^X4, ffun)
table(filt4)
expr4 <- expr4[filt4, ]
# 
#expr4 <- t(scale(t(expr4), center = T, scale = T))

######################
## expr5
head(rownames(expr5))

rownames(expr5) <- Dataset5$keys
expr5 <- expr5[!is.na(rownames(expr5)), ]
dim(expr5)
# 
X5 <- expr5
filt5 <- genefilter(2^X5, ffun)
table(filt5)
expr5 <- expr5[filt5, ]

#expr5 <- t(scale(t(expr5), center = T, scale = T))

#####################
## expr6
# processed

dim(expr6)
#X6 <- expr6
#filt6 <- genefilter(2^X6, ffun)
#expr6 <- expr6[filt6, ]

#expr6 <- t(scale(t(expr6), center = T, scale = T))

#####################
## Expr7
# Processed
#range(expr7)
#expr7 <- expr7+2

# dim(expr7)
# X7 <- expr7
# filt7 <- genefilter(2^X7, ffun)
# expr7 <- expr7[filt7, ]

expr7 <- expr7+2

#expr7 <- t(scale(t(expr7), center = T, scale = T))

#######################
## Expr8 
# Processed

#expr8 <- t(scale(t(expr8), center = T, scale = T))

#######################
## Expr9 (Processed)
#rownames(expr9) <- Dataset9$keys
#summary(is.na(rownames(expr9)))

#expr9 <- expr9[!is.na(rownames(expr9)), ]
#dim(expr9)

rownames(expr9) <- gsub("\\,.+", "", rownames(expr9))
rownames(expr9) <- gsub("\\-.+", "", rownames(expr9))

expr9 <- aggregate(expr9[,], list(Gene = rownames(expr9)), FUN = mean)
rownames(expr9) <- expr9$Gene
expr9$Gene <- NULL
expr9 <- as.matrix(expr9)

#expr9 <- expr9+4

dim(expr9)
X9 <- expr9
#ffun <- filterfun(cv(a = 0.1, b = 10))
filt9 <- genefilter(2^X9, ffun)
table(filt9)
expr9 <- expr9[filt9, ]

#expr9 <- t(scale(t(expr9), center = T, scale = T))

########################
## Expr10 (Processed)
# rownames(expr10) <- Dataset10$keys
# summary(is.na(rownames(expr10)))
# 
# expr10 <- expr10[!is.na(rownames(expr10)), ]
dim(expr10)
sel = which(apply(expr10, 1, function(x) all(is.finite(x)) ))
expr10 <- expr10[sel,]

dim(expr10)
X10 <- expr10
filt10 <- genefilter(2^X10, ffun)
table(filt10)
expr10 <- expr10[filt10, ]

#expr10 <- t(scale(t(expr10), center = T, scale = T))

#############################################
### Common genes
# AllExpr <- list(expr1, expr2, expr3, expr4, expr5, expr6)
# commonGenes <- Reduce("intersect", lapply(AllExpr, rownames))
# 
# expr1 <- expr1[commonGenes, ]
# expr2 <- expr2[commonGenes, ]
# expr3 <- expr3[commonGenes, ]
# expr4 <- expr4[commonGenes, ]
# expr5 <- expr5[commonGenes, ]
# expr6 <- expr6[commonGenes, ]
# 
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

Dataset8$expr <- expr8
Dataset8$pheno <- pheno8
Dataset8$keys <- rownames(expr8)
###########
## Pheno9 
# Processed
#pheno9$Metastasis <- as.factor(pheno9$`recurrence status:ch1`)
#table(pheno9$Metastasis)
#levels(pheno9$Metastasis) <- c("No_Mets", "Mets")

#all(rownames(pheno9) == colnames(expr9))

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
Dataset6 <- classFunction(Dataset6, column = "Metastasis", diseaseTerms = c("Mets"))
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

Prostate_meta <- geneSymbolCorrection(Prostate_meta)

## Run the meta analysis
Prostate_metaanalysis <- runMetaAnalysis(Prostate_meta, runLeaveOneOutAnalysis = F, maxCores = 3)


## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
Prostate_metaanalysis <- filterGenes(Prostate_metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.15, numberStudiesThresh = 6)

## Assigning a name to the filter
filter <- Prostate_metaanalysis$filterResults[[1]]
filter

PositiveGenes <- filter$posGeneNames
NegativeGenes <- filter$negGeneNames
save(PositiveGns, NegativeGns, file = "./Objs/NewSigGenes_Pre.rda")

## Summarize filter results
filter_summary <- summarizeFilterResults(metaObject = Prostate_metaanalysis, getMostRecentFilter(Prostate_metaanalysis))

## Save the filter
save(filter, file = "./Objs/Meta/filter.rda")

## Save a table of the positive genes and negative genes
write.table(filter_summary$pos, file = "./Objs/Meta/Positive_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")
write.table(filter_summary$neg, file = "./Objs/Meta/Negative_genes_filter.csv", quote = TRUE, sep = "\t", col.names = TRUE, row.names = TRUE, dec = ".")

## Modify the gene signature for more accuracy and AUC
# Using forward search 
New_filter <- forwardSearch(metaObject = Prostate_metaanalysis, filterObject = filter)
  
save(New_filter, file = "./Objs/Meta/NewFilter.rda")

## Replace the old filter with the new smaller one
Prostate_metaanalysis$filterResults$FDR0.15_es0_nStudies6_looaFALSE_hetero0$posGeneNames <- New_filter$posGeneNames
Prostate_metaanalysis$filterResults$FDR0.15_es0_nStudies6_looaFALSE_hetero0$negGeneNames <- New_filter$negGeneNames

New_filter <- Prostate_metaanalysis$filterResults[[1]]
New_filter_summary <- summarizeFilterResults(metaObject = Prostate_metaanalysis, getMostRecentFilter(Prostate_metaanalysis))

## Save the tables of positive and negative genes
write.table(New_filter_summary$pos, file = "./Objs/NewFilter_Positive_genes.csv", quote = T, sep = "\t", col.names = T, row.names = F)
write.table(New_filter_summary$neg, file = "./Objs/NewFilter_Negative_genes.csv", quote = T, sep = "\t", col.names = T, row.names = F)

## Gene names
PositiveGenes <- New_filter$posGeneNames
PositiveGenes
NegativeGenes <- New_filter$negGeneNames
NegativeGenes

New_SignatureGenes <- c(PositiveGenes, NegativeGenes)
save(New_SignatureGenes, file = "./Objs/NewSigGenes.rda")

# Compare with OncoDx
load("./Objs/OncoDxGenes.rda")
intersect(New_SignatureGenes, OncoDx_Genes)


## Create a summary ROC curve (Training data sets)
set.seed(333)
png(filename = "./Figs/Meta/Pooled_ROC_TrainingDatasets.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = Prostate_metaanalysis, filterObject = New_filter)
dev.off()

## Effect size visualization in each study
# For the up-regulated genes
pdf(file = "./Figs/Meta/UpGenes_effect_size.pdf", title = "Effect size of up-regulated genes across the discovery data sets", width = 20, height = 10)
par(mfrow=c(2,2))

forestPlot(metaObject = Prostate_metaanalysis, geneName = "CAMK2N1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "ASNS", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "NCOA2", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "PTPN9", textColor = "black")
#forestPlot(metaObject = Prostate_metaanalysis, geneName = "GNPTAB", textColor = "black")
#forestPlot(metaObject = Prostate_metaanalysis, geneName = "IQGAP3", textColor = "black")
#forestPlot(metaObject = Prostate_metaanalysis, geneName = "RPRML", textColor = "black")
#forestPlot(metaObject = Prostate_metaanalysis, geneName = "NCOA2", textColor = "black")
dev.off()

# For the down-regulated genes
pdf(file = "./Figs/Meta/DownGenes_effect_size.pdf", title = "Effect size of Down-regulated genes across the discovery data sets", width = 20, height = 12)
par(mfrow=c(3,4))

forestPlot(metaObject = Prostate_metaanalysis, geneName = "AZGP1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "CPA3", textColor = "black")
#forestPlot(metaObject = Prostate_metaanalysis, geneName = "BRE", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "UFM1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "DPT", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "PART1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "CBLL1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "LUZP2", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "CHRNA2", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "WNT8B", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "NT5E", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "DCXR", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "EDN3", textColor = "black")
#forestPlot(metaObject = Prostate_metaanalysis, geneName = "NT5DC1", textColor = "black")
#forestPlot(metaObject = Prostate_metaanalysis, geneName = "COL4A5", textColor = "black")

dev.off()

## Heatmap of the effect sizes of the signature genes

png("./Figs/Meta/Heatmap_EffectSizes.png", width = 4000, height = 2000, res = 300)
heatmapPlot(Prostate_metaanalysis, filterObject = New_filter)
dev.off()

###########################################################################################
#############################################################################

## The next step is validation on indepndent data set (GSE116918, and GSE16560)

# Get the validation data set
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

#################################
## Modify val_expr1
rownames(val_expr1) <- valDataSet1$keys
summary(is.na(rownames(val_expr1)))
val_expr1 <- val_expr1[!is.na(rownames(val_expr1)), ]
dim(val_expr1)

#val_expr1 <- val_expr1[commonGenes, ]

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
## val_Pheno2 is already processed

## Label the samples
valDataSet1 <- classFunction(valDataSet1, column = "Metastasis", diseaseTerms = c("Mets"))
valDataSet2 <- classFunction(valDataSet2, column = "Metastasis", diseaseTerms = c("Mets"))
valDataSet3 <- classFunction(valDataSet3, column = "Metastasis", diseaseTerms = c("Mets"))
valDataSet4 <- classFunction(valDataSet4, column = "Metastasis", diseaseTerms = c("Mets"))
#valDataSet5 <- classFunction(valDataSet5, column = "Metastasis", diseaseTerms = c("Mets"))
#valDataSet6 <- classFunction(valDataSet6, column = "Metastasis", diseaseTerms = c("Mets"))

##############################################

### Now lets examine the performance of our filter on the independent data set

#### Using ROC curve

# Validation dataset1
png(filename = "./Figs/Meta/ROC_Test1.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = valDataSet1, filterObject = filter)
dev.off()

# Validation dataset2
valDataSet2$keys <- rownames(val_expr2)

png(filename = "./Figs/Meta/ROC_Test2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = valDataSet2, filterObject = filter)
dev.off()

# Validation dataset3
valDataSet3$keys <- rownames(val_expr3)

png(filename = "./Figs/Meta/ROC_Test3.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = valDataSet3, filterObject = filter)
dev.off()

# Validation dataset4
valDataSet4$keys <- rownames(val_expr4)

png(filename = "./Figs/Meta/ROC_Test4.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = valDataSet4, filterObject = filter)
dev.off()

# Validation dataset5
#valDataSet5$keys <- rownames(val_expr5)

#png(filename = "./Figs/Meta/ROC_Test5.png", width = 2000, height = 2000, res = 300)
#rocPlot(datasetObject = valDataSet5, filterObject = New_filter)
#dev.off()

# Validation dataset6
#valDataSet6$keys <- rownames(val_expr6)

#png(filename = "./Figs/Meta/ROC_Test6.png", width = 2000, height = 2000, res = 300)
#rocPlot(datasetObject = valDataSet6, filterObject = New_filter)
#dev.off()


################################

#### Using PRC plot

# val dataset1
prcPlot(datasetObject = valDataSet1, filterObject = New_filter)

# val dataset2
prcPlot(datasetObject = valDataSet2, filterObject =filter)

# val dataset3
prcPlot(datasetObject = valDataSet3, filterObject =filter)

# val dataset4
prcPlot(datasetObject = valDataSet4, filterObject =filter)

# val dataset5
#prcPlot(datasetObject = valDataSet5, filterObject =New_filter)

# val dataset6
#prcPlot(datasetObject = valDataSet6, filterObject =New_filter)

##################################

#### Using violin plot

# Val Dataset1 
png(filename = "./Figs/Meta/ViolinPlot_Test1.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = filter, datasetObject = valDataSet1, labelColumn = "Metastasis")
dev.off()

# Val Dataset2
png(filename = "./Figs/Meta/ViolinPlot_Test2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = filter, datasetObject = valDataSet2, labelColumn = "Metastasis")
dev.off()

# Val Dataset3
png(filename = "./Figs/Meta/ViolinPlot_Test3.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = filter, datasetObject = valDataSet3, labelColumn = "Metastasis")
dev.off()

# Val Dataset4
png(filename = "./Figs/Meta/ViolinPlot_Test4.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = filter, datasetObject = valDataSet4, labelColumn = "Metastasis")
dev.off()

# Val Dataset5
#png(filename = "./Figs/Meta/ViolinPlot_Test5.png", width = 2000, height = 2000, res = 300)
#violinPlot(filterObject = New_filter, datasetObject = valDataSet5, labelColumn = "Metastasis")
#dev.off()

# Val Dataset6
#png(filename = "./Figs/Meta/ViolinPlot_Test6.png", width = 2000, height = 2000, res = 300)
#violinPlot(filterObject = New_filter, datasetObject = valDataSet6, labelColumn = "Metastasis")
#dev.off()

###################################################
# pdf(file = "./Figs/Meta/Multi.pdf", width = 10, height = 10)
# par(mfrow = c(2,2))
# summaryROCPlot(metaObject = Prostate_metaanalysis, filterObject = New_filter, alphaBetaPlots = FALSE)
# rocPlot(datasetObject = valDataSet, filterObject = New_filter)
# violinPlot(filterObject = New_filter, datasetObject = valDataSet, labelColumn = "Metastasis")
# dev.off()
# #####################################################

##### Calculate ROC data

# Val Dataset1
#scoreRs1 <- calculateScore(filterObject = New_filter, datasetObject = valDataSet1)
#rocRes1 <- calculateROC(predictions = scoreRs1, labels = valDataSet1$class)

# Val Dataset2
#scoreRs2 <- calculateScore(filterObject = New_filter, datasetObject = valDataSet2)
#rocRes2 <- calculateROC(predictions = scoreRs2, labels = valDataSet2$class)

##############
#### Calculate a signature score (Z score) and add it to the phenotype table

# Val Dataset1
val_pheno1$score <- calculateScore(filterObject = filter, datasetObject = valDataSet1)

# Val Dataset2
val_pheno2$score <- calculateScore(filterObject = filter, datasetObject = valDataSet2)

# Val Dataset3
val_pheno3$score <- calculateScore(filterObject = filter, datasetObject = valDataSet3)

# Val Dataset4
val_pheno4$score <- calculateScore(filterObject = filter, datasetObject = valDataSet4)

# Val Dataset5
#val_pheno5$score <- calculateScore(filterObject = New_filter, datasetObject = valDataSet5)

# Val Dataset6
#val_pheno6$score <- calculateScore(filterObject = New_filter, datasetObject = valDataSet6)

## save val_pheno1 with z-score for futher survival analysis
save(val_pheno1, val_pheno2, val_pheno3, val_pheno4, PositiveGenes, NegativeGenes, thr_test1, thr_test2, thr_test3, thr_test4, file = "./Objs/val_pheno1_Zscore.rda")
#########################################

## ROC curve in the Testing data using pROC

#png(filename = "./Figs/Meta/ROC_Test_pROC.png", width = 2000, height = 2000, res = 300)
#roc(val_pheno1$Metastasis, val_pheno1$score, plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), col="blue", lwd=2, grid=TRUE, auc = TRUE, ci = TRUE)
#dev.off()

#### Find the best threshold (for further use for the classification in the Test data)

# Val Dataset1
thr_test1 <- coords(roc(val_pheno1$Metastasis, val_pheno1$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
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

#### Find the optimal trade off between the sensitivity and specificity (We want good sensitivity with preservation of a decent accuracy)

# Val Dataset1
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

#############################
##### Predictions

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

################################
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

Meta_Test <- geneSymbolCorrection(metaObject = Meta_Test)

#set.seed(333)
#summaryROCPlot(metaObject = Meta_Test, filterObject = filter)

set.seed(333)
png(filename = "./Figs/Meta/PooledROC_InitialSig_Testing.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = Meta_Test, filterObject = filter)
dev.off()

set.seed(333)
png(filename = "./Figs/Meta/PooledROC_RefinedSig_Testing.png", width = 2000, height = 2000, res = 300)
pooledROCPlot(metaObject = Meta_Test, filterObject = New_filter)
dev.off()

###############################################################################
## Gene set enrichment analysis
library(enrichR)
dbs <- listEnrichrDbs()
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018" , "ChEA_2016" ,"KEGG_2019_Human")
Enriched_PositiveGns <- enrichr(genes = PositiveGenes, databases = dbs)
Enriched_NegativeGns <- enrichr(genes = NegativeGenes, databases = dbs)
printEnrich(Enriched_PositiveGns, "PosGnsEnrichment.txt" , sep = "\t", columns = c(1:9))
printEnrich(Enriched_NegativeGns, "NegGnsEnrichment.txt" , sep = "\t", columns = c(1:9))

Pos_GO_BP <- Enriched_PositiveGns["GO_Biological_Process_2018"]
Pos_GO_BP <- Pos_GO_BP$GO_Biological_Process_2018
Pos_GO_BP <- Pos_GO_BP[Pos_GO_BP$P.value <= 0.05, ]

Neg_GO_BP <- Enriched_NegativeGns["GO_Biological_Process_2018"]
Neg_GO_BP <- Neg_GO_BP$GO_Biological_Process_2018
Neg_GO_BP <- Neg_GO_BP[Neg_GO_BP$P.value <= 0.05, ]


Pos_KEGG <- Enriched_PositiveGns["KEGG_2019_Human"]
Pos_KEGG <- Pos_KEGG$KEGG_2019_Human
Pos_KEGG <- Pos_KEGG[Pos_KEGG$P.value <= 0.05, ]

Neg_KEGG <- Enriched_NegativeGns["KEGG_2019_Human"]
Neg_KEGG <- Neg_KEGG$KEGG_2019_Human
Neg_KEGG <- Neg_KEGG[Neg_KEGG$P.value <= 0.05, ]

########################################################################
## Estimating the immune cell proportions based on the gene expression profiles
immuno_meta <- immunoStatesMeta(metaObject = Prostate_metaanalysis)

## Correct the underlying gene expression data based on the cell proportions
immuno_meta_corrected <- immunoStatesDecov(metaObject = Prostate_metaanalysis)

#########################################################################
## Using lincsCorrelate to look for a drug with a gene expression profile that reverses the metastasis profile
lincProstate <- lincsCorrelate(metaObject = Prostate_metaanalysis, filterObject = New_filter, dataset = "CP", direction = "reverse")

####
val_expr1 <- t(scale(t(val_expr1), center = T, scale = T))
val_expr2 <- t(scale(t(val_expr2), center = T, scale = T))
val_expr3 <- t(scale(t(val_expr3), center = T, scale = T))
val_expr4 <- t(scale(t(val_expr4), center = T, scale = T))
val_expr5 <- t(scale(t(val_expr5), center = T, scale = T))
val_expr6 <- t(scale(t(val_expr6), center = T, scale = T))



AllValExpr <- list(val_expr1, val_expr2, val_expr3, val_expr6)
CommonGns <- Reduce("intersect", lapply(AllValExpr, rownames))

AllValExpr <- mapply(x=AllValExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=CommonGns))



AllValExpr <- do.call("cbind", AllValExpr)
AllValExpr <- normalizeBetweenArrays(AllValExpr, method = "quantile")


val_pheno1 <- as.data.frame(val_pheno1[,"Metastasis"])
val_pheno2 <- as.data.frame(val_pheno2[,"Metastasis"])
val_pheno3 <- as.data.frame(val_pheno3[,"Metastasis"])
val_pheno4 <- as.data.frame(val_pheno4[,"Metastasis"])
val_pheno5 <- as.data.frame(val_pheno5[,"Metastasis"])
val_pheno6 <- as.data.frame(val_pheno6[,"Metastasis"])

colnames(val_pheno1) <- "Metastasis"
colnames(val_pheno2) <- "Metastasis"
colnames(val_pheno3) <- "Metastasis"
colnames(val_pheno4) <- "Metastasis"
colnames(val_pheno5) <- "Metastasis"
colnames(val_pheno6) <- "Metastasis"


AllValPheno <- rbind(val_pheno1, val_pheno2, val_pheno3, val_pheno6)

all(rownames(AllValPheno) == colnames(AllValExpr))
rownames(AllValPheno) <- colnames(AllValExpr)

ValDataSet <- list()
ValDataSet$expr <- AllValExpr
ValDataSet$pheno <- AllValPheno
ValDataSet$keys <- rownames(AllValExpr)
ValDataSet$formattedName <- "ValDataSet"

ValDataSet <- classFunction(ValDataSet, column = "Metastasis", diseaseTerms = c("Mets"))

##############################################

### Now lets examine the performance of our filter on the independent data set

rocPlot(datasetObject = ValDataSet, filterObject = OncoDx_Filter)

X <- getSampleLevelGeneData(datasetObject = valDataSet1, geneNames = c(PositiveGenes))

X <- cleanUpPheno(valDataSet1)