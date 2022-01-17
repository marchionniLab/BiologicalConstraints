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

#load("./Data/Dataset7.rda")
#######################################################

## Getting the phenotype data for each data set
#pheno4 <- pData(Dataset4)
#ProstateData$originalData$GSE46691$pheno <- pheno4
pheno2 <- Dataset2$pheno
pheno3 <- Dataset3$pheno
pheno4 <- Dataset4$pheno
Pheno5 <- Dataset5$pheno
Pheno6 <- Dataset6$pheno
#Pheno7 <- Dataset7$pheno
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
#expr7 <- Dataset7$expr
## Checking if the expression data are normalized and log2 transformed
# boxplot(expr2[,1:15], outline= FALSE)
# boxplot(expr3[,1:15], outline = FALSE)
# boxplot(expr4[,1:15], outline= FALSE)

#################################################################
## Create a list containing training data sets
AllDataSets <- list(Dataset2, Dataset3, Dataset4, Dataset5)
names(AllDataSets) <- c("GSE55935", "GSE51066", "GSE46691", "GSE41408")

###################################################################

## Annotate expression

## Expr1
head(rownames(expr1))
rownames(expr1) <- Dataset1$keys
expr1 <- expr1[!is.na(rownames(expr1)), ]
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
## Expr7
# head(rownames(expr7))
# rownames(expr7) <- Dataset7$keys
# summary(is.na(rownames(expr7)))
# expr7 <- expr7[!is.na(rownames(expr7)), ]
# dim(expr7)

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
# Pheno7$Metastasis <- Pheno7$`pathology stage:ch1`
# 
# Pheno7$Metastasis <- gsub("T.+M", "M", Pheno7$Metastasis)
# Pheno7$Metastasis <- gsub("N.+", "", Pheno7$Metastasis)
# 
# table(Pheno7$Metastasis)
# 
# Pheno7 <- Pheno7[!(Pheno7$Metastasis == "Mx"), ]
# 
# Pheno7 <- Pheno7[!(Pheno7$Metastasis == "pMx"), ]
# Pheno7 <- Pheno7[!(Pheno7$Metastasis == "U"), ]
# 
# Pheno7$Metastasis <- as.factor(Pheno7$Metastasis)
# levels(Pheno7$Metastasis) <- c("No_Mets", "No_Mets", "Mets", "Mets")
# table(Pheno7$Metastasis)
# 
# expr7 <- expr7[, colnames(expr7) %in% rownames(Pheno7)]
# all(rownames(Pheno7) == colnames(expr7))
# 
# Dataset7$pheno <- Pheno7
# Dataset7$expr <- expr7
# Dataset7$keys <- rownames(expr7)
#########################################################################
############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘disease’ or 0 for ‘control’)
Dataset2 <- classFunction(Dataset2, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset3 <- classFunction(Dataset3, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset4 <- classFunction(Dataset4, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset5 <- classFunction(Dataset5, column = "Metastasis", diseaseTerms = c("Mets"))
#Dataset7 <- classFunction(Dataset7, column = "Metastasis", diseaseTerms = c("Mets"))

############################################################################
############################################################################
#########################################################################
## The metaanalysis

## Creating the meta object
AllDataSets <- list(Dataset2, Dataset3, Dataset4, Dataset5)
names(AllDataSets) <- c("GSE55935", "GSE51066", "GSE46691", "GSE41408")

Prostate_meta <- list()
Prostate_meta$originalData <- AllDataSets

## Replace keys within each data set
Prostate_meta$originalData$GSE55935$keys <- rownames(expr2)
Prostate_meta$originalData$GSE51066$keys <- rownames(expr3)
Prostate_meta$originalData$GSE46691$keys <- rownames(expr4)
Prostate_meta$originalData$GSE41408$keys <- rownames(expr5)

## Check the meta object before the metaanalysis
checkDataObject(Prostate_meta, "Meta", "Pre-Analysis") ## If true, Proceed to the meta analysis


## Run the meta analysis
Prostate_metaanalysis <- runMetaAnalysis(Prostate_meta, runLeaveOneOutAnalysis = F, maxCores = 3)

## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
Prostate_metaanalysis <- filterGenes(Prostate_metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.15, numberStudiesThresh = 4, heterogeneityPvalThresh = 0.05)

## Assigning a name to the filter
filter <- Prostate_metaanalysis$filterResults[[1]]

PositiveGns <- filter$posGeneNames
NegativeGns <- filter$negGeneNames
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
Prostate_metaanalysis$filterResults$FDR0.15_es0_nStudies4_looaFALSE_hetero0.05$posGeneNames <- New_filter$posGeneNames
Prostate_metaanalysis$filterResults$FDR0.15_es0_nStudies4_looaFALSE_hetero0.05$negGeneNames <- New_filter$negGeneNames

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
png(filename = "./Figs/Meta/Summary_ROC_TrainingDatasets.png", width = 2000, height = 2000, res = 300)
summaryROCPlot(metaObject = Prostate_metaanalysis, filterObject = New_filter)
dev.off()

## Effect size visualization in each study
# For the up-regulated genes
pdf(file = "./Figs/Meta/UpGenes_effect_size.pdf", title = "Effect size of up-regulated genes across the discovery data sets", width = 20, height = 10)
par(mfrow=c(2,4))

forestPlot(metaObject = Prostate_metaanalysis, geneName = "TMSB10", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "ASPN", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "P4HA2", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "PTDSS1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "GNPTAB", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "IQGAP3", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "RPRML", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "NCOA2", textColor = "black")
dev.off()

# For the down-regulated genes
pdf(file = "./Figs/Meta/DownGenes_effect_size.pdf", title = "Effect size of Down-regulated genes across the discovery data sets", width = 20, height = 12)
par(mfrow=c(3,5))

forestPlot(metaObject = Prostate_metaanalysis, geneName = "AZGP1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "CPA3", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "BRE", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "UFM1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "DPT", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "PART1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "CBLL1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "SIDT1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "CHRNA2", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "WNT8B", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "NT5E", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "LTBR", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "EDN3", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "NT5DC1", textColor = "black")
forestPlot(metaObject = Prostate_metaanalysis, geneName = "COL4A5", textColor = "black")

dev.off()

## Heatmap of the effect sizes of the signature genes
png("./Figs/Meta/Heatmap_EffectSizes.png", width = 4000, height = 2000, res = 300)
heatmapPlot(Prostate_metaanalysis, filterObject = New_filter)
dev.off()

###########################################################################################
#############################################################################

## The next step is validation on indepndent data set (GSE116918, and GSE16560)

# Get the validation data set
#valDataSet <- getGEOData("GSE116918")
valDataSet1 <- Dataset1
valDataSet2 <- Dataset6
#valDataSet3 <- Dataset7

# Acess the phenotype and expression data
val_pheno1 <- Dataset1$pheno
val_expr1 <- Dataset1$expr

val_pheno2 <- Pheno6
val_expr2 <- expr6

#val_pheno3 <- Pheno7
#val_expr3 <- expr7

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
#valDataSet3 <- classFunction(valDataSet3, column = "Metastasis", diseaseTerms = c("Mets"))

##############################################

### Now lets examine the performance of our filter on the independent data set

#### Using ROC curve

# Validation dataset1
png(filename = "./Figs/Meta/ROC_Test1.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = valDataSet1, filterObject = New_filter)
dev.off()

# Validation dataset2
valDataSet2$keys <- rownames(val_expr2)

png(filename = "./Figs/Meta/ROC_Test2.png", width = 2000, height = 2000, res = 300)
rocPlot(datasetObject = valDataSet2, filterObject = New_filter)
dev.off()

# Validation dataset3
#valDataSet3$keys <- rownames(val_expr3)

#png(filename = "./Figs/Meta/ROC_Test2.png", width = 2000, height = 2000, res = 300)
#rocPlot(datasetObject = valDataSet3, filterObject = New_filter)
#dev.off()

################################

#### Using PRC plot

# val dataset1
prcPlot(datasetObject = valDataSet1, filterObject = New_filter)

# val dataset2
prcPlot(datasetObject = valDataSet2, filterObject =New_filter)

# val dataset3
prcPlot(datasetObject = valDataSet3, filterObject =New_filter)
##################################

#### Using violin plot

# Val Dataset1 
png(filename = "./Figs/Meta/ViolinPlot_Test1.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter, datasetObject = valDataSet1, labelColumn = "Metastasis")
dev.off()

# Val Dataset2
png(filename = "./Figs/Meta/ViolinPlot_Test2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter, datasetObject = valDataSet2, labelColumn = "Metastasis")
dev.off()

# Val Dataset3
png(filename = "./Figs/Meta/ViolinPlot_Test2.png", width = 2000, height = 2000, res = 300)
violinPlot(filterObject = New_filter, datasetObject = valDataSet3, labelColumn = "Metastasis")
dev.off()


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
scoreRs1 <- calculateScore(filterObject = New_filter, datasetObject = valDataSet1)
rocRes1 <- calculateROC(predictions = scoreRs1, labels = valDataSet1$class)

# Val Dataset2
scoreRs2 <- calculateScore(filterObject = New_filter, datasetObject = valDataSet2)
rocRes2 <- calculateROC(predictions = scoreRs2, labels = valDataSet2$class)

##############
#### Calculate a signature score (Z score) and add it to the phenotype table

# Val Dataset1
val_pheno1$score <- calculateScore(filterObject = New_filter, datasetObject = valDataSet1)

# Val Dataset2
val_pheno2$score <- calculateScore(filterObject = New_filter, datasetObject = valDataSet2)

## save val_pheno1 with z-score for futher survival analysis
save(val_pheno1, val_pheno2, PositiveGenes, NegativeGenes, thr_test1, thr_test2, file = "./Objs/val_pheno1_Zscore.rda")
#########################################

## ROC curve in the Testing data using pROC

png(filename = "./Figs/Meta/ROC_Test_pROC.png", width = 2000, height = 2000, res = 300)
roc(val_pheno1$Metastasis, val_pheno1$score, plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), col="blue", lwd=2, grid=TRUE, auc = TRUE, ci = TRUE)
dev.off()

#### Find the best threshold (for further use for the classification in the Test data)

# Val Dataset1
thr_test1 <- coords(roc(val_pheno1$Metastasis, val_pheno1$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
thr_test1 

# Val Dataset2
thr_test2 <- coords(roc(val_pheno2$Metastasis, val_pheno2$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
thr_test2 

#### Find the optimal trade off between the sensitivity and specificity (We want good sensitivity with preservation of a decent accuracy)

# Val Dataset1
coords(roc(val_pheno1$Metastasis, val_pheno1$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")

# Val Dataset2
coords(roc(val_pheno2$Metastasis, val_pheno2$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")

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

################################

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

