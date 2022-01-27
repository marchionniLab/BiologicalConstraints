############################################################################
rm(list = ls())


library(GEOquery)
library(Biobase)
library(sampling)
library(limma)
library(genefilter)
library(edgeR)
library(caret)
library(sampling)
library(org.Hs.eg.db)
library(AnnotationDbi)


### Getting the data 
# dataset_GSE78220 <- getGEO("GSE78220", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset_GSE78220 <- dataset_GSE78220$GSE78220_series_matrix.txt.gz
# 
# #dataset_GSE93157 <- getGEO("GSE93157", GSEMatrix = TRUE, AnnotGPL = TRUE)
# #dataset_GSE93157 <- dataset_GSE93157$GSE93157_series_matrix.txt.gz
# 
# dataset_GSE91061 <- getGEO("GSE91061", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset_GSE91061 <- dataset_GSE91061$GSE91061_series_matrix.txt.gz
# 
# dataset_GSE115821 <- getGEO("GSE115821", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset_GSE115821_GPL11154 <- dataset_GSE115821$`GSE115821-GPL11154_series_matrix.txt.gz`
# dataset_GSE115821_GPL18573 <- dataset_GSE115821$`GSE115821-GPL18573_series_matrix.txt.gz`
# 
# save(dataset_GSE78220, dataset_GSE93157, dataset_GSE91061, dataset_GSE115821_GPL11154, dataset_GSE115821_GPL18573, file = "./Data/datasets.rda")

load("./Data/datasets.rda")
###################################################

## Get the phenotype
Pheno_GSE78220 <- pData(dataset_GSE78220)
#Pheno_GSE93157 <- pData(dataset_GSE93157)
Pheno_GSE91061 <- pData(dataset_GSE91061)
Pheno_GSE115821_GPL11154 <- pData(dataset_GSE115821_GPL11154)
Pheno_GSE115821_GPL18573 <- pData(dataset_GSE115821_GPL18573)
Pheno_GSE115821 <- rbind(Pheno_GSE115821_GPL11154, Pheno_GSE115821_GPL18573)
Pheno_TCGA <- read.delim("./Data/skcm_mskcc_2014/skcm_mskcc_2014_clinical_data.tsv")
Pheno_VanAllen <- read.delim("./Data/VanAllen/skcm_dfci_2015_clinical_data.tsv")

## Get the expression matrices
Expr_GSE78220 <- read.delim("./Data/HUGO/HUGO_Expr.csv", header = TRUE, row.names = 1, sep = ",")
#Expr_GSE93157 <- exprs(dataset_GSE93157)
Expr_GSE91061 <- read.delim("./Data/RIAZ_Expr2.txt", header = TRUE, row.names = 1, sep = "")
Expr_GSE115821 <- read.delim("./Data/MGH_Expr.csv", header = TRUE, sep = ",")
Expr_TCGA <-  read.delim("./Data/skcm_mskcc_2014/data_RNA_Seq_expression_median.txt")
Expr_VanAllen <-  read.delim("./Data/VanAllen/data_RNA_Seq_expression_median.txt")

## Get the feature data
#FeatData_GSE93157 <- fData(dataset_GSE93157)

##############################################################
## Process expression matrices

################
# Process Expr_GSE78220

# Remove genes with zero expression
keep <- rowSums(Expr_GSE78220) > 0
table(keep)

Expr_GSE78220 <- Expr_GSE78220[keep,]

# Modify the column names to be consistent with the sample names in the phenotype
colnames(Expr_GSE78220) <- gsub("\\..+", "", colnames(Expr_GSE78220))

# remove miRNA
MiRNA <- grep("^MIR", rownames(Expr_GSE78220))
Expr_GSE78220 <- Expr_GSE78220[-MiRNA, ]

# filtering
# mycpm <- cpm(Expr_GSE78220)
# thresh <- mycpm > 1
# keep <- rowSums(thresh) >= round(ncol(Expr_GSE78220)/2)
# table(keep)
# 
# Expr_GSE78220 <- Expr_GSE78220[keep,]
# dim(Expr_GSE78220)

# Normalize
#Expr_GSE78220 <- DGEList(Expr_GSE78220)
#Expr_GSE78220 <- calcNormFactors(Expr_GSE78220, method = c("TMM"))
#Expr_GSE78220 <- cpm(Expr_GSE78220, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# remove genes with missing values and empty gene symbols
sel <- which(apply(Expr_GSE78220, 1, function(x) all(is.finite(x)) ))
Expr_GSE78220 <- Expr_GSE78220[sel, ]
Expr_GSE78220 <- Expr_GSE78220[!is.na(rownames(Expr_GSE78220)),]
Expr_GSE78220 <- Expr_GSE78220[!(rownames(Expr_GSE78220) == ""), ] 
dim(Expr_GSE78220)

## Filtering
# ffun <- kOverA(round(ncol(Expr_GSE78220)*0.2), 10)
# X <- Expr_GSE78220
# Filt <- genefilter(X, ffun)
# summary(Filt)
# Expr_GSE78220 <- Expr_GSE78220[Filt, ]
# dim(Expr_GSE78220)
# 
# # Log2-transform
range(Expr_GSE78220)
Expr_GSE78220 <- log2(Expr_GSE78220 + 1)


Expr_GSE78220 <- t(scale(t(Expr_GSE78220), center = TRUE, scale = TRUE))

################
# Process Expr_GSE91061

# Mapping Entrez ID to gene symbols
tmp <- as.character(rownames(Expr_GSE91061))
Expr_GSE91061$GeneSymbol<- mapIds(org.Hs.eg.db,
                                  keys=tmp,
                                  column="SYMBOL",
                                  keytype="ENTREZID",
                                  multiVals="first")

Expr_GSE91061 <- Expr_GSE91061[!duplicated(Expr_GSE91061$GeneSymbol), ] 
Expr_GSE91061 <- Expr_GSE91061[!is.na(Expr_GSE91061$GeneSymbol), ] 
rownames(Expr_GSE91061) <- Expr_GSE91061$GeneSymbol
Expr_GSE91061$GeneSymbol <- NULL
Expr_GSE91061$HUGO <- NULL

H19 <- grep("^H19", rownames(Expr_GSE91061))
Expr_GSE91061 <- Expr_GSE91061[-H19, ]
MiR <- grep("^MIR", rownames(Expr_GSE91061))
Expr_GSE91061 <- Expr_GSE91061[-MiR, ]

# filtering
# mycpm <- cpm(Expr_GSE91061)
# thresh <- mycpm > 1
# keep <- rowSums(thresh) >= round(ncol(Expr_GSE91061)/2)
# table(keep)

#Expr_GSE91061 <- Expr_GSE91061[keep,]
#dim(Expr_GSE91061)

# Normalize
#Expr_GSE91061 <- DGEList(Expr_GSE91061)
#Expr_GSE91061 <- calcNormFactors(Expr_GSE91061, method = c("TMM"))
#Expr_GSE91061 <- cpm(Expr_GSE91061, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# remove genes with missing values and empty gene symbols
sel <- which(apply(Expr_GSE91061, 1, function(x) all(is.finite(x)) ))
Expr_GSE91061 <- Expr_GSE91061[sel, ]
Expr_GSE91061 <- Expr_GSE91061[!is.na(rownames(Expr_GSE91061)), ]
Expr_GSE91061 <- Expr_GSE91061[!(rownames(Expr_GSE91061) == ""), ] 
dim(Expr_GSE91061)

# ## Filtering
# ffun <- kOverA(round(ncol(Expr_GSE91061)*0.2), 200)
# X <- Expr_GSE91061
# Filt <- genefilter(X, ffun)
# summary(Filt)
# Expr_GSE91061 <- Expr_GSE91061[Filt, ]
# dim(Expr_GSE91061)
# 
# # Log2_transform
range(Expr_GSE91061)
Expr_GSE91061 <- log2(Expr_GSE91061 + 1)

Expr_GSE91061 <- t(scale(t(Expr_GSE91061), center = TRUE, scale = TRUE))

################
# Process Expr_GSE115821

# Process and normalize
Expr_GSE115821 <- Expr_GSE115821[!duplicated(Expr_GSE115821$Geneid), ]
rownames(Expr_GSE115821) <- Expr_GSE115821$Geneid
Expr_GSE115821 <- Expr_GSE115821[, -c(1:6)]

miR <- grep("^MIR", rownames(Expr_GSE115821))
Expr_GSE115821 <- Expr_GSE115821[-miR, ]

# filtering
# mycpm <- cpm(Expr_GSE115821)
# thresh <- mycpm > 1
# keep <- rowSums(thresh) >= round(ncol(Expr_GSE115821)/2)
# table(keep)
# 
# Expr_GSE115821 <- Expr_GSE115821[keep,]
# dim(Expr_GSE115821)

# Normalize
# Expr_GSE115821 <- DGEList(Expr_GSE115821)
# Expr_GSE115821 <- calcNormFactors(Expr_GSE115821, method = c("TMM"))
# Expr_GSE115821 <- cpm(Expr_GSE115821, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# # Log2_transform
range(Expr_GSE115821)
Expr_GSE115821 <- log2(Expr_GSE115821 + 1)

# z-score
Expr_GSE115821 <- t(scale(t(Expr_GSE115821), center = TRUE, scale = TRUE))

#boxplot(Expr_GSE115821)

## MAke the sample names consistent 
colnames(Expr_GSE115821) <- gsub("X", "", colnames(Expr_GSE115821))
colnames(Expr_GSE115821) <- gsub(".bam", "", colnames(Expr_GSE115821))
colnames(Expr_GSE115821) <- gsub("\\.", "-", colnames(Expr_GSE115821))

# same for pheno
rownames(Pheno_GSE115821) <- Pheno_GSE115821$title
rownames(Pheno_GSE115821) <- gsub(".bam", "", rownames(Pheno_GSE115821))

#####################
## Process Expr_TCGA

# Mapping the ENTREZ IDs to Gene Symbols
tmp <- as.character(Expr_TCGA$Entrez_Gene_Id)
Expr_TCGA$GeneSymbol<- mapIds(org.Hs.eg.db,
                              keys=tmp,
                              column="SYMBOL",
                              keytype="ENTREZID",
                              multiVals="first")

Expr_TCGA <- Expr_TCGA[!duplicated(Expr_TCGA$GeneSymbol), ] 
Expr_TCGA <- Expr_TCGA[!is.na(Expr_TCGA$GeneSymbol), ] 
rownames(Expr_TCGA) <- Expr_TCGA$GeneSymbol
Expr_TCGA$Entrez_Gene_Id <- NULL
Expr_TCGA$GeneSymbol <- NULL

H19 <- grep("^H19", rownames(Expr_TCGA))
Expr_TCGA <- Expr_TCGA[-H19, ]

# filtering
# mycpm <- cpm(Expr_TCGA)
# thresh <- mycpm > 1
# keep <- rowSums(thresh) >= round(ncol(Expr_TCGA)/2)
# table(keep)
# 
# Expr_TCGA <- Expr_TCGA[keep,]
# dim(Expr_TCGA)
# 
# # Normalize
# Expr_TCGA <- DGEList(Expr_TCGA)
# Expr_TCGA <- calcNormFactors(Expr_TCGA, method = c("TMM"))
# Expr_TCGA <- cpm(Expr_TCGA, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# Log2transform
range(Expr_TCGA)
Expr_TCGA <- log2(Expr_TCGA+1)

## Filtering
# ffun <- pOverA(p = 0.1, A = 10)
# X <- Expr_TCGA
# Filt <- genefilter(2^X, ffun)
# summary(Filt)
# Expr_TCGA <- Expr_TCGA[Filt, ]
# dim(Expr_TCGA)


Expr_TCGA <- t(scale(t(Expr_TCGA), center = TRUE, scale = TRUE))

#######################
## Process Expr_VanAllen

# Mapping the ENTREZ IDs to Gene Symbols
tmp <- as.character(Expr_VanAllen$Entrez_Gene_Id)
Expr_VanAllen$GeneSymbol<- mapIds(org.Hs.eg.db,
                                  keys=tmp,
                                  column="SYMBOL",
                                  keytype="ENTREZID",
                                  multiVals="first")

Expr_VanAllen <- Expr_VanAllen[!duplicated(Expr_VanAllen$GeneSymbol), ] 
Expr_VanAllen <- Expr_VanAllen[!is.na(Expr_VanAllen$GeneSymbol), ] 
rownames(Expr_VanAllen) <- Expr_VanAllen$GeneSymbol
Expr_VanAllen$Entrez_Gene_Id <- NULL
Expr_VanAllen$GeneSymbol <- NULL

H19 <- grep("^H19", rownames(Expr_VanAllen))
Expr_VanAllen <- Expr_VanAllen[-H19, ]
MiR <- grep("^MIR", rownames(Expr_VanAllen))
Expr_VanAllen <- Expr_VanAllen[-MiR, ]

rownames(Expr_VanAllen) <- gsub("-","", rownames(Expr_VanAllen))
rownames(Expr_VanAllen) <- gsub("_","",rownames(Expr_VanAllen))
rownames(Expr_VanAllen) <- gsub("\\.","",rownames(Expr_VanAllen))

keep <- rowSums(Expr_VanAllen) > 0
table(keep)

Expr_VanAllen <- Expr_VanAllen[keep,]
dim(Expr_VanAllen)

# filtering
#mycpm <- cpm(Expr_VanAllen)
#thresh <- mycpm > 1
#keep <- rowSums(thresh) >= round(ncol(Expr_VanAllen)/2)
#table(keep)

#Expr_VanAllen <- Expr_VanAllen[keep,]
#dim(Expr_VanAllen)

# Normalize
#Expr_VanAllen <- DGEList(Expr_VanAllen)
#Expr_VanAllen <- calcNormFactors(Expr_VanAllen, method = c("TMM"))
#Expr_VanAllen <- cpm(Expr_VanAllen, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# Log2transform
range(Expr_VanAllen)
Expr_VanAllen <- log2(Expr_VanAllen + 1)
 
# ## Filtering
# ffun <- pOverA(p = 0.1, A = 10)
# X <- Expr_VanAllen
# Filt <- genefilter(2^X, ffun)
# summary(Filt)
# Expr_VanAllen <- Expr_VanAllen[Filt, ]
# dim(Expr_VanAllen)

Expr_VanAllen <- t(scale(t(Expr_VanAllen), center = TRUE, scale = TRUE))

###########################################
### Modify the phenotypes

## Modify Pheno_GSE78220

rownames(Pheno_GSE78220) <- Pheno_GSE78220$title

# Check for consistency
all(rownames(Pheno_GSE78220) == colnames(Expr_GSE78220))

# Keep only pre-treatment samples
Pheno_GSE78220 <- Pheno_GSE78220[Pheno_GSE78220$`biopsy time:ch1` == 'pre-treatment', ]

# Get the response variable and convert to factor
Pheno_GSE78220$Response <- as.factor(Pheno_GSE78220$`anti-pd-1 response:ch1`)
table(Pheno_GSE78220$Response)
levels(Pheno_GSE78220$Response) <- c("R", "R", "NR")
Pheno_GSE78220$Response <- factor(Pheno_GSE78220$Response, levels = c("NR", "R"))

# subset the expression
Expr_GSE78220 <- Expr_GSE78220[, colnames(Expr_GSE78220) %in% rownames(Pheno_GSE78220)]
# Check for consistency
all(rownames(Pheno_GSE78220) == colnames(Expr_GSE78220))

#####################
## Modify Pheno_GSE91061

# Modify the rownames to match the expression column names
Pheno_GSE91061$title <- gsub('_A.+', '', Pheno_GSE91061$title)
Pheno_GSE91061$title <- gsub('_E.+', '', Pheno_GSE91061$title)
Pheno_GSE91061 <- Pheno_GSE91061[!(Pheno_GSE91061$geo_accession == 'GSM2420320'), ]

rownames(Pheno_GSE91061) <- Pheno_GSE91061$title

# Keep only pre-treatment samples
Pheno_GSE91061 <- Pheno_GSE91061[Pheno_GSE91061$`visit (pre or on treatment):ch1` == 'Pre', ]

# remove samples with unknown response
Pheno_GSE91061 <- Pheno_GSE91061[!(Pheno_GSE91061$`response:ch1` == 'UNK'), ]

# Get the response variable
Pheno_GSE91061$Response <- as.factor(Pheno_GSE91061$`response:ch1`)
table(Pheno_GSE91061$Response)

# Re-name the levels 
levels(Pheno_GSE91061$Response) <- c("NR", "R", "NR") 
Pheno_GSE91061$Response <- factor(Pheno_GSE91061$Response, levels = c("NR", "R"))

#Remove unwanted samples from the expression
Keep <- intersect(rownames(Pheno_GSE91061), colnames(Expr_GSE91061))
Expr_GSE91061 <- Expr_GSE91061[, Keep]
# Check if sample names are identical
all(rownames(Pheno_GSE91061) == colnames(Expr_GSE91061))

#####################
## Modify Pheno_GSE115821

rownames(Pheno_GSE115821) <- Pheno_GSE115821$title
rownames(Pheno_GSE115821) <- gsub(".bam", "", rownames(Pheno_GSE115821))

## Reorder based on rownames and colnames
Expr_GSE115821 <- Expr_GSE115821[, order(colnames(Expr_GSE115821))]
Pheno_GSE115821 <- Pheno_GSE115821[order(rownames(Pheno_GSE115821)), ]

# Check for consistency
all(rownames(Pheno_GSE115821) == colnames(Expr_GSE115821))

# Keep only pre-treatment samples
Pheno_GSE115821 <- Pheno_GSE115821[Pheno_GSE115821$`treatment state:ch1` %in% c('PRE Immune Checkpoint Blockade Therapy', 'PRE Immune Checkpoint Blockade Therapy (On dabrafenib+trametinib)'), ]

# Remove duplicate samples
Pheno_GSE115821 <- Pheno_GSE115821[!(Pheno_GSE115821$geo_accession %in% c('GSM3190470', 'GSM3190484', 'GSM3190485', 'GSM3190487', 'GSM3190494', 'GSM3190495')), ]

# Get the response variable
Pheno_GSE115821$Response <- as.factor(Pheno_GSE115821$`response:ch1`)
table(Pheno_GSE115821$Response)
levels(Pheno_GSE115821$Response) <- c("NR", "R")

# subset the expression
Expr_GSE115821 <- Expr_GSE115821[, colnames(Expr_GSE115821) %in% rownames(Pheno_GSE115821)]

# Check for consistency
all(rownames(Pheno_GSE115821) == colnames(Expr_GSE115821))

#####################
## Modify Pheno_TCGA

# Remove the samples with no corresponding RNA-Seq
rownames(Pheno_TCGA) <- Pheno_TCGA$Sample.ID
Keep <- intersect(rownames(Pheno_TCGA), colnames(Expr_TCGA))
Pheno_TCGA <- Pheno_TCGA[Keep, ]

# Check consistency
all(rownames(Pheno_TCGA) == colnames(Expr_TCGA))

# Keep only pre-treatment samples
Pheno_TCGA <- Pheno_TCGA[Pheno_TCGA$Biopsy.Time == 'pre', ]
Pheno_TCGA <- Pheno_TCGA[!is.na(Pheno_TCGA$Treatment.Response), ]

# Get the response variable
Pheno_TCGA$Response <- as.factor(Pheno_TCGA$Treatment.Response)
table(Pheno_TCGA$Response)
levels(Pheno_TCGA$Response) <- c("NR", "R")

# subset the expression
Expr_TCGA <- Expr_TCGA[, colnames(Expr_TCGA) %in% rownames(Pheno_TCGA)]

# Check for consistency
all(rownames(Pheno_TCGA) == colnames(Expr_TCGA))

#####################
## Modify Pheno_VanAllen

# Remove the samples with no corresponding RNA-Seq
rownames(Pheno_VanAllen) <- Pheno_VanAllen$Sample.ID
Keep <- intersect(rownames(Pheno_VanAllen), colnames(Expr_VanAllen))
Pheno_VanAllen <- Pheno_VanAllen[Keep, ]

# Check consistency
all(rownames(Pheno_VanAllen) == colnames(Expr_VanAllen))

# all samples are pre-treatment

# Get the response variable
Pheno_VanAllen$Response <- as.factor(Pheno_VanAllen$Durable.Clinical.Benefit)
table(Pheno_VanAllen$Response)
levels(Pheno_VanAllen$Response) <- c("R", "NR", "R", "R", NA)
Pheno_VanAllen$Response <- factor(Pheno_VanAllen$Response, levels = c("NR", "R"))
Pheno_VanAllen <- Pheno_VanAllen[!is.na(Pheno_VanAllen$Response), ]

# subset the expression
Expr_VanAllen <- Expr_VanAllen[, colnames(Expr_VanAllen) %in% rownames(Pheno_VanAllen)]

# Check for consistency
all(rownames(Pheno_VanAllen) == colnames(Expr_VanAllen))

#########################################################################################
#########################################################################################
## Aggregate

allPheno <- list(Pheno_GSE78220, Pheno_GSE91061, Pheno_GSE115821, Pheno_TCGA, Pheno_VanAllen)
names(allPheno) <- c("GSE78220", "GSE91061", "GSE115821", "TCGA", "VanAllen")

ICB_Response_ALL <- list(Pheno_GSE78220$Response, Pheno_GSE91061$Response, Pheno_GSE115821$Response, Pheno_TCGA$Response, Pheno_VanAllen$Response)
names(ICB_Response_ALL) <- c("GSE78220", "GSE91061", "GSE115821", "TCGA", "VanAllen")

allExpr <- list(Expr_GSE78220, Expr_GSE91061, Expr_GSE115821, Expr_TCGA, Expr_VanAllen)

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(allExpr, rownames))

### Filter expression for the required samples
exprsICB <- mapply(x=allExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

names(exprsICB) <- c("GSE78220", "GSE91061", "GSE115821", "TCGA", "VanAllen")

### Check
all(names(exprsICB) == names(ICB_Response_ALL))

### Check order
all(rownames(allPheno$GSE78220) == colnames(exprsICB$GSE78220))
all(rownames(allPheno$GSE91061) == colnames(exprsICB$GSE91061))
all(rownames(allPheno$GSE115821) == colnames(exprsICB$GSE115821))
all(rownames(allPheno$TCGA) == colnames(exprsICB$TCGA))
all(rownames(allPheno$VanAllen) == colnames(exprsICB$VanAllen))

###################################################################
#############
## Cross-study validation

# Leave GSE78220 out
train <- c("GSE91061", "GSE115821", "TCGA", "VanAllen")
test <- c("GSE78220")

## Training
trainMat <- do.call("cbind", exprsICB[train])
trainGroup <- factor(do.call("c", ICB_Response_ALL[train]))
levels(trainGroup) <- c("NR", "R")
table(trainGroup)

## Testing
testMat <- exprsICB$GSE78220
testGroup <- factor(do.call("c", ICB_Response_ALL[test]))
levels(testGroup) <- c("NR", "R")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/icbData_GSE78220Out_pre.rda")

##########################
# Leave GSE91061 out
train <- c("GSE78220", "GSE115821", "TCGA", "VanAllen")
test <- c("GSE91061")

## Training
trainMat <- do.call("cbind", exprsICB[train])
trainGroup <- factor(do.call("c", ICB_Response_ALL[train]))
levels(trainGroup) <- c("NR", "R")
table(trainGroup)

## Testing
testMat <- exprsICB$GSE91061
testGroup <- factor(do.call("c", ICB_Response_ALL[test]))
levels(testGroup) <- c("NR", "R")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/icbData_GSE91061Out_pre.rda")

##########################
# Leave GSE115821 out
train <- c("GSE78220", "GSE91061", "TCGA", "VanAllen")
test <- c("GSE115821")

## Training
trainMat <- do.call("cbind", exprsICB[train])
trainGroup <- factor(do.call("c", ICB_Response_ALL[train]))
levels(trainGroup) <- c("NR", "R")
table(trainGroup)

## Testing
testMat <- exprsICB$GSE115821
testGroup <- factor(do.call("c", ICB_Response_ALL[test]))
levels(testGroup) <- c("NR", "R")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/icbData_GSE115821Out_pre.rda")

##########################
# Leave TCGA out
train <- c("GSE78220", "GSE91061", "GSE115821", "VanAllen")
test <- c("TCGA")

## Training
trainMat <- do.call("cbind", exprsICB[train])
trainGroup <- factor(do.call("c", ICB_Response_ALL[train]))
levels(trainGroup) <- c("NR", "R")
table(trainGroup)

## Testing
testMat <- exprsICB$TCGA
testGroup <- factor(do.call("c", ICB_Response_ALL[test]))
levels(testGroup) <- c("NR", "R")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/icbData_TCGAOut_pre.rda")

##########################
# Leave VanAllen out
train <- c("GSE78220", "GSE91061", "GSE115821", "TCGA")
test <- c("VanAllen")

## Training
trainMat <- do.call("cbind", exprsICB[train])
trainGroup <- factor(do.call("c", ICB_Response_ALL[train]))
levels(trainGroup) <- c("NR", "R")
table(trainGroup)

## Testing
testMat <- exprsICB$VanAllen
testGroup <- factor(do.call("c", ICB_Response_ALL[test]))
levels(testGroup) <- c("NR", "R")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/icbData_VanAllenOut_pre.rda")

##################################################################################
#####################
## All combined
allMat <- do.call("cbind", exprsICB)
allGroup <- unlist(ICB_Response_ALL)
allStudies <- gsub("\\..+", "", names(allGroup))

names(allGroup) <- colnames(allMat)
all(colnames(allMat) == names(allGroup))

####################################################################################
### Covariates for stratified sampling

#############################################################
### Pre or on

# allPheno$GSE115821$BiopsyTiming <- as.factor(allPheno$GSE115821$`treatment state:ch1`)
# levels(allPheno$GSE115821$BiopsyTiming) <- c('on', 'pre', 'pre')
# 
# allPheno$GSE91061$BiopsyTiming <- as.factor(allPheno$GSE91061$`visit (pre or on treatment):ch1`)
# levels(allPheno$GSE91061$BiopsyTiming) <- c('on', 'pre')
# 
# allPheno$GSE78220$BiopsyTiming <- as.factor(allPheno$GSE78220$`biopsy time:ch1`)
# levels(allPheno$GSE78220$BiopsyTiming) <- c('on', 'pre')
# 
# allPheno$TCGA$BiopsyTiming <- as.factor(allPheno$TCGA$Biopsy.Time)
# levels(allPheno$TCGA$BiopsyTiming) <- c('pre')
# 
# allPheno$VanAllen$BiopsyTiming <- as.factor(rep('pre', nrow(allPheno$VanAllen)))
# 
# ### Covariates of relevance select complete cases: GRADE
# allBioposyTiming <- lapply(allPheno, function(x) {
#   i <- grep("BiopsyTiming", colnames(x))
#   if (length(i) == 0) out <- factor(rep("NA", nrow(x)))
#   else x <- factor(x[, i  ])
# })
# 
# allBioposyTiming <- factor(unlist(allBioposyTiming))
# table(allBioposyTiming)

#############################################################
### gender

# Pheno_GSE91061: no gender info
# Pheno_GSE115821: no gender info

allPheno$GSE78220$gender <- as.factor(allPheno$GSE78220$`gender:ch1`)
allPheno$GSE78220$`gender:ch1` <- NULL
levels(allPheno$GSE78220$gender)

allPheno$TCGA$gender <- as.factor(allPheno$TCGA$Sex)
levels(allPheno$TCGA$gender) <- c('F', 'M')

allPheno$VanAllen$gender <- as.factor(allPheno$VanAllen$Sex)
levels(allPheno$VanAllen$gender) <- c('F', 'M')

### Covariates of relevance select complete cases: GRADE
allGender <- lapply(allPheno, function(x) {
  i <- grep("gender", colnames(x))
  if (length(i) == 0) out <- factor(rep("NA", nrow(x)))
  else x <- factor(x[, i  ])
})

allGender <- factor(unlist(allGender))
table(allGender)

#################################################################################
###############################
### Age

# Pheno_GSE91061: no age info

allPheno$GSE78220$AGE <- as.numeric(allPheno$GSE78220$`age (yrs):ch1`)
allPheno$GSE115821$AGE <- as.numeric(allPheno$GSE115821$`age at the baseline:ch1`)
allPheno$TCGA$AGE <- as.numeric(allPheno$TCGA$Age.at.Diagnosis)
allPheno$VanAllen$AGE <- as.numeric(allPheno$VanAllen$Age.at.Diagnosis)


### Covariates of relevance select complete cases: AGE
allAGE <- lapply(allPheno, function(x) {
  i <- grep("AGE", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})

allAGE <- unlist(allAGE)
summary(allAGE)

#########################################################################

### Assemble in one data.frame and turn numeric
covs <- data.frame(STUDIES=allStudies,
                   GENDER=allGender,
                   AGE=allAGE)

### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))) )

###########################################################################
###SAMPLING

### Balanced stratification
set.seed(333)
trainingOrTesting <- balancedstratification(
  covs[ , , drop=FALSE], strata=1*(allGroup == "R"),
  pik=inclusionprobabilities(1:nrow(covs), nrow(covs) * 0.25),
  comment=TRUE, method=1)

### Show
apply(covs[, -ncol(covs),drop=FALSE], 2, table, allGroup, trainingOrTesting)

### Subset Training
mixTrainMat <- allMat[ , trainingOrTesting == 0]
mixTrainGroup <- allGroup[ trainingOrTesting == 0]
mixTrainStudy <- allStudies[ trainingOrTesting == 0]

### Subset Testing
mixTestMat <- allMat[ , trainingOrTesting == 1]
mixTestGroup <- allGroup[ trainingOrTesting == 1]
mixTestStudy <- allStudies[ trainingOrTesting == 1]

table(mixTrainGroup)
table(mixTestGroup)

###########################################################################
### Save
save(exprsICB, 
     mixTrainMat, mixTrainGroup, mixTrainStudy,
     mixTestMat, mixTestGroup, mixTestStudy,
     file="./Objs/icbData_pre.rda")

