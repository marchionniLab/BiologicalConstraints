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
library(illuminaHumanv4.db)

### Getting the data 
# dataset_GSE78220 <- getGEO("GSE78220", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset_GSE78220 <- dataset_GSE78220$GSE78220_series_matrix.txt.gz
# 
# #dataset_GSE93157 <- getGEO("GSE93157", GSEMatrix = TRUE, AnnotGPL = TRUE)
# #dataset_GSE93157 <- dataset_GSE93157$GSE93157_series_matrix.txt.gz
# 
# dataset_GSE91061 <- getGEO("GSE91061", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset_GSE91061 <- dataset_GSE91061$GSE91061_series_matrix.txt.gz

# dataset_GSE115821 <- getGEO("GSE115821", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset_GSE115821_GPL11154 <- dataset_GSE115821$`GSE115821-GPL11154_series_matrix.txt.gz`
# dataset_GSE115821_GPL18573 <- dataset_GSE115821$`GSE115821-GPL18573_series_matrix.txt.gz`

#dataset_GSE122220 <- getGEO('GSE122220', GSEMatrix = TRUE, AnnotGPL = TRUE)
#dataset_GSE122220 <- dataset_GSE122220$GSE122220_series_matrix.txt.gz

#save(dataset_GSE78220, dataset_GSE93157, dataset_GSE91061, dataset_GSE115821_GPL11154, dataset_GSE115821_GPL18573, dataset_GSE122220, file = "./Data/datasets.rda")

load("./Data/datasets.rda")

load("./Data/NP_Dataset.rda")

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
Pheno_GSE122220 <- pData(dataset_GSE122220)
# load the pheno of the NP dataset
Pheno_NP <- pData(NP_Dataset)

## Get the expression matrices
Expr_GSE78220 <- read.delim("./Data/HUGO/HUGO_Expr.csv", header = TRUE, row.names = 1, sep = ",")
#Expr_GSE93157 <- exprs(dataset_GSE93157)
Expr_GSE91061 <- read.delim("./Data/RIAZ_Expr2.txt", header = TRUE, row.names = 1, sep = "")
Expr_GSE115821 <- read.delim("./Data/MGH_Expr.csv", header = TRUE, sep = ",")
Expr_TCGA <-  read.delim("./Data/skcm_mskcc_2014/data_RNA_Seq_expression_median.txt")
Expr_VanAllen <-  read.delim("./Data/VanAllen/data_RNA_Seq_expression_median.txt")
Expr_GSE122220 <- exprs(dataset_GSE122220)
# read the expression of NP dataset
Expr_NP <- read.delim("./Data/NP_Expr/GSE49711_SEQC_NB_TUC_G_log2.txt")

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

#Expr_GSE78220 <- Expr_GSE78220[keep,]
#dim(Expr_GSE78220)

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

# Normalize
# Expr_GSE78220 <- DGEList(Expr_GSE78220)
# Expr_GSE78220 <- calcNormFactors(Expr_GSE78220, method = c("TMM"))
# Expr_GSE78220 <- cpm(Expr_GSE78220, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

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
#mycpm <- cpm(Expr_GSE91061)
#thresh <- mycpm > 1
#keep <- rowSums(thresh) >= round(ncol(Expr_GSE91061)/2)
#table(keep)

# Expr_GSE91061 <- Expr_GSE91061[keep,]
# dim(Expr_GSE91061)

# Normalize
# Expr_GSE91061 <- DGEList(Expr_GSE91061)
# Expr_GSE91061 <- calcNormFactors(Expr_GSE91061, method = c("TMM"))
# Expr_GSE91061 <- cpm(Expr_GSE91061, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# remove genes with missing values and empty gene symbols
sel <- which(apply(Expr_GSE91061, 1, function(x) all(is.finite(x)) ))
Expr_GSE91061 <- Expr_GSE91061[sel, ]
Expr_GSE91061 <- Expr_GSE91061[!is.na(rownames(Expr_GSE91061)), ]
Expr_GSE91061 <- Expr_GSE91061[!(rownames(Expr_GSE91061) == ""), ] 
dim(Expr_GSE91061)

# Normalize
# Expr_GSE91061 <- DGEList(Expr_GSE91061)
# Expr_GSE91061 <- calcNormFactors(Expr_GSE91061, method = c("TMM"))
# Expr_GSE91061 <- cpm(Expr_GSE91061, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

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

#Expr_GSE115821 <- Expr_GSE115821[keep,]
#dim(Expr_GSE115821)

# Normalize
Expr_GSE115821 <- DGEList(Expr_GSE115821)
Expr_GSE115821 <- calcNormFactors(Expr_GSE115821, method = c("TMM"))
Expr_GSE115821 <- cpm(Expr_GSE115821, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

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

# Normalize
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
# mycpm <- cpm(Expr_VanAllen)
# thresh <- mycpm > 1
# keep <- rowSums(thresh) >= round(ncol(Expr_VanAllen)/2)
# table(keep)
# 
# Expr_VanAllen <- Expr_VanAllen[keep,]
# dim(Expr_VanAllen)

# Normalize
# Expr_VanAllen <- DGEList(Expr_VanAllen)
# Expr_VanAllen <- calcNormFactors(Expr_VanAllen, method = c("TMM"))
# Expr_VanAllen <- cpm(Expr_VanAllen, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

# Log2transform
range(Expr_VanAllen)
Expr_VanAllen <- log2(Expr_VanAllen + 1)
# 
# ## Filtering
# ffun <- pOverA(p = 0.1, A = 10)
# X <- Expr_VanAllen
# Filt <- genefilter(2^X, ffun)
# summary(Filt)
# Expr_VanAllen <- Expr_VanAllen[Filt, ]
# dim(Expr_VanAllen)

Expr_VanAllen <- t(scale(t(Expr_VanAllen), center = TRUE, scale = TRUE))

#######################
## Process Expr_NP

Expr_NP$X00gene_id <- gsub("-", "", Expr_NP$X00gene_id)
# Remove duplicate and missing Gene symbols
Expr_NP <- Expr_NP[!duplicated(Expr_NP$X00gene_id),]
Expr_NP <- Expr_NP[!is.na(Expr_NP$X00gene_id), ]

# Change the rownames to the gene symbol
rownames(Expr_NP) <- Expr_NP$X00gene_id

# Remove columns with useless gene information
Expr_NP$X00gene_id <- NULL

# Convert to matrix
Expr_NP <- as.matrix(Expr_NP)
dim(Expr_NP)
range(Expr_NP)

# Remove Genes with missing expression
sel <- which(apply(Expr_NP, 1, function(x) all(is.finite(x)) ))
Expr_NP <- Expr_NP[sel, ]

Expr_NP <- t(scale(t(Expr_NP), center = TRUE, scale = TRUE))

##############################################
###########################################
### Modify the phenotypes

## Modify Pheno_GSE78220

rownames(Pheno_GSE78220) <- Pheno_GSE78220$title

# Check for consistency
all(rownames(Pheno_GSE78220) == colnames(Expr_GSE78220))

# Keep only pre-treatment samples
#Pheno_GSE78220 <- Pheno_GSE78220[Pheno_GSE78220$`biopsy time:ch1` == 'pre-treatment', ]

# Get the response variable and convert to factor
Pheno_GSE78220$Response <- as.factor(Pheno_GSE78220$`anti-pd-1 response:ch1`)
table(Pheno_GSE78220$Response)
levels(Pheno_GSE78220$Response) <- c("R", "R", "NR")
Pheno_GSE78220$Response <- ordered(Pheno_GSE78220$Response, levels = c("NR", "R"))

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
#Pheno_GSE91061 <- Pheno_GSE91061[Pheno_GSE91061$`visit (pre or on treatment):ch1` == 'Pre', ]

# remove samples with unknown response
Pheno_GSE91061 <- Pheno_GSE91061[!(Pheno_GSE91061$`response:ch1` == 'UNK'), ]

# Check if sample names are identical
all(rownames(Pheno_GSE91061) == colnames(Expr_GSE91061))

# Get the response variable
Pheno_GSE91061$Response <- as.factor(Pheno_GSE91061$`response:ch1`)
table(Pheno_GSE91061$Response)

# Re-name the levels 
levels(Pheno_GSE91061$Response) <- c("NR", "R", "R") 
Pheno_GSE91061$Response <- ordered(Pheno_GSE91061$Response, levels = c("NR", "R"))

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
#Pheno_GSE115821 <- Pheno_GSE115821[Pheno_GSE115821$`treatment state:ch1` %in% c('PRE Immune Checkpoint Blockade Therapy', 'PRE Immune Checkpoint Blockade Therapy (On dabrafenib+trametinib)'), ]

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
Pheno_VanAllen$Response <- ordered(Pheno_VanAllen$Response, levels = c("NR", "R"))
Pheno_VanAllen <- Pheno_VanAllen[!is.na(Pheno_VanAllen$Response), ]

# subset the expression
Expr_VanAllen <- Expr_VanAllen[, colnames(Expr_VanAllen) %in% rownames(Pheno_VanAllen)]

# Check for consistency
all(rownames(Pheno_VanAllen) == colnames(Expr_VanAllen))

###########################
## Modify Pheno_NP

rownames(Pheno_NP) <- Pheno_NP$title

colnames(Expr_NP) <- rownames(Pheno_NP)

# Keep only samples where high_risk = Progression (Either 0 or 1)
Keep <- which(Pheno_NP$`progression:ch1` == Pheno_NP$`high risk:ch1`)

# Subset the phenotype table
Pheno_NP <- Pheno_NP[Keep, ]

Pheno_NP$SR <- as.factor(Pheno_NP$`progression:ch1`)

Pheno_NP$`age at diagnosis:ch1` <- as.numeric(Pheno_NP$`age at diagnosis:ch1`)
Pheno_NP <- Pheno_NP[!(Pheno_NP$`age at diagnosis:ch1` == 0), ] 
Pheno_NP$Age_Month <- Pheno_NP$`age at diagnosis:ch1`/30

Pheno_NP <- Pheno_NP[Pheno_NP$Age_Month < 18, ] 

## Convert to factor
table(Pheno_NP$SR)
Pheno_NP$SR <- factor(Pheno_NP$SR, levels = c(1,0))
levels(Pheno_NP$SR) <- c("NR", "R")  
Pheno_NP$SR <- factor(Pheno_NP$SR, levels = c('NR','R'))


# Keep the 156 samples from RUppen
load("./Data/NP_SampleNames.rda")
KeepSampls <- intersect(rownames(Pheno_NP), SampleNames)
Pheno_NP <- Pheno_NP[KeepSampls, ]

## Remove filtered samples from the expression
Expr_NP <- Expr_NP[, KeepSampls]

# Check consistency
all(rownames(Pheno_NP) == colnames(Expr_NP))

#########################################################################################
#########################################################################################
## Aggregate

allPheno_test <- list(Pheno_GSE78220, Pheno_GSE91061, Pheno_GSE115821, Pheno_TCGA, Pheno_VanAllen)
names(allPheno_test) <- c("GSE78220", "GSE91061", "GSE115821", "TCGA", "VanAllen")

ICB_Response_test <- list(Pheno_GSE78220$Response, Pheno_GSE91061$Response, Pheno_GSE115821$Response, Pheno_TCGA$Response, Pheno_VanAllen$Response)
names(ICB_Response_test) <- c("GSE78220", "GSE91061", "GSE115821", "TCGA", "VanAllen")

allExpr_test <- list(Expr_GSE78220, Expr_GSE91061, Expr_GSE115821, Expr_TCGA, Expr_VanAllen)

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(allExpr_test, rownames))

### Find commom subset of genes
commonGenes2 <- intersect(commonGenes, rownames(Expr_NP))

### Filter expression for the required samples
exprsICB_test <- mapply(x=allExpr_test, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes2))

names(exprsICB_test) <- c("GSE78220", "GSE91061", "GSE115821", "TCGA", "VanAllen")

# also subset the np dataset
Expr_NP <- Expr_NP[commonGenes2, ]

### Check
all(names(exprsICB_test) == names(ICB_Response_test))

### Check order
all(rownames(allPheno_test$GSE78220) == colnames(exprsICB_test$GSE78220))
all(rownames(allPheno_test$GSE91061) == colnames(exprsICB_test$GSE91061))
all(rownames(allPheno_test$GSE115821) == colnames(exprsICB_test$GSE115821))
all(rownames(allPheno_test$TCGA) == colnames(exprsICB_test$TCGA))
all(rownames(allPheno_test$VanAllen) == colnames(exprsICB_test$VanAllen))

##################################################################################
#####################
# training: NP
trainMat <- Expr_NP
trainGroup <- Pheno_NP$SR

##combined (Testing)
testMat <- do.call("cbind", exprsICB_test)
testGroup <- unlist(ICB_Response_test)

names(testGroup) <- colnames(testMat)
all(colnames(testMat) == names(testGroup))


###########################################################################
### Save
save(trainMat, trainGroup,
     testMat, testGroup,
     file="./Objs/icbData_withNP.rda")

