
rm(list = ls())

setwd("/Volumes/Macintosh/Research/Projects/ICB_Response")

library(GEOquery)
library(Biobase)
library(sampling)
library(limma)
library(genefilter)
library(edgeR)
library(sampling)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicFeatures)
require(reshape)
library(annotate)

########################################3

# Get the NP dataset

# NP_Dataset <- getGEO("GSE49711", GSEMatrix = TRUE, AnnotGPL = TRUE)
# NP_Dataset <- NP_Dataset$GSE49711_series_matrix.txt.gz
# save(NP_Dataset, file = "./Data/NP_Dataset.rda")

load("./Data/NP_Dataset.rda")

# read the expression
Expr_NP <- read.delim("./Data/NP_Expr/GSE49711_SEQC_NB_TUC_G_log2.txt")

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
boxplot(Expr_NP[,1:50])

rownames(Expr_NP) <- gsub("-","", rownames(Expr_NP))
rownames(Expr_NP) <- gsub("_","",rownames(Expr_NP))
rownames(Expr_NP) <- gsub("\\.","",rownames(Expr_NP))

# Remove Genes with missing expression
sel <- which(apply(Expr_NP, 1, function(x) all(is.finite(x)) ))
Expr_NP <- Expr_NP[sel, ]

###############

### Get the pheno
Pheno_NP <- pData(NP_Dataset)

rownames(Pheno_NP) <- colnames(Expr_NP)
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
levels(Pheno_NP$SR) <- c("no", "yes")  
  
# Change the rownames (Sample IDs) to match those in the Expression matrix
#rownames(Pheno_NP) <- Pheno_NP$title
  
## Remove filtered samples from the expression
keepSamples <- rownames(Pheno_NP)
Expr_NP <- Expr_NP[, keepSamples]

# Check consistency
all(rownames(Pheno_NP) == colnames(Expr_NP))

# convert to Z-score
keep <- rowSums(Expr_NP) > 0
table(keep)

Expr_NP <- Expr_NP[keep,]
dim(Expr_NP)

#Expr_NP <- t(scale(t(Expr_NP), center = T, scale = T))
###################################################################
###################################################################
# Get the Chen dataset

Chen <- read.delim("./Data/Chen.csv", header = TRUE, sep = ",", row.names = 1)


# Separate into 3 different datasets (Pre_aCTLA4, Pre_aPD1, and On_aPD1)
Chen_Pre_aCTLA4 <- Chen[1:16, ]
Chen_Pre_aPD1 <- Chen[22:37, ] 
Chen_On_aPD1 <- Chen[38:47, ]

####################################################
### Saprate each one into phenotype and expression

## Chen_Pre_aCTLA4
Pheno_Chen_Pre_aCTLA4 <- Chen_Pre_aCTLA4[, 1:4]
Expr_Chen_Pre_aCTLA4 <- t(Chen_Pre_aCTLA4[, 5:length(colnames(Chen_Pre_aCTLA4))])

rownames(Expr_Chen_Pre_aCTLA4) <- gsub("-","", rownames(Expr_Chen_Pre_aCTLA4))
rownames(Expr_Chen_Pre_aCTLA4) <- gsub("_","",rownames(Expr_Chen_Pre_aCTLA4))
rownames(Expr_Chen_Pre_aCTLA4) <- gsub("\\.","",rownames(Expr_Chen_Pre_aCTLA4))

# Check consistent Sample IDs
all(rownames(Pheno_Chen_Pre_aCTLA4) == colnames(Expr_Chen_Pre_aCTLA4))

# convert to Z-score
#Expr_Chen_Pre_aCTLA4 <- t(scale(t(Expr_Chen_Pre_aCTLA4), center = T, scale = T))

#######
## Chen_Pre_aPD1
Pheno_Chen_Pre_aPD1 <- Chen_Pre_aPD1[, 1:4]
Expr_Chen_Pre_aPD1 <- t(Chen_Pre_aPD1[, 5:length(colnames(Chen_Pre_aPD1))])
rownames(Expr_Chen_Pre_aPD1) <- gsub("-","", rownames(Expr_Chen_Pre_aPD1))
rownames(Expr_Chen_Pre_aPD1) <- gsub("_","",rownames(Expr_Chen_Pre_aPD1))
rownames(Expr_Chen_Pre_aPD1) <- gsub("\\.","",rownames(Expr_Chen_Pre_aPD1))

all(rownames(Pheno_Chen_Pre_aPD1) == colnames(Expr_Chen_Pre_aPD1))

# convert to Z-score
#Expr_Chen_Pre_aPD1 <- t(scale(t(Expr_Chen_Pre_aPD1), center = T, scale = T))

#######
## Chen_On_aPD1
Pheno_Chen_On_aPD1 <- Chen_On_aPD1[, 1:4]
Expr_Chen_On_aPD1 <- t(Chen_On_aPD1[, 5:length(colnames(Chen_On_aPD1))])
rownames(Expr_Chen_On_aPD1) <- gsub("-","", rownames(Expr_Chen_On_aPD1))
rownames(Expr_Chen_On_aPD1) <- gsub("_","",rownames(Expr_Chen_On_aPD1))
rownames(Expr_Chen_On_aPD1) <- gsub("\\.","",rownames(Expr_Chen_On_aPD1))

all(rownames(Pheno_Chen_On_aPD1) == colnames(Expr_Chen_On_aPD1))

# convert to Z-score
#Expr_Chen_On_aPD1 <- t(scale(t(Expr_Chen_On_aPD1), center = T, scale = T))

#######################################################
### Process the phenotypes

## Chen_Pre_aCTLA4 
Pheno_Chen_Pre_aCTLA4$Response <- as.factor(Pheno_Chen_Pre_aCTLA4$anti.CTLA.4.response)
table(Pheno_Chen_Pre_aCTLA4$Response)
levels(Pheno_Chen_Pre_aCTLA4$Response) <- c(NA,"NR", "R")
# Get the class of interest
Group_Chen_Pre_aCTLA4 <- Pheno_Chen_Pre_aCTLA4$Response
names(Group_Chen_Pre_aCTLA4) <- rownames(Pheno_Chen_Pre_aCTLA4)

########
## Chen_Pre_aPD1 
Pheno_Chen_Pre_aPD1$Response <- as.factor(Pheno_Chen_Pre_aPD1$anti.PD.1.response)
table(Pheno_Chen_Pre_aPD1$Response)
levels(Pheno_Chen_Pre_aPD1$Response) <- c(NA, "NR", "R")
# Get the class of interest
Group_Chen_Pre_aPD1 <- Pheno_Chen_Pre_aPD1$Response
names(Group_Chen_Pre_aPD1) <- rownames(Pheno_Chen_Pre_aPD1)

########
## Chen_On_aPD1
Pheno_Chen_On_aPD1$Response <- as.factor(Pheno_Chen_On_aPD1$anti.PD.1.response)
table(Pheno_Chen_On_aPD1$Response)
levels(Pheno_Chen_On_aPD1$Response) <- c(NA, "NR", "R")
# Get the class of interest
Group_Chen_On_aPD1 <- Pheno_Chen_On_aPD1$Response
names(Group_Chen_On_aPD1) <- rownames(Pheno_Chen_On_aPD1)

###############################################################
###############################################################
## Get the HUGO Dataset (All is : Anti-PD-1)

# Hugo_Dataset <- getGEO("GSE78220", GSEMatrix = TRUE, AnnotGPL = TRUE)
# Hugo_Dataset <- Hugo_Dataset$GSE78220_series_matrix.txt.gz
# save(Hugo_Dataset, file = "./Data/HUGO_Dataset.rda")

load("./Data/HUGO_Dataset.rda")

# Load the expression file
Expr_HUGO <- read.delim("./Data/HUGO/HUGO_Expr.csv", header = TRUE, row.names = 1, sep = ",")

# Remove genes with zero expression
keep <- rowSums(Expr_HUGO) > 0
table(keep)

Expr_HUGO <- Expr_HUGO[keep,]
dim(Expr_HUGO)

# Log2-transform
Expr_HUGO <- log2(Expr_HUGO + 1)

# Modify the column names to be consistent with the sample names in the phenotype
colnames(Expr_HUGO) <- gsub("\\..+", "", colnames(Expr_HUGO))

# Modify gene symbols
# remove miRNA
MiRNA <- grep("^MIR", rownames(Expr_HUGO))
Expr_HUGO <- Expr_HUGO[-MiRNA, ]

rownames(Expr_HUGO) <- gsub("-","", rownames(Expr_HUGO))
rownames(Expr_HUGO) <- gsub("_","",rownames(Expr_HUGO))
rownames(Expr_HUGO) <- gsub("\\.","",rownames(Expr_HUGO))


# Get the phenotype
Pheno_HUGO <- pData(Hugo_Dataset)
rownames(Pheno_HUGO) <- Pheno_HUGO$title

# Check for consistency
all(rownames(Pheno_HUGO) == colnames(Expr_HUGO))

# Get the response variable and convert to factor
Pheno_HUGO$Response <- as.factor(Pheno_HUGO$`anti-pd-1 response:ch1`)
table(Pheno_HUGO$Response)
levels(Pheno_HUGO$Response) <- c("R", "R", "NR")
Pheno_HUGO$Response <- ordered(Pheno_HUGO$Response, levels = c("NR", "R"))

# Get the class of interest
Group_HUGO <- Pheno_HUGO$Response  
names(Group_HUGO) <- rownames(Pheno_HUGO)  
  
# convert to Z-score
#Expr_HUGO <- t(scale(t(Expr_HUGO), center = T, scale = T))

###################################################################
###################################################################
## Get the TCGA dataset

## Get the Expr
Expr_TCGA <-  read.delim("./Data/skcm_mskcc_2014/data_RNA_Seq_expression_median.txt")

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
rownames(Expr_TCGA) <- gsub("-","", rownames(Expr_TCGA))
rownames(Expr_TCGA) <- gsub("_","",rownames(Expr_TCGA))
rownames(Expr_TCGA) <- gsub("\\.","",rownames(Expr_TCGA))

keep <- rowSums(Expr_TCGA) > 0
table(keep)

Expr_TCGA <- Expr_TCGA[keep,]
dim(Expr_TCGA)

# Log2transform
Expr_TCGA <- log2(Expr_TCGA+1)

## Get the Pheno
Pheno_TCGA <- read.delim("./Data/skcm_mskcc_2014/skcm_mskcc_2014_clinical_data.tsv")

# Remove the samples with no corresponding RNA-Seq
rownames(Pheno_TCGA) <- Pheno_TCGA$Sample.ID
Keep <- intersect(rownames(Pheno_TCGA), colnames(Expr_TCGA))
Pheno_TCGA <- Pheno_TCGA[Keep, ]

# Check consistency
all(rownames(Pheno_TCGA) == colnames(Expr_TCGA))

# Get the response variable
Pheno_TCGA$Response <- as.factor(Pheno_TCGA$Treatment.Response)
table(Pheno_TCGA$Response)
levels(Pheno_TCGA$Response) <- c("NR", "R")

# Remove missing
Pheno_TCGA <- Pheno_TCGA[!is.na(Pheno_TCGA$Response), ]

# the same for the expression
Expr_TCGA <- Expr_TCGA[,colnames(Expr_TCGA) %in% rownames(Pheno_TCGA)]

# Check consistency
all(rownames(Pheno_TCGA) == colnames(Expr_TCGA))

# Get the class of interest
Group_TCGA <- Pheno_TCGA$Response  
names(Group_TCGA) <- rownames(Pheno_TCGA)  

# convert to Z-score
#Expr_TCGA <- t(scale(t(Expr_TCGA), center = T, scale = T))

###################################################################
###################################################################
## Get the PRAT Dataset

# PRAT_Dataset <- getGEO("GSE93157", GSEMatrix = TRUE, AnnotGPL = TRUE)
# PRAT_Dataset <- PRAT_Dataset$GSE93157_series_matrix.txt.gz
# save(PRAT_Dataset, file = "./Data/PRAT_Dataset.rda")

load("./Data/PRAT_Dataset.rda")

## Get the Expression
Expr_PRAT <- exprs(PRAT_Dataset)

# Remove Genes with missing expression
sel <- which(apply(Expr_PRAT, 1, function(x) all(is.finite(x)) ))
Expr_PRAT <- Expr_PRAT[sel, ]

rownames(Expr_PRAT) <- gsub("-","", rownames(Expr_PRAT))
rownames(Expr_PRAT) <- gsub("_","",rownames(Expr_PRAT))
rownames(Expr_PRAT) <- gsub("\\.","",rownames(Expr_PRAT))

# Get the Pheno
Pheno_PRAT <- pData(PRAT_Dataset)

## Keep only the melanoma samples
Pheno_PRAT <- Pheno_PRAT[Pheno_PRAT$source_name_ch1 == "MELANOMA", ]

# The same in expression
Expr_PRAT <- Expr_PRAT[, colnames(Expr_PRAT) %in% rownames(Pheno_PRAT)]

# Check consistency
all(rownames(Pheno_PRAT) == colnames(Expr_PRAT))

# Get the response variable
Pheno_PRAT$Response <- as.factor(Pheno_PRAT$characteristics_ch1.11)
table(Pheno_PRAT$Response)
levels(Pheno_PRAT$Response) <- c("R", "NR", "R", "R")
Pheno_PRAT$Response <- ordered(Pheno_PRAT$Response, levels = c("NR", "R"))
# Get the class of interest
Group_PRAT <- Pheno_PRAT$Response  
names(Group_PRAT) <- rownames(Pheno_PRAT)  

# convert to Z-score
#Expr_PRAT <- t(scale(t(Expr_PRAT), center = T, scale = T))

##################################################################
##################################################################
# Get the RIAZ dataset

#RIAZ_Dataset <- getGEO("GSE91061", GSEMatrix = TRUE, AnnotGPL = TRUE)
#RIAZ_Dataset <- RIAZ_Dataset$GSE91061_series_matrix.txt.gz
#save(RIAZ_Dataset, file = "./Data/RIAZ_Dataset.rda")
load("./Data/RIAZ_Dataset.rda")

# Read the expression file
Expr_RIAZ <- read.delim("./Data/RIAZ_Expr.csv", header = TRUE, row.names = 1, sep = ",")

# Mapping Entrez ID to gene symbols
tmp <- as.character(rownames(Expr_RIAZ))
Expr_RIAZ$GeneSymbol<- mapIds(org.Hs.eg.db,
                              keys=tmp,
                              column="SYMBOL",
                              keytype="ENTREZID",
                              multiVals="first")

Expr_RIAZ <- Expr_RIAZ[!duplicated(Expr_RIAZ$GeneSymbol), ] 
Expr_RIAZ <- Expr_RIAZ[!is.na(Expr_RIAZ$GeneSymbol), ] 
rownames(Expr_RIAZ) <- Expr_RIAZ$GeneSymbol
Expr_RIAZ$GeneSymbol <- NULL

H19 <- grep("^H19", rownames(Expr_RIAZ))
Expr_RIAZ <- Expr_RIAZ[-H19, ]
MiR <- grep("^MIR", rownames(Expr_RIAZ))
Expr_RIAZ <- Expr_RIAZ[-MiR, ]

rownames(Expr_RIAZ) <- gsub("-","", rownames(Expr_RIAZ))
rownames(Expr_RIAZ) <- gsub("_","",rownames(Expr_RIAZ))
rownames(Expr_RIAZ) <- gsub("\\.","",rownames(Expr_RIAZ))

keep <- rowSums(Expr_RIAZ) > 0
table(keep)

Expr_RIAZ <- Expr_RIAZ[keep,]
dim(Expr_RIAZ)

# Log2_transform
Expr_RIAZ <- log2(Expr_RIAZ + 1)

## Get the phenotype
Pheno_RIAZ <- pData(RIAZ_Dataset)

# Change the rownames to be consistent with that in the expression matrix
rownames(Pheno_RIAZ) <- Pheno_RIAZ$title

# Check if sample names are identical
all(rownames(Pheno_RIAZ) == colnames(Expr_RIAZ))

# Modify colnames of the expression
colnames(Expr_RIAZ) <- gsub("\\.", "-", colnames(Expr_RIAZ))

# reorder the column names of expr and rownames of pheno
Expr_RIAZ <- Expr_RIAZ[, order(colnames(Expr_RIAZ))]
Pheno_RIAZ <- Pheno_RIAZ[order(rownames(Pheno_RIAZ)), ]

# Now re-check
all(rownames(Pheno_RIAZ) == colnames(Expr_RIAZ))

# Get the response variable
Pheno_RIAZ$Response <- as.factor(Pheno_RIAZ$`response:ch1`)
table(Pheno_RIAZ$Response)

#Remove unknown TTT
Pheno_RIAZ <- Pheno_RIAZ[!(Pheno_RIAZ$Response == "UNK"), ]

# The same for the Expression
Expr_RIAZ <- Expr_RIAZ[, colnames(Expr_RIAZ) %in% rownames(Pheno_RIAZ)]

# Re-name the levels 
levels(Pheno_RIAZ$Response) <- c("NR", "R", "NR", NA) 

# convert to Z-score
#Expr_RIAZ <- t(scale(t(Expr_RIAZ), center = T, scale = T))

#########
### Divide into 2 datasets (Pre_aPD1 and On_aPD1)

## Pre_aPD1
Pheno_RIAZ_Pre_aPD1 <- Pheno_RIAZ[Pheno_RIAZ$`visit (pre or on treatment):ch1` == "Pre", ]
Expr_RIAZ_Pre_aPD1 <- Expr_RIAZ[,colnames(Expr_RIAZ) %in% rownames(Pheno_RIAZ_Pre_aPD1)]
all(rownames(Pheno_RIAZ_Pre_aPD1) == colnames(Expr_RIAZ_Pre_aPD1))
# Get the class of interest
Group_RIAZ_Pre_aPD1 <- Pheno_RIAZ_Pre_aPD1$Response
table(Group_RIAZ_Pre_aPD1)
names(Group_RIAZ_Pre_aPD1) <- rownames(Pheno_RIAZ_Pre_aPD1)

## On_aPD1
Pheno_RIAZ_On_aPD1 <- Pheno_RIAZ[Pheno_RIAZ$`visit (pre or on treatment):ch1` == "On", ]
Expr_RIAZ_On_aPD1 <- Expr_RIAZ[,colnames(Expr_RIAZ) %in% rownames(Pheno_RIAZ_On_aPD1)]
all(rownames(Pheno_RIAZ_On_aPD1) == colnames(Expr_RIAZ_On_aPD1))
# Get the class of interest
Group_RIAZ_On_aPD1 <- Pheno_RIAZ_On_aPD1$Response
table(Group_RIAZ_On_aPD1)
names(Group_RIAZ_On_aPD1) <- rownames(Pheno_RIAZ_On_aPD1)

################################################################
#################################################################

## Get the MGH dataset

# MGH_Dataset <- getGEO("GSE115821", GSEMatrix = TRUE, AnnotGPL = TRUE)
# MGH_Dataset1 <- MGH_Dataset$`GSE115821-GPL11154_series_matrix.txt.gz`
# MGH_Dataset2 <- MGH_Dataset$`GSE115821-GPL18573_series_matrix.txt.gz`
# save(MGH_Dataset1, MGH_Dataset2, file = "./Data/MGH_Dataset.rda")

load("./Data/MGH_Dataset.rda")

## Get the 2 phenotypes and merge them in one table
tmpPheno1 <- pData(MGH_Dataset1)
tmpPheno2 <- pData(MGH_Dataset2)

Pheno_MGH <- rbind(tmpPheno1, tmpPheno2)

########

## Read the expression file
Expr_MGH <- read.delim("./Data/MGH_Expr.csv", header = TRUE, sep = ",")

# Process and normalize
Expr_MGH <- Expr_MGH[!duplicated(Expr_MGH$Geneid), ]
rownames(Expr_MGH) <- Expr_MGH$Geneid
Expr_MGH <- Expr_MGH[, -c(1:6)]

miR <- grep("^MIR", rownames(Expr_MGH))
Expr_MGH <- Expr_MGH[-miR, ]

rownames(Expr_MGH) <- gsub("-","", rownames(Expr_MGH))
rownames(Expr_MGH) <- gsub("_","",rownames(Expr_MGH))
rownames(Expr_MGH) <- gsub("\\.","",rownames(Expr_MGH))

## Filter expr3 by keeping only the genes with cpm > 1 (raw counts > 15) in at least 50% of samples
#mycpm <- cpm(Expr_MGH)
#plot(expr3_raw[,1],mycpm[,1],xlim=c(0,40),ylim=c(0,3))
#abline(v=15,col=2)
#abline(h=1,col=4)

#thresh <- mycpm > 1
#keep <- rowSums(thresh) >= 18
#table(keep)

#Expr_MGH <- Expr_MGH[keep,]
#dim(Expr_MGH)

Expr_MGH <- DGEList(Expr_MGH)

Expr_MGH <- calcNormFactors(Expr_MGH, method = c("TMM"))
#expr3$samples

plotMD(Expr_MGH,column=2)
abline(h=0,col="grey")

Expr_MGH <- cpm(Expr_MGH, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
boxplot(Expr_MGH)

#########3
## MAke the sample names consistent 
colnames(Expr_MGH) <- gsub("X", "", colnames(Expr_MGH))
colnames(Expr_MGH) <- gsub(".bam", "", colnames(Expr_MGH))
colnames(Expr_MGH) <- gsub("\\.", "-", colnames(Expr_MGH))

rownames(Pheno_MGH) <- Pheno_MGH$title
rownames(Pheno_MGH) <- gsub(".bam", "", rownames(Pheno_MGH))

## Reorder based on rownames and colnames
Expr_MGH <- Expr_MGH[,order(colnames(Expr_MGH))]
Pheno_MGH <- Pheno_MGH[order(rownames(Pheno_MGH)), ]

# Check for consistency
all(rownames(Pheno_MGH) == colnames(Expr_MGH))

# convert to Z-score
#Expr_MGH <- t(scale(t(Expr_MGH), center = T, scale = T))

#########
# Process the phenotype
Pheno_MGH$Response <- as.factor(Pheno_MGH$`response:ch1`)
table(Pheno_MGH$Response)
levels(Pheno_MGH$Response) <- c("NR", "R")

####################
### Divide into 2 datasets (On_aPD1 & On_aCTLA4)

## MGH_on_aPD1
Pheno_MGH_On_aPD1 <- Pheno_MGH[Pheno_MGH$`antibody:ch1` == "anti-PD-1" | Pheno_MGH$`antibody:ch1` == "anti-PD-1+anti-CTLA-4",  ]
Expr_MGH_On_aPD1 <- Expr_MGH[, colnames(Expr_MGH) %in% rownames(Pheno_MGH_On_aPD1)]
all(rownames(Pheno_MGH_On_aPD1) == colnames(Expr_MGH_On_aPD1))
# Get the class of interest
Group_MGH_On_aPD1 <- Pheno_MGH_On_aPD1$Response
table(Group_MGH_On_aPD1)
names(Group_MGH_On_aPD1) <- rownames(Pheno_MGH_On_aPD1)

## MGH_on_aCTLA4
Pheno_MGH_On_aCTLA4 <- Pheno_MGH[Pheno_MGH$`antibody:ch1` == "anti-CTLA-4" | Pheno_MGH$`antibody:ch1` == "anti-PD-1+anti-CTLA-4",  ]
Expr_MGH_On_aCTLA4 <- Expr_MGH[, colnames(Expr_MGH) %in% rownames(Pheno_MGH_On_aCTLA4)]
all(rownames(Pheno_MGH_On_aCTLA4) == colnames(Expr_MGH_On_aCTLA4))
# Get the class of interest
Group_MGH_On_aCTLA4 <- Pheno_MGH_On_aCTLA4$Response
table(Group_MGH_On_aCTLA4)
names(Group_MGH_On_aCTLA4) <- rownames(Pheno_MGH_On_aCTLA4)

##########################################################################
#########################################################################



## Combine Expr
allExpr <- list(Expr_NP, Expr_HUGO, Expr_TCGA,
                Expr_RIAZ_Pre_aPD1, Expr_RIAZ_On_aPD1, Expr_MGH_On_aPD1,
                Expr_MGH_On_aCTLA4)

# names(allExpr) <- c("Expr_NP", "Expr_Chen_Pre_aCTLA4", "Expr_Chen_Pre_aPD1", 
#                     "Expr_Chen_On_aPD1", "Expr_HUGO", "Expr_TCGA", "Expr_PRAT",
#                     "Expr_RIAZ_Pre_aPD1", "Expr_RIAZ_On_aPD1", "Expr_MGH_On_aPD1",
#                     "Expr_MGH_On_aCTLA4")



### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(allExpr, rownames))


### Filter expression for the required samples
exprsAll <- mapply(x=allExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))


# Get the Expression matrices with common genes
Expr_NP <- Expr_NP[commonGenes, ]
#Expr_Chen_Pre_aCTLA4 <- Expr_Chen_Pre_aCTLA4[commonGenes, ]
#Expr_Chen_Pre_aPD1 <- Expr_Chen_Pre_aPD1[commonGenes, ]
#Expr_Chen_On_aPD1 <- Expr_Chen_On_aPD1[commonGenes, ]
Expr_HUGO <- Expr_HUGO[commonGenes, ]
Expr_TCGA <- Expr_TCGA[commonGenes, ]
#Expr_PRAT <- Expr_PRAT[commonGenes, ]
Expr_RIAZ_Pre_aPD1 <- Expr_RIAZ_Pre_aPD1[commonGenes, ]
Expr_RIAZ_On_aPD1 <- Expr_RIAZ_On_aPD1[commonGenes, ]
Expr_MGH_On_aPD1 <- Expr_MGH_On_aPD1[commonGenes, ]
Expr_MGH_On_aCTLA4 <- Expr_MGH_On_aCTLA4[commonGenes, ]

#########3 
##All combined
Melanoma_Mat <- do.call("cbind", exprsAll[-1])

# Normalize between arrays
Melanoma_Mat <- normalizeBetweenArrays(Melanoma_Mat, method = "quantile")


Melanoma_Group <- c(Group_HUGO, Group_TCGA, 
                    Group_RIAZ_Pre_aPD1, Group_RIAZ_On_aPD1,
                    Group_MGH_On_aPD1, Group_MGH_On_aCTLA4)

all(names(Melanoma_Group) == colnames(Melanoma_Mat))
table(Melanoma_Group)
Melanoma_Group <- as.factor(Melanoma_Group)
levels(Melanoma_Group) <- c("NR", "R")
####################################################################
####################################################################

## Save the processed Datasets
save(Expr_NP, Pheno_NP, 
     Expr_Chen_Pre_aCTLA4, Pheno_Chen_Pre_aCTLA4,
     Expr_Chen_Pre_aPD1, Pheno_Chen_Pre_aPD1, 
     Expr_Chen_On_aPD1, Pheno_Chen_On_aPD1, 
     Expr_HUGO, Pheno_HUGO, 
     Expr_TCGA, Pheno_TCGA,
     Expr_PRAT, Pheno_PRAT, 
     Expr_RIAZ_Pre_aPD1, Pheno_RIAZ_Pre_aPD1,
     Expr_RIAZ_On_aPD1, Pheno_RIAZ_On_aPD1,
     Expr_MGH_On_aPD1, Pheno_MGH_On_aPD1,
     Expr_MGH_On_aCTLA4, Pheno_MGH_On_aCTLA4,
     Melanoma_Mat, Melanoma_Group,
     file = "./Objs/ProcessedDatasets.rda")









