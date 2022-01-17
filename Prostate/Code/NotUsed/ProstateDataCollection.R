###############################################################################
### Mohamed Omar
### 01/07/2019
### Goal : Collecting Prostate Cancer Metastasis data 
#################################################################################

rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

############################################################################
### Load library
require(GEOquery)
require(Biobase)
require(limma)
require(caret)
require(reshape)
require(sampling)
require(genefilter)
library(MetaIntegrator)
library(data.table)

###########################################################################
## Download data sets

#Dataset1 <- getGEO("GSE116918", GSEMatrix = TRUE, AnnotGPL = TRUE)
#Dataset1 <- Dataset1$GSE116918_series_matrix.txt.gz

#Dataset1 <- getGEOData("GSE116918")
#Dataset1 <- Dataset1$originalData$GSE116918

#Dataset2 <- getGEO("GSE55935", GSEMatrix = TRUE, AnnotGPL = TRUE)
#Dataset2 <- Dataset2$GSE55935_series_matrix.txt.gz

#Dataset3 <- getGEOData("GSE51066")
#Dataset3 <- Dataset3$originalData$GSE51066

#save(Dataset1, Dataset2, Dataset3, file = "./Data/Datasets.rda")

load("./Data/Datasets.rda")

#Dataset4 <- getGEO("GSE46691", GSEMatrix = TRUE)
#Dataset4 <- Dataset4$GSE46691_series_matrix.txt.gz

## Get the phenotype
#pheno1 <- pData(Dataset1)
pheno1 <- Dataset1$pheno
pheno2 <- pData(Dataset2)
pheno3 <- Dataset3$pheno
pheno4 <- pData(Dataset4)

# pheno5 <- read.table("./Data/prad_eururol_2017/data_clinical_patient.txt", header = TRUE, stringsAsFactors = FALSE)
# rownames(pheno5) <- pheno5$PATIENT_ID
# pheno5 <- pheno5[,-1]
## Get feature data
#FeatData1 <- fData(Dataset1)
FeatData2 <- fData(Dataset2)

## Get expression
#expr1 <- exprs(Dataset1)
expr1 <- Dataset1$expr
expr2 <- exprs(Dataset2)
expr3 <- Dataset3$expr

#expr4 <- fread("./Data/GSE46691_quantile_normalized.txt.gz", header=TRUE, stringsAsFactors=FALSE, verbose = TRUE)

# Load annotated expr4 
load("/Users/mohamedomar/Documents/Research/Projects/Prostate/Data/expr4.rda")

# expr5 <- read.table("./Data/prad_eururol_2017/data_mRNA_seq_fpkm.txt", header = TRUE, stringsAsFactors = FALSE)
# expr5 <- as.matrix(expr5)
# ## Check if expression is normalized and log-scaled
#boxplot(expr1[,1:50], outline = FALSE)
#boxplot(expr2[,1:50], outline = FALSE) #not log-scaled
#boxplot(expr3[,1:50], outline = FALSE)
######################################

## Annotate expr1
rownames(expr1) <- Dataset1$keys
#summary(is.na(rownames(expr1)))
#rownames(expr1) <- gsub(".+///", "", rownames(expr1))
rownames(expr1) <- gsub("-", "", rownames(expr1))

#expr1 <- aggregate(expr1[,], list(Gene = rownames(expr1)), FUN = mean)
#rownames(expr1) <- expr1$Gene

#expr1$Gene <- NULL
expr1 <- expr1[!(rownames(expr1) == ""), ]
expr1 <- expr1[!is.na(rownames(expr1)), ]

dim(expr1)

# X1 <- expr1
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt1 <- genefilter(2^X1,ffun)
# expr1 <- expr1[filt1,]

expr1 <- t(scale(t(expr1), scale = T, center = T))

##############################################

## Annotate expr2
rownames(expr2) <- FeatData2$`Gene symbol`
#summary(is.na(rownames(expr2)))
expr2 <- expr2[!(rownames(expr2) == ""), ]
rownames(expr2) <- gsub("-", "", rownames(expr2))

dim(expr2)
## Take the log2 of expr2
expr2 <- log2(expr2)

# X2 <- expr2
# filt2 <- genefilter(2^X2,ffun)
# expr2 <- expr2[filt2,]

expr2 <- t(scale(t(expr2), scale = T, center = T))
##################################################

rownames(expr3) <- Dataset3$keys
#summary(is.na(rownames(expr3)))
expr3 <- expr3[!(rownames(expr3) == ""), ]
expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)

#expr3 <- aggregate(expr3[,], list(Gene = rownames(expr3)), FUN = median)
#rownames(expr3) <- expr3$Gene

#expr3$Gene <- NULL

# dim(expr3)
# X3 <- expr3
# filt3 <- genefilter(2^X3, ffun)
# expr3 <- expr3[filt3, ]

rownames(expr3) <- gsub("-", "", rownames(expr3))

expr3 <- t(scale(t(expr3), scale = T, center = T))

####################################################
## Annotate expr4
#rownames(expr4) <- expr4$ID_REF
#expr4$ID_REF <- NULL
#expr4 <- as.matrix(expr4)
#rownames(expr4) <- Dataset3$keys
#summary(is.na(rownames(expr4)))
#expr4 <- expr4[!is.na(rownames(expr4)), ]
#dim(expr4)

#boxplot(expr4[,1:25])

## Save annotated expr4 for further faster use
#save(expr4, file = "./Data/expr4.rda")

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

#expr4 <- expr4[ , match(pheno4$description, colnames(expr4))]
#colnames(expr4) <- rownames(pheno4)
# X4 <- expr4
# filt4 <- genefilter(2^X4, ffun)
# expr4 <- expr4[filt4, ]

dim(expr4)
expr4 <- t(scale(t(expr4), scale = T, center = T))


##### 
# rownames(expr5) <- expr5[,"Hugo_Symbol"]
# expr5 <- expr5[,-c(1:2)]
# expr5 <- expr5[!is.na(rownames(expr5)), ]
# expr5 <- expr5[!(rownames(expr5) == ""), ]
# class(expr5) <- "numeric"
# 
# dim(expr5)
# X5 <- expr5
# ffun <- filterfun(pOverA(p = 1, A = 1))
# filt5 <- genefilter(X5,ffun)
# expr5 <- expr5[filt5, ]
# 
# expr5 <- log2(expr5+0.1)
# 
# head(colnames(expr5))
# expr5 <- expr5[,order(colnames(expr5))]
# expr5 <- t(scale(t(expr5)))
# 
# ######################################################################

## Modify pheno1
pheno1$Metastasis <- pheno1$`met event (1=yes, 0=no):ch1`
pheno1$Metastasis[pheno1$Metastasis == 0] <- "No_Mets"
pheno1$Metastasis[pheno1$Metastasis == 1] <- "Mets"
pheno1$Metastasis <- as.factor(pheno1$Metastasis)
table(pheno1$Metastasis)
all(rownames(pheno1) == colnames(expr1))

## Modify pheno2
# Remove cell lines
pheno2 <- pheno2[-c(47:54), ]
pheno2$Metastasis <- pheno2$`lymph node metastasis status:ch1`
pheno2$Metastasis[pheno2$Metastasis == 0] <- "No_Mets"
pheno2$Metastasis[pheno2$Metastasis == 1] <- "Mets"
pheno2 <- pheno2[!(pheno2$Metastasis == "NA"), ]

pheno2$Metastasis <- as.factor(pheno2$Metastasis)
table(pheno2$Metastasis)
levels(pheno2$Metastasis)
# Modify expr2
expr2 <- expr2[,colnames(expr2) %in% rownames(pheno2)]
all(rownames(pheno2) == colnames(expr2))
####

## Modify pheno3
pheno3$Metastasis <- pheno3$`metastatic event:ch1`
pheno3$Metastasis[pheno3$Metastasis == 0] <- "No_Mets"
pheno3$Metastasis[pheno3$Metastasis == 1] <- "Mets"

pheno3$Metastasis <- as.factor(pheno3$Metastasis)
table(pheno3$Metastasis)
levels(pheno3$Metastasis)

all(rownames(pheno3) == colnames(expr3))
#####

## Modify pheno4
pheno4$Metastasis <- pheno4$`metastatic event:ch1`
pheno4$Metastasis[pheno4$Metastasis == 0] <- "No_Mets"
pheno4$Metastasis[pheno4$Metastasis == 1] <- "Mets"
pheno4 <- pheno4[!(pheno4$Metastasis == "NA"), ]

pheno4$Metastasis <- as.factor(pheno4$Metastasis)
table(pheno4$Metastasis)
levels(pheno4$Metastasis)

all(rownames(pheno4) == colnames(expr4))
#####

## Modify pheno5
# head(rownames(pheno5))
# pheno5 <- pheno5[order(rownames(pheno5)), ]
# colnames(expr5) <- rownames(pheno5)
# 
# pheno5$Metastasis <- pheno5$LYMPH_NODE_METASTASIS
# pheno5$Metastasis[pheno5$Metastasis == "No"] <- "No_Mets"
# pheno5$Metastasis[pheno5$Metastasis == "Yes"] <- "Mets"
# pheno5 <- pheno5[!is.na(pheno5$Metastasis), ]
# 
# pheno5$Metastasis <- as.factor(pheno5$Metastasis)
# table(pheno5$Metastasis)
# 
# expr5 <- expr5[,colnames(expr5) %in% rownames(pheno5)]
# all(rownames(pheno5) == colnames(expr5))
# 
############
#Sub <- subset(pheno4, pheno4$Metastasis == "Mets")

#pheno4_Test <- Sub[1:106, ]  # To be added to the testing set
#expr4_Test <- expr4[, colnames(expr4) %in% rownames(pheno4_Test)]  # To be added to the testing set
#all(rownames(pheno4_Test) == colnames(expr4_Test))

#pheno4_Train <- pheno4[!(rownames(pheno4) %in% rownames(pheno4_Test)), ]
#expr4_Train <- expr4[,!(colnames(expr4) %in% colnames(expr4_Test))]
#all(rownames(pheno4_Train) == colnames(expr4_Train))
#########################################################################

## Combine expression and phenotype

AllPheno <- list(pheno1, pheno2, pheno3, pheno4)
names(AllPheno) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691")

AllExpr <- list(expr1, expr2, expr3, expr4)
names(AllExpr) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691")

GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis)
names(GroupMets) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691")
GroupMets <- unlist(GroupMets)
table(GroupMets)

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(AllExpr, rownames))

### Filter expression for the common genes
exprsMetastasis <- mapply(x=AllExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

##########
GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis)
names(GroupMets) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691")


Train <- c("GSE55935", "GSE51066","GSE46691")
Test <- c("GSE116918")

### Pheno
Pheno_Train <- AllPheno[Train]
Pheno_Test <- AllPheno$GSE116918
# 
## Training Mat and Group
trainMat <- exprsMetastasis[Train]
trainMat <- do.call("cbind", trainMat)
# 
trainGroup <- GroupMets[Train]
trainGroup <- unlist(trainGroup)
table(trainGroup)
# 
names(trainGroup) <- colnames(trainMat)
# 
# ## Testing Mat and Group
testMat <- exprsMetastasis$GSE116918
##mixTestMat <- mixTestMat$GSE116918
# 
testGroup <- GroupMets$GSE116918
# mixTestGroup <- unlist(mixTestGroup)
table(testGroup)
# 
names(testGroup) <- colnames(testMat)


## Assemble in one data frame

GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis)
names(GroupMets) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691")
GroupMets <- unlist(GroupMets)
table(GroupMets)


AllMat <- do.call("cbind", exprsMetastasis)

## Check if sample names are identical
all(colnames(AllMat) == names(GroupMets))
names(GroupMets) <- colnames(AllMat)
#########################################################################
## Stage

AllPheno$GSE116918$STAGE <- as.factor(AllPheno$GSE116918$`t-stage:ch1`)
AllPheno$GSE55935$STAGE <- as.factor(AllPheno$GSE55935$`clinical t-stage (ct):ch1`)

### Covariates of relevance select complete cases: STAGE
AllStage <- lapply(AllPheno, function(x) {
  i <- grep("STAGE", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

levels(AllStage$GSE116918) <- c("T1", "T2", "T3", "T4","TX")
levels(AllStage$GSE55935) <- c("T1", "T2", "T3", "TX")

AllStage <- factor(unlist(AllStage))
levels(AllStage)[levels(AllStage) == "TX"] <- NA
AllStage[AllStage == ""] <- NA
######################################################################
## Age

AllPheno$GSE116918$Age <- as.character(AllPheno$GSE116918$`patient age (years):ch1`)

### Covariates of relevance select complete cases: AGE
allAGE <- lapply(AllPheno, function(x) {
  i <- grep("^Age", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})

allAGE <- unlist(allAGE)

######################################################################
## Gleason
AllPheno$GSE116918$Gleason <- as.character(AllPheno$GSE116918$`gleason grade:ch1`)
AllPheno$GSE46691$Gleason <- as.character(AllPheno$GSE46691$`gleason score:ch1`)

AllGleason <- lapply(AllPheno, function(x) {
  i <- grep("Gleason", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

AllGleason <- unlist(AllGleason)

AllGleason[AllGleason == ""] <- NA

### Assemble in one data.frame and turn numeric
covs <- data.frame(Stage = AllStage, Age = allAGE, Gleason = AllGleason)

### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))) )

###########################################################################
###SAMPLING

### Balanced stratification
set.seed(333)
trainingOrTesting <- balancedstratification(
  covs[ , , drop=FALSE], strata=1*(GroupMets == "Mets"),
  pik=inclusionprobabilities(1:nrow(covs), nrow(covs)*0.2),
  comment=TRUE, method=1)

### Show
apply(covs[, ,drop=FALSE], 2, table, GroupMets, trainingOrTesting)

### Subset Training
mixTrainMat <- AllMat[ , trainingOrTesting == 0]
mixTrainGroup <- GroupMets[ trainingOrTesting == 0]
table(mixTrainGroup)
### Subset Testing
mixTestMat <- AllMat[ , trainingOrTesting == 1]
mixTestGroup <- GroupMets[ trainingOrTesting == 1]
table(mixTestGroup)

##########################################################################
## Divide data into training and testing data
#set.seed(333)
#ind <- createDataPartition(y = GroupMets, p = 0.5, list = FALSE)

#mixTrainGroup <- GroupMets[ind]
#mixTestGroup <- GroupMets[-ind]

#mixTrainMat <- AllMat[, ind]
#mixTestMat <- AllMat[, -ind]

## Save
save(mixTestMat, mixTrainMat, mixTrainGroup, mixTestGroup, trainMat, trainGroup, Pheno_Train, testMat, testGroup, Pheno_Test, file = "./Objs/MetastasisData.rda")
