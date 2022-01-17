##############################################################################
### Mohamed Omar
### 01/07/2019
### Goal : Collecting Prostate Cancer Metastasis data 
#################################################################################

rm(list = ls())

setwd("/Volumes/Macintosh/Research/Projects/Prostate")

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

load("./Data/ProstateData.rda")

## Load Testing data set
load("./Data/Dataset1.rda")
load("./Data/Dataset6.rda")

#Dataset4 <- getGEO("GSE46691", GSEMatrix = TRUE)
#Dataset4 <- Dataset4$GSE46691_series_matrix.txt.gz

## Get the phenotype
#pheno1 <- pData(Dataset1)
pheno1 <- Dataset1$pheno
pheno2 <- Dataset2$pheno
pheno3 <- Dataset3$pheno
pheno4 <- Dataset4$pheno
Pheno5 <- Dataset5$pheno
Pheno6 <- Dataset6$pheno

# pheno5 <- read.table("./Data/prad_eururol_2017/data_clinical_patient.txt", header = TRUE, stringsAsFactors = FALSE)
# rownames(pheno5) <- pheno5$PATIENT_ID
# pheno5 <- pheno5[,-1]
## Get feature data
#FeatData1 <- fData(Dataset1)
#FeatData2 <- fData(Dataset2)

## Get expression
expr1 <- Dataset1$expr
expr2 <- Dataset2$expr
expr3 <- Dataset3$expr
expr4 <- Dataset4$expr
expr5 <- Dataset5$expr
expr6 <- Dataset6$expr

###################################################################

## Annotate expression

## Expr1
head(rownames(expr1))
#rownames(expr1) <- Dataset1$keys
expr1 <- expr1[!is.na(rownames(expr1)), ]
range(expr1)

expr1 <- t(scale(t(expr1), center = TRUE, scale = TRUE))

#####################
## expr2
head(rownames(expr2))
#rownames(expr2) <- Dataset2$keys
expr2 <- expr2[!is.na(rownames(expr2)), ]
dim(expr2)
range(expr2)

# X2 <- expr2
# ffun <- filterfun(pOverA(p = 0.8, A = 100))
# filt2 <- genefilter(2^X2, ffun)
# expr2 <- expr2[filt2, ]

expr2 <- t(scale(t(expr2), center = TRUE, scale = TRUE))


#####################
## expr3
#rownames(expr3) <- Dataset3$keys
#expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)
range(expr3)

# X3 <- expr3
# filt3 <- genefilter(2^X3, ffun)
# expr3 <- expr3[filt3, ]

expr3 <- t(scale(t(expr3), center = TRUE, scale = TRUE))

# #######################
dim(expr4)
range(expr4)

# X4 <- expr4
# filt4 <- genefilter(2^X4, ffun)
# expr4 <- expr4[filt4, ]

expr4 <- t(scale(t(expr4), center = TRUE, scale = TRUE))

######################
## expr5
head(rownames(expr5))

rownames(expr5) <- Dataset5$keys
expr5 <- expr5[!is.na(rownames(expr5)), ]
dim(expr5)
range(expr5)

expr5 <- t(scale(t(expr5), center = TRUE, scale = TRUE))

#####################
## expr6
# processed
range(expr6)

expr6 <- t(scale(t(expr6), center = TRUE, scale = TRUE))

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
# ######################################################################

## Modify pheno1
pheno1$Metastasis <- pheno1$`met event (1=yes, 0=no):ch1`
pheno1$Metastasis[pheno1$Metastasis == 0] <- "No_Mets"
pheno1$Metastasis[pheno1$Metastasis == 1] <- "Mets"
pheno1$Metastasis <- as.factor(pheno1$Metastasis)
table(pheno1$Metastasis)
all(rownames(pheno1) == colnames(expr1))

## Modify pheno2
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

#########################################################################

## Combine expression and phenotype

AllPheno <- list(pheno1, pheno2, pheno3, pheno4, Pheno5, Pheno6)
names(AllPheno) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560")

AllExpr <- list(expr1, expr2, expr3, expr4, expr5, expr6)
names(AllExpr) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560")

GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis, Pheno5$Metastasis, Pheno6$Metastasis)
names(GroupMets) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560")
GroupMets <- unlist(GroupMets)
table(GroupMets)

### Find commom subset of genes
AllExpr2 <- list(expr1, expr2, expr3, expr4, expr5, expr6)

commonGenes <- Reduce("intersect", lapply(AllExpr2, rownames))

### Filter expression for the common genes
exprsMetastasis <- mapply(x=AllExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

#save(exprsMetastasis, file = "./Objs/Correlation/exprsMetastasis.rda")

##########
GroupMets2 <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis, Pheno5$Metastasis, Pheno6$Metastasis)
names(GroupMets2) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560")

###################################################################
#########################3
### Leave one study out

# GSE116918 Out
Train <- c("GSE55935", "GSE51066","GSE46691", "GSE41408", "GSE16560")
Test <- c("GSE116918")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets2[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE116918
testGroup <- factor(do.call("c", GroupMets2[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE116918Out.rda")

############################
# GSE55935 Out
Train <- c("GSE116918", "GSE51066","GSE46691", "GSE41408", "GSE16560")
Test <- c("GSE55935")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets2[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE116918
testGroup <- factor(do.call("c", GroupMets2[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE55935Out.rda")

############################
# GSE51066 Out
Train <- c("GSE116918", "GSE55935","GSE46691", "GSE41408", "GSE16560")
Test <- c("GSE51066")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets2[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE116918
testGroup <- factor(do.call("c", GroupMets2[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE51066Out.rda")

############################
# GSE46691 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE41408", "GSE16560")
Test <- c("GSE46691")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets2[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE116918
testGroup <- factor(do.call("c", GroupMets2[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE46691Out.rda")

############################
# GSE41408 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE16560")
Test <- c("GSE41408")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets2[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE116918
testGroup <- factor(do.call("c", GroupMets2[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE41408Out.rda")

############################
# GSE16560 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE41408")
Test <- c("GSE16560")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets2[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE116918
testGroup <- factor(do.call("c", GroupMets2[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE16560Out.rda")



############################################################################
##########################3
### Stratified Sampling

## Assemble in one data frame

GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis, Pheno5$Metastasis, P
                  
                  
                  