###############################################################################
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

load("./Objs/Dataset7.rda")
load("./Objs/Dataset8.rda")

load("./Objs/Dataset9.rda")
load("./Objs/Dataset10.rda")
#
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
Pheno7 <- Dataset7$pheno
pheno8 <- Dataset8$pheno
pheno9 <- Dataset9$pheno
pheno10 <- Dataset10$pheno

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
expr7 <- Dataset7$expr
expr8 <- Dataset8$expr
expr9 <- Dataset9$expr
expr10 <- Dataset10$expr

###################################################################

## Annotate expression

## Expr1
head(rownames(expr1))
# rownames(expr1) <- Dataset1$keys
expr1 <- expr1[!is.na(rownames(expr1)), ]
# 
# rownames(expr1) <- gsub("-", "", rownames(expr1))
# rownames(expr1) <- gsub("_", "", rownames(expr1))
# 
# X1 <- expr1
# ffun <- filterfun(pOverA(p = 0.2, A = 100))
# filt1 <- genefilter(2^X1, ffun)
# expr1 <- expr1[filt1, ]
dim(expr1)

expr1 <- t(scale(t(expr1), center = TRUE, scale = TRUE))

#####################
## expr2
head(rownames(expr2))
#rownames(expr2) <- Dataset2$keys
expr2 <- expr2[!is.na(rownames(expr2)), ]
dim(expr2)

#rownames(expr2) <- gsub("-", "", rownames(expr2))
#rownames(expr2) <- gsub("_", "", rownames(expr2))

# X2 <- expr2
# filt2 <- genefilter(2^X2, ffun)
# expr2 <- expr2[filt2, ]
# dim(expr2)

expr2 <- t(scale(t(expr2), center = TRUE, scale = TRUE))


#####################
## expr3
#rownames(expr3) <- Dataset3$keys
#expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)

#rownames(expr3) <- gsub("-", "", rownames(expr3))
#rownames(expr3) <- gsub("_", "", rownames(expr3))

# X3 <- expr3
# filt3 <- genefilter(2^X3, ffun)
# expr3 <- expr3[filt3, ]
# dim(expr3)

expr3 <- t(scale(t(expr3), center = TRUE, scale = TRUE))

# #######################
dim(expr4)

#rownames(expr4) <- gsub("-", "", rownames(expr4))
#rownames(expr4) <- gsub("_", "", rownames(expr4))

# X4 <- expr4
# filt4 <- genefilter(2^X4, ffun)
# expr4 <- expr4[filt4, ]
# dim(expr4)
# 
expr4 <- t(scale(t(expr4), center = TRUE, scale = TRUE))

######################
## expr5
head(rownames(expr5))

#rownames(expr5) <- gsub("-", "", rownames(expr5))
#rownames(expr5) <- gsub("_", "", rownames(expr5))

rownames(expr5) <- Dataset5$keys
expr5 <- expr5[!is.na(rownames(expr5)), ]
dim(expr5)

# X5 <- expr5
# filt5 <- genefilter(2^X5, ffun)
# expr5 <- expr5[filt5, ]
# dim(expr5)

expr5 <- t(scale(t(expr5), center = TRUE, scale = TRUE))

#####################
## expr6
# processed

#rownames(expr6) <- gsub("-", "", rownames(expr6))
#rownames(expr6) <- gsub("_", "", rownames(expr6))

# X6 <- expr6
# filt6 <- genefilter(2^X6, ffun)
# expr6 <- expr6[filt6, ]
# dim(expr6)
# 

expr6 <- t(scale(t(expr6), center = TRUE, scale = TRUE))

#####################
## Expr7
# Processed

#rownames(expr7) <- gsub("-", "", rownames(expr7))
#rownames(expr7) <- gsub("_", "", rownames(expr7))

# X7 <- expr7
# filt7 <- genefilter(2^X7, ffun)
# expr7 <- expr7[filt7, ]
# dim(expr7)


expr7 <- t(scale(t(expr7), center = TRUE, scale = TRUE))

#######################
## Expr8 
# Processed

expr8 <- t(scale(t(expr8), center = T, scale = T))

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

expr9 <- t(scale(t(expr9), center = T, scale = T))

########################
## Expr10 (Processed)
# rownames(expr10) <- Dataset10$keys
# summary(is.na(rownames(expr10)))
# 
# expr10 <- expr10[!is.na(rownames(expr10)), ]
dim(expr10)
sel = which(apply(expr10, 1, function(x) all(is.finite(x)) ))
expr10 <- expr10[sel,]

expr10 <- t(scale(t(expr10), center = T, scale = T))

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

########
# Pheno7
# Processed
Pheno7$Metastasis <- as.factor(Pheno7$Metastasis) 

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

#########################################################################

## Combine expression and phenotype

AllPheno <- list(pheno1, pheno2, pheno3, pheno4, Pheno5, Pheno6, Pheno7, pheno8, pheno9, pheno10)
names(AllPheno) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")

AllExpr <- list(expr1, expr2, expr3, expr4, expr5, expr6, expr7, expr8, expr9, expr10)
names(AllExpr) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")

GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis, Pheno5$Metastasis, Pheno6$Metastasis, Pheno7$Metastasis, pheno8$Metastasis, pheno9$Metastasis, pheno10$Metastasis)
names(GroupMets) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")
GroupMets <- unlist(GroupMets)
table(GroupMets)

### Find commom subset of genes
#AllExpr2 <- list(expr1, expr2, expr3, expr4, expr5, expr7)

commonGenes <- Reduce("intersect", lapply(AllExpr, rownames))

### Filter expression for the common genes
exprsMetastasis <- mapply(x=AllExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

#save(exprsMetastasis, file = "./Objs/Correlation/exprsMetastasis.rda")


##########
GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis, Pheno5$Metastasis, Pheno6$Metastasis, Pheno7$Metastasis, pheno8$Metastasis, pheno9$Metastasis, pheno10$Metastasis)
names(GroupMets) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")


###################################################################
#########################3
### Leave one study out

# GSE116918 Out
Train <- c("GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")
Test <- c("GSE116918")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE116918
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE116918Out.rda")

############################
# GSE55935 Out
Train <- c("GSE116918", "GSE51066","GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")
Test <- c("GSE55935")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE55935
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE55935Out.rda")

############################
# GSE51066 Out
Train <- c("GSE116918", "GSE55935","GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")
Test <- c("GSE51066")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE51066
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE51066Out.rda")

############################
# GSE46691 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")
Test <- c("GSE46691")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE46691
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE46691Out.rda")

############################
# GSE41408 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")
Test <- c("GSE41408")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE41408
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE41408Out.rda")

############################
# GSE16560 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE41408", "JHU", "HPFS", "GSE25136", "GSE70769")
Test <- c("GSE16560")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE16560
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE16560Out.rda")


############################
# JHU Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE41408", "GSE16560", "HPFS", "GSE25136", "GSE70769")
Test <- c("JHU")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$JHU
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_JHUOut.rda")


############################
# HPFS Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "GSE25136", "GSE70769")
Test <- c("HPFS")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$HPFS
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_HPFSOut.rda")


############################
# GSE25136 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE70769")
Test <- c("GSE25136")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE25136
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE25136Out.rda")

############################
# GSE70769 Out
Train <- c("GSE116918", "GSE55935","GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136")
Test <- c("GSE70769")

## Training
trainMat <- do.call("cbind", exprsMetastasis[Train])
trainGroup <- factor(do.call("c", GroupMets[Train]))
levels(trainGroup) <- c("Mets", "No_Mets")
table(trainGroup)

## Testing
testMat <- exprsMetastasis$GSE70769
testGroup <- factor(do.call("c", GroupMets[Test]))
levels(testGroup) <- c("Mets", "No_Mets")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/LOO/MetastasisData_GSE70769Out.rda")


###########################################################################
###########################################################################
### Stratified Sampling

## Assemble in one data frame

GroupMets <- list(pheno1$Metastasis, pheno2$Metastasis, pheno3$Metastasis, pheno4$Metastasis, Pheno5$Metastasis, Pheno6$Metastasis, Pheno7$Metastasis, pheno8$Metastasis, pheno9$Metastasis, pheno10$Metastasis)
names(GroupMets) <- c("GSE116918", "GSE55935", "GSE51066", "GSE46691", "GSE41408", "GSE16560", "JHU", "HPFS", "GSE25136", "GSE70769")
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
AllPheno$GSE41408$STAGE <- as.factor(AllPheno$GSE41408$`pt-stage:ch1`)
AllPheno$JHU$STAGE <- as.factor(AllPheno$JHU$cstage)
AllPheno$HPFS$STAGE <- as.factor(AllPheno$HPFS$cTNM)
AllPheno$GSE70769$STAGE <- as.factor(AllPheno$GSE70769$`clinical stage:ch1`)


### Covariates of relevance select complete cases: STAGE
AllStage <- lapply(AllPheno, function(x) {
  i <- grep("STAGE", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

levels(AllStage$GSE116918) <- c("T1", "T2", "T3", "T4","TX")
levels(AllStage$GSE55935) <- c("T1", "T2", "T3", "TX")
levels(AllStage$GSE41408) <- c("T2", "T2", "T2", "T2", "T3", "T3", "T4")
levels(AllStage$JHU) <- c("T1", "T1", "T2", "T2", "T2", "T3", "T3", "TX")
levels(AllStage$HPFS) <- c("T1", "T2", "T3")
levels(AllStage$GSE70769) <- c("T1", "T2", "T2", "T3", "TX")


AllStage <- factor(unlist(AllStage))
levels(AllStage)[levels(AllStage) == "TX"] <- NA
levels(AllStage)[levels(AllStage) == ""] <- NA
AllStage[AllStage == ""] <- NA
table(AllStage)

######################################################################
## Age

AllPheno$GSE116918$Age <- as.character(AllPheno$GSE116918$`patient age (years):ch1`)
AllPheno$GSE16560$Age <- as.character(AllPheno$GSE16560$`age:ch1`)
AllPheno$JHU$Age <- as.character(AllPheno$JHU$age.x)
AllPheno$HPFS$Age <- as.character(AllPheno$HPFS$agedx)


### Covariates of relevance select complete cases: AGE
allAGE <- lapply(AllPheno, function(x) {
  i <- grep("^Age", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})

allAGE <- unlist(allAGE)
summary(allAGE)

######################################################################
## Gleason
AllPheno$GSE116918$Gleason <- as.character(AllPheno$GSE116918$`gleason grade:ch1`)
AllPheno$GSE46691$Gleason <- as.character(AllPheno$GSE46691$`gleason score:ch1`)
AllPheno$GSE41408$Gleason <- c(6,6,8,6,6,7,8,8,6,8,6,7,11,6,8,8,6,7,7,6,8,7,7,6,7,6,6,7,7,6,7,7,6,7,6,7,6,6,6,7,7,7,6,6,6,6,8,6)
AllPheno$GSE16560$Gleason <- as.character(AllPheno$GSE16560$`gleason:ch1`)
AllPheno$JHU$Gleason <- as.character(AllPheno$JHU$clings)
AllPheno$HPFS$Gleason <- as.character(AllPheno$HPFS$gleason_sum)

AllPheno$GSE70769$characteristics_ch1 <- gsub("\\=.+", "", AllPheno$GSE70769$characteristics_ch1)
AllPheno$GSE70769$characteristics_ch1 <- gsub("tumour gleason: ", "", AllPheno$GSE70769$characteristics_ch1)
AllPheno$GSE70769$Gleason <- as.character(AllPheno$GSE70769$characteristics_ch1)


AllGleason <- lapply(AllPheno, function(x) {
  i <- grep("^Gleason", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

AllGleason <- unlist(AllGleason)

AllGleason[AllGleason == ""] <- NA
AllGleason[AllGleason == "unknown"] <- NA

table(AllGleason)

#########################################################
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
  pik=inclusionprobabilities(1:nrow(covs), nrow(covs)*0.3),
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
## Second validation
#usedTestMat2 <- expr6[commonGenes, ]
#usedTestGroup2 <- Pheno6$Metastasis
#levels(usedTestGroup2)
#table(usedTestGroup2)


## Save
save(mixTestMat, mixTrainMat, mixTrainGroup, mixTestGroup, file = "./Objs/MetastasisDataGood.rda")
