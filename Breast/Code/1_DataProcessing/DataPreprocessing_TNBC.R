#############################################################################
#############################################################################
## Mohamed Omar
## December 2, 2019
## Goal: Getting and organizing the data for the next step
############################################################################
rm(list = ls())


library(GEOquery)
library(Biobase)
library(sampling)
library(limma)
library(genefilter)
library(edgeR)
library(caret)

### Getting the data 
# dataset1 <- getGEO("GSE25055", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset1 <- dataset1$GSE25055_series_matrix.txt.gz
# 
# dataset2 <- getGEO("GSE25065", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset2 <- dataset2$GSE25065_series_matrix.txt.gz

#save(dataset1, dataset2, file = "./Data/datasets.rda")

load("./Data/datasets_new.rda")
###################################################

## Get the phenotype
Pheno1 <- pData(dataset1)
Pheno2 <- pData(dataset2)
Pheno3 <- pData(dataset3)
Pheno4 <- pData(dataset4)
#Pheno5 <- pData(dataset5)
Pheno6 <- pData(dataset6)
Pheno7 <- pData(dataset7)
Pheno8 <- pData(dataset8)

## Get the expression matrices
Expr1 <- exprs(dataset1)
Expr2 <- exprs(dataset2)
Expr3 <- exprs(dataset3)
Expr4 <- exprs(dataset4)
#Expr5 <- exprs(dataset5)
Expr6 <- exprs(dataset6)
Expr7 <- exprs(dataset7)
Expr8 <- exprs(dataset8)

## Get the feature data
FeatData1 <- fData(dataset1)
FeatData2 <- fData(dataset2)
FeatData3 <- fData(dataset3)
FeatData4 <- fData(dataset4)
#FeatData5 <- fData(dataset5)
FeatData6 <- fData(dataset6)
FeatData7 <- fData(dataset7)
FeatData8 <- fData(dataset8)

## Annotation of Expr1
rownames(Expr1) <- FeatData1$`Gene symbol`
summary(is.na(rownames(Expr1)))
rownames(Expr1) <- gsub("-","", rownames(Expr1))
rownames(Expr1) <- gsub("_","",rownames(Expr1))
sel <- which(apply(Expr1, 1, function(x) all(is.finite(x)) ))
Expr1 <- Expr1[sel, ]
Expr1 <- Expr1[!is.na(rownames(Expr1)),]
Expr1 <- Expr1[!(rownames(Expr1) == ""), ] 
dim(Expr1)

## Filtering
# X1 <- Expr1
# ffun <- pOverA(p = 1, A = 50)
# Filt1 <- genefilter(2^X1, ffun)
# Expr1 <- Expr1[Filt1, ]
# dim(Expr1)

#Expr1 <- t(scale(t(Expr1), center = TRUE, scale = TRUE))

## Annotation of Expr2
rownames(Expr2) <- FeatData2$`Gene symbol`
summary(is.na(rownames(Expr2)))
rownames(Expr2) <- gsub("-","", rownames(Expr2))
rownames(Expr2) <- gsub("_","",rownames(Expr2))
sel <- which(apply(Expr2, 1, function(x) all(is.finite(x)) ))
Expr2 <- Expr2[sel, ]
Expr2 <- Expr2[!is.na(rownames(Expr2)),]
Expr2 <- Expr2[!(rownames(Expr2) == ""), ] 
dim(Expr2)

## Filtering
# X2 <- Expr2
# Filt2 <- genefilter(2^X2, ffun)
# Expr2 <- Expr2[Filt2, ]
# dim(Expr2)
# 
#Expr2 <- t(scale(t(Expr2), center = TRUE, scale = TRUE))

## Annotation of Expr3
rownames(Expr3) <- FeatData3$`Gene symbol`
summary(is.na(rownames(Expr3)))
rownames(Expr3) <- gsub("-","", rownames(Expr3))
rownames(Expr3) <- gsub("_","",rownames(Expr3))
sel <- which(apply(Expr3, 1, function(x) all(is.finite(x)) ))
Expr3 <- Expr3[sel, ]
Expr3 <- Expr3[!is.na(rownames(Expr3)),]
Expr3 <- Expr3[!(rownames(Expr3) == ""), ] 
dim(Expr3)

## Filtering
# X3 <- Expr3
# Filt3 <- genefilter(2^X3, ffun)
# Expr3 <- Expr3[Filt3, ]
# dim(Expr3)

## Annotation of Expr4
rownames(Expr4) <- FeatData4$`Gene symbol`
summary(is.na(rownames(Expr4)))
rownames(Expr4) <- gsub("-","", rownames(Expr4))
rownames(Expr4) <- gsub("_","",rownames(Expr4))
sel <- which(apply(Expr4, 1, function(x) all(is.finite(x)) ))
Expr4 <- Expr4[sel, ]
Expr4 <- Expr4[!is.na(rownames(Expr4)),]
Expr4 <- Expr4[!(rownames(Expr4) == ""), ] 
dim(Expr4)

## Filtering
# X4 <- Expr4
# Filt4 <- genefilter(2^X4, ffun)
# Expr4 <- Expr4[Filt4, ]
# dim(Expr4)

## Annotation of Expr6
rownames(Expr6) <- FeatData6$`Gene symbol`
summary(is.na(rownames(Expr6)))
rownames(Expr6) <- gsub("-","", rownames(Expr6))
rownames(Expr6) <- gsub("_","",rownames(Expr6))
sel <- which(apply(Expr6, 1, function(x) all(is.finite(x)) ))
Expr6 <- Expr6[sel, ]
Expr6 <- Expr6[!is.na(rownames(Expr6)),]
Expr6 <- Expr6[!(rownames(Expr6) == ""), ] 
dim(Expr6)

## Filtering
# X6 <- Expr6
# Filt6 <- genefilter(2^X6, ffun)
# Expr6 <- Expr6[Filt6, ]
# dim(Expr6)

## Annotation of Expr7
rownames(Expr7) <- FeatData7$`Gene symbol`
summary(is.na(rownames(Expr7)))
rownames(Expr7) <- gsub("-","", rownames(Expr7))
rownames(Expr7) <- gsub("_","",rownames(Expr7))
sel <- which(apply(Expr7, 1, function(x) all(is.finite(x)) ))
Expr7 <- Expr7[sel, ]
Expr7 <- Expr7[!is.na(rownames(Expr7)),]
Expr7 <- Expr7[!(rownames(Expr7) == ""), ] 
dim(Expr7)

# Expression values are not logged
range(Expr7)
Expr7 <- log2(Expr7 + 1)

## Filtering
# X7 <- Expr7
# Filt7 <- genefilter(2^X7, ffun)
# Expr7 <- Expr7[Filt7, ]
# dim(Expr7)

## Annotation of Expr8
rownames(Expr8) <- FeatData8$`Gene symbol`
summary(is.na(rownames(Expr8)))
rownames(Expr8) <- gsub("-","", rownames(Expr8))
rownames(Expr8) <- gsub("_","",rownames(Expr8))
sel <- which(apply(Expr8, 1, function(x) all(is.finite(x)) ))
Expr8 <- Expr8[sel, ]
Expr8 <- Expr8[!is.na(rownames(Expr8)),]
Expr8 <- Expr8[!(rownames(Expr8) == ""), ] 
dim(Expr8)

## Filtering
# X8 <- Expr8
# Filt8 <- genefilter(2^X8, ffun)
# Expr8 <- Expr8[Filt8, ]
# dim(Expr8)


###########################################
### Modify the phenotypes

#Keep only triple-negative Breast cancer TNBC
table(Pheno1$`er_status_ihc:ch1`)
Pheno1_TNBC <- Pheno1[Pheno1$`er_status_ihc:ch1` == "N", ]

table(Pheno1_TNBC$`pr_status_ihc:ch1`)
Pheno1_TNBC <- Pheno1_TNBC[Pheno1_TNBC$`pr_status_ihc:ch1` == "N", ]

table(Pheno1_TNBC$`her2_status:ch1`)
Pheno1_TNBC <- Pheno1_TNBC[Pheno1_TNBC$`her2_status:ch1` == "N", ]

# Convert to factor 
table(Pheno1_TNBC$`pathologic_response_pcr_rd:ch1`)
Pheno1_TNBC <- Pheno1_TNBC[!(Pheno1_TNBC$`pathologic_response_pcr_rd:ch1` == "NA"), ]
Pheno1_TNBC$ChemoResponse <- Pheno1_TNBC$`pathologic_response_pcr_rd:ch1`
Pheno1_TNBC$ChemoResponse[Pheno1_TNBC$ChemoResponse == "RD"] <- "Resistant"
Pheno1_TNBC$ChemoResponse[Pheno1_TNBC$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno1_TNBC$ChemoResponse)
Pheno1_TNBC$ChemoResponse <- factor(Pheno1_TNBC$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno1_TNBC$ChemoResponse)

## Clean pheno1
# colnames(Pheno1)
# keep <- c("characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4",                   
#           "characteristics_ch1.5", "characteristics_ch1.7",                    
#           "characteristics_ch1.8", "characteristics_ch1.9", "characteristics_ch1.10", 
#           "characteristics_ch1.13", "characteristics_ch1.14", "characteristics_ch1.15", 
#           "characteristics_ch1.16", "characteristics_ch1.17",                  
#           "characteristics_ch1.19", "characteristics_ch1.20", "ChemoResponse")

# Pheno1 <- Pheno1[, keep]
# colnames(Pheno1) <- c("Age_Years", "ER_Status", "PR_Status", "HER2_Status", "cT_Stage", 
#                       "cN_Stage", "AJCC_Stage", "Grade", "DRFS_Event", "DRFS_Year", "ESR1_Status",
#                       "ERBB2_Status", "SET_Class", "GGI_Class", "PAM50_Class", "ChemoResponse")
# 

Expr1_TNBC <- Expr1[, colnames(Expr1)%in%rownames(Pheno1_TNBC)]
all(rownames(Pheno1_TNBC) == colnames(Expr1_TNBC))

#####################
## Pheno2

# Keep only TNBC
table(Pheno2$`er_status_ihc:ch1`)
Pheno2_TNBC <- Pheno2[Pheno2$`er_status_ihc:ch1` == "N", ]

table(Pheno2_TNBC$`pr_status_ihc:ch1`)
Pheno2_TNBC <- Pheno2_TNBC[Pheno2_TNBC$`pr_status_ihc:ch1` == "N", ]

table(Pheno2_TNBC$`her2_status:ch1`)
Pheno2_TNBC <- Pheno2_TNBC[Pheno2_TNBC$`her2_status:ch1` == "N", ]

table(Pheno2_TNBC$`pathologic_response_pcr_rd:ch1`)
Pheno2_TNBC <- Pheno2_TNBC[!(Pheno2_TNBC$`pathologic_response_pcr_rd:ch1` == "NA"), ]
Pheno2_TNBC$ChemoResponse <- Pheno2_TNBC$`pathologic_response_pcr_rd:ch1`
Pheno2_TNBC$ChemoResponse[Pheno2_TNBC$ChemoResponse == "RD"] <- "Resistant"
Pheno2_TNBC$ChemoResponse[Pheno2_TNBC$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno2_TNBC$ChemoResponse)
Pheno2_TNBC$ChemoResponse <- factor(Pheno2_TNBC$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno2_TNBC$ChemoResponse)

## Clean pheno2
# colnames(Pheno2)
# keep <- c("characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4",                   
#           "characteristics_ch1.5", "characteristics_ch1.6", "characteristics_ch1.7",                    
#           "characteristics_ch1.8", "characteristics_ch1.9", "characteristics_ch1.12",
#           "characteristics_ch1.13", "characteristics_ch1.15", 
#           "characteristics_ch1.16", "characteristics_ch1.17",                  
#           "characteristics_ch1.19", "characteristics_ch1.20", "ChemoResponse")
# 
# Pheno2 <- Pheno2[, keep]
# colnames(Pheno2) <- c("Age_Years", "ER_Status", "PR_Status", "HER2_Status", "cT_Stage", 
#                       "cN_Stage", "AJCC_Stage", "Grade", "DRFS_Event", "DRFS_Year", "ESR1_Status",
#                       "ERBB2_Status", "SET_Class", "GGI_Class",  "PAM50_Class", "ChemoResponse")

Expr2_TNBC <- Expr2[, colnames(Expr2)%in%rownames(Pheno2_TNBC)]
all(rownames(Pheno2_TNBC) == colnames(Expr2_TNBC))

##############
# Modify Pheno3
#Keep only triple-negative Breast cancer TNBC
table(Pheno3$`er status:ch1`)
Pheno3_TNBC <- Pheno3[Pheno3$`er status:ch1` == "0", ]

table(Pheno3_TNBC$`pr status:ch1`)
Pheno3_TNBC <- Pheno3_TNBC[Pheno3_TNBC$`pr status:ch1` == "0", ]

table(Pheno3_TNBC$`her2 status:ch1`)
Pheno3_TNBC <- Pheno3_TNBC[Pheno3_TNBC$`her2 status:ch1` == "0", ]

# Convert to factor 
table(Pheno3_TNBC$`pathological response:ch1`)
Pheno3_TNBC$ChemoResponse <- Pheno3_TNBC$`pathological response:ch1`
Pheno3_TNBC$ChemoResponse[Pheno3_TNBC$ChemoResponse %in% c("pNC", "pPR")] <- "Resistant"
Pheno3_TNBC$ChemoResponse[Pheno3_TNBC$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno3_TNBC$ChemoResponse)
Pheno3_TNBC$ChemoResponse <- factor(Pheno3_TNBC$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno3_TNBC$ChemoResponse)

Expr3_TNBC <- Expr3[, colnames(Expr3) %in% rownames(Pheno3_TNBC)]
all(rownames(Pheno3_TNBC) == colnames(Expr3_TNBC))

################
# Modify pheno4
# Convert to factor 
table(Pheno4$`response:ch1`)
Pheno4$ChemoResponse <- Pheno4$`response:ch1`
Pheno4$ChemoResponse[Pheno4$ChemoResponse == "0"] <- "Resistant"
Pheno4$ChemoResponse[Pheno4$ChemoResponse == "1"] <- "Sensitive"
table(Pheno4$ChemoResponse)
Pheno4$ChemoResponse <- factor(Pheno4$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno4$ChemoResponse)

#Expr4 <- Expr4[, colnames(Expr4) %in% rownames(Pheno4)]
all(rownames(Pheno4) == colnames(Expr4))

###################
# Modify pheno6

# Keep Triple negative
table(Pheno6$`er_status:ch1`)
Pheno6_TNBC <- Pheno6[Pheno6$`er_status:ch1`== "N", ]

table(Pheno6_TNBC$`pr_status:ch1`)
Pheno6_TNBC <- Pheno6_TNBC[Pheno6_TNBC$`pr_status:ch1`== "N", ]

table(Pheno6_TNBC$`her2 status:ch1`)
Pheno6_TNBC <- Pheno6_TNBC[Pheno6_TNBC$`her2 status:ch1`== "N", ]

#Pheno6_TNBC <- Pheno6_TNBC[!(Pheno6_TNBC$`pcr.v.rd:ch1` == "NA"), ]

# Convert to factor 
table(Pheno6_TNBC$`pcr_vs_rd:ch1`)
Pheno6_TNBC$ChemoResponse <- Pheno6_TNBC$`pcr_vs_rd:ch1`
Pheno6_TNBC$ChemoResponse[Pheno6_TNBC$ChemoResponse == "RD" ] <- "Resistant"
Pheno6_TNBC$ChemoResponse[Pheno6_TNBC$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno6_TNBC$ChemoResponse)
Pheno6_TNBC$ChemoResponse <- factor(Pheno6_TNBC$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno6_TNBC$ChemoResponse)

Expr6_TNBC <- Expr6[, colnames(Expr6) %in% rownames(Pheno6_TNBC)]
all(rownames(Pheno6_TNBC) == colnames(Expr6_TNBC))

###############
# Modify pheno7

# Keep Triple negative
table(Pheno7$`er status:ch1`)
Pheno7_TNBC <- Pheno7[Pheno7$`er status:ch1`== "N", ]

table(Pheno7_TNBC$`pr status:ch1`)
Pheno7_TNBC <- Pheno7_TNBC[Pheno7_TNBC$`pr status:ch1`== "N", ]

table(Pheno7_TNBC$`her 2 status:ch1`)
Pheno7_TNBC <- Pheno7_TNBC[Pheno7_TNBC$`her 2 status:ch1`== "N", ]

#Pheno7_TNBC <- Pheno7_TNBC[!(Pheno7_TNBC$`pcr.v.rd:ch1` == "NA"), ]

# Convert to factor 
table(Pheno7_TNBC$`pcr or rd:ch1`)
Pheno7_TNBC$ChemoResponse <- Pheno7_TNBC$`pcr or rd:ch1`
Pheno7_TNBC$ChemoResponse[Pheno7_TNBC$ChemoResponse == "RD" ] <- "Resistant"
Pheno7_TNBC$ChemoResponse[Pheno7_TNBC$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno7_TNBC$ChemoResponse)
Pheno7_TNBC$ChemoResponse <- factor(Pheno7_TNBC$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno7_TNBC$ChemoResponse)

## 
Expr7_TNBC <- Expr7[, colnames(Expr7) %in% rownames(Pheno7_TNBC)]
all(rownames(Pheno7_TNBC) == colnames(Expr7_TNBC))

##########################
# Modify pheno8

# Keep Triple negative
table(Pheno8$`er status ihc:ch1`)
Pheno8_TNBC <- Pheno8[Pheno8$`er status ihc:ch1`== "negative", ]

table(Pheno8_TNBC$`pr status ihc:ch1`)

table(Pheno8_TNBC$`her2 status fish:ch1`)
Pheno8_TNBC <- Pheno8_TNBC[Pheno8_TNBC$`her2 status fish:ch1`== "negative", ]

#Pheno8_TNBC <- Pheno8_TNBC[!(Pheno8_TNBC$`pcr.v.rd:ch1` == "NA"), ]

# Convert to factor 
table(Pheno8_TNBC$`pathologic response pcr ncr:ch1`)
Pheno8_TNBC$ChemoResponse <- Pheno8_TNBC$`pathologic response pcr ncr:ch1`
Pheno8_TNBC$ChemoResponse[Pheno8_TNBC$ChemoResponse == "nCR" ] <- "Resistant"
Pheno8_TNBC$ChemoResponse[Pheno8_TNBC$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno8_TNBC$ChemoResponse)
Pheno8_TNBC$ChemoResponse <- factor(Pheno8_TNBC$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno8_TNBC$ChemoResponse)

## 
Expr8_TNBC <- Expr8[, colnames(Expr8) %in% rownames(Pheno8_TNBC)]
all(rownames(Pheno8_TNBC) == colnames(Expr8_TNBC))


##################################
## Aggregate

allPheno_TNBC <- list(Pheno1_TNBC, Pheno2_TNBC, Pheno3_TNBC, Pheno4, Pheno6_TNBC, Pheno7_TNBC, Pheno8_TNBC)
names(allPheno_TNBC) <- c("GSE25055", "GSE25065", "GSE140494", "GSE103668", "GSE20194", "GSE20271", "GSE32646")

ChemoSensitivity_ALL_TNBC <- list(Pheno1_TNBC$ChemoResponse, Pheno2_TNBC$ChemoResponse, Pheno3_TNBC$ChemoResponse, Pheno4$ChemoResponse, Pheno6_TNBC$ChemoResponse, Pheno7_TNBC$ChemoResponse, Pheno8_TNBC$ChemoResponse)
names(ChemoSensitivity_ALL_TNBC) <- c("GSE25055", "GSE25065", "GSE140494", "GSE103668", "GSE20194", "GSE20271", "GSE32646")

allExpr_TNBC <- list(Expr1_TNBC, Expr2_TNBC, Expr3_TNBC, Expr4, Expr6_TNBC, Expr7_TNBC, Expr8_TNBC)

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(allExpr_TNBC, rownames))

### Filter expression for the required samples
exprsChemo_TNBC <- mapply(x=allExpr_TNBC, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

names(exprsChemo_TNBC) <- c("GSE25055", "GSE25065", "GSE140494", "GSE103668", "GSE20194", "GSE20271", "GSE32646")

### Check
all(names(exprsChemo_TNBC) == names(ChemoSensitivity_ALL_TNBC))

### Check order
all(rownames(allPheno_TNBC$GSE25055) == colnames(exprsChemo_TNBC$GSE25055))
all(rownames(allPheno_TNBC$GSE25065) == colnames(exprsChemo_TNBC$GSE25065))
all(rownames(allPheno_TNBC$GSE140494) == colnames(exprsChemo_TNBC$GSE140494))
all(rownames(allPheno_TNBC$GSE103668) == colnames(exprsChemo_TNBC$GSE103668))
all(rownames(allPheno_TNBC$GSE20194) == colnames(exprsChemo_TNBC$GSE20194))
all(rownames(allPheno_TNBC$GSE20271) == colnames(exprsChemo_TNBC$GSE20271))
all(rownames(allPheno_TNBC$GSE32646) == colnames(exprsChemo_TNBC$GSE32646))

###################################################################
#############
## Cross-study validation

# Leave GSE25055 out
train <- c("GSE25065", "GSE140494", "GSE103668", "GSE20194", "GSE20271", "GSE32646")
test <- c("GSE25055")

## Training
trainMat <- do.call("cbind", exprsChemo_TNBC[train])
trainGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[train]))
levels(trainGroup) <- c("Sensitive", "Resistant")
table(trainGroup)

## Testing
testMat <- exprsChemo_TNBC$GSE25055
testGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[test]))
levels(testGroup) <- c("Sensitive", "Resistant")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/ChemoDataNew_GSE25055Out.rda")

##########################
# Leave GSE25065 out
train <- c("GSE25055", "GSE140494", "GSE103668", "GSE20194", "GSE20271", "GSE32646")
test <- c("GSE25065")

## Training
trainMat <- do.call("cbind", exprsChemo_TNBC[train])
trainGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[train]))
levels(trainGroup) <- c("Sensitive", "Resistant")
table(trainGroup)

## Testing
testMat <- exprsChemo_TNBC$GSE25065
testGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[test]))
levels(testGroup) <- c("Sensitive", "Resistant")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/ChemoDataNew_GSE25065Out.rda")

##########################
# Leave GSE140494 out
train <- c("GSE25055", "GSE25065", "GSE103668", "GSE20194", "GSE20271", "GSE32646")
test <- c("GSE140494")

## Training
trainMat <- do.call("cbind", exprsChemo_TNBC[train])
trainGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[train]))
levels(trainGroup) <- c("Sensitive", "Resistant")
table(trainGroup)

## Testing
testMat <- exprsChemo_TNBC$GSE140494
testGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[test]))
levels(testGroup) <- c("Sensitive", "Resistant")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/ChemoDataNew_GSE140494Out.rda")

##########################
# Leave GSE103668 out
train <- c("GSE25055", "GSE25065", "GSE140494", "GSE20194", "GSE20271", "GSE32646")
test <- c("GSE103668")

## Training
trainMat <- do.call("cbind", exprsChemo_TNBC[train])
trainGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[train]))
levels(trainGroup) <- c("Sensitive", "Resistant")
table(trainGroup)

## Testing
testMat <- exprsChemo_TNBC$GSE103668
testGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[test]))
levels(testGroup) <- c("Sensitive", "Resistant")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/ChemoDataNew_GSE103668Out.rda")

##########################
# Leave GSE20194 out
train <- c("GSE25055", "GSE25065", "GSE140494", "GSE103668", "GSE20271", "GSE32646")
test <- c("GSE20194")

## Training
trainMat <- do.call("cbind", exprsChemo_TNBC[train])
trainGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[train]))
levels(trainGroup) <- c("Sensitive", "Resistant")
table(trainGroup)

## Testing
testMat <- exprsChemo_TNBC$GSE20194
testGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[test]))
levels(testGroup) <- c("Sensitive", "Resistant")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/ChemoDataNew_GSE20194Out.rda")

##########################
# Leave GSE20271 out
train <- c("GSE25055", "GSE25065", "GSE140494", "GSE103668", "GSE20194", "GSE32646")
test <- c("GSE20271")

## Training
trainMat <- do.call("cbind", exprsChemo_TNBC[train])
trainGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[train]))
levels(trainGroup) <- c("Sensitive", "Resistant")
table(trainGroup)

## Testing
testMat <- exprsChemo_TNBC$GSE20271
testGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[test]))
levels(testGroup) <- c("Sensitive", "Resistant")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/ChemoDataNew_GSE20271Out.rda")

##########################
# Leave GSE32646 out
train <- c("GSE25055", "GSE25065", "GSE140494", "GSE103668", "GSE20194", "GSE20271")
test <- c("GSE32646")

## Training
trainMat <- do.call("cbind", exprsChemo_TNBC[train])
trainGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[train]))
levels(trainGroup) <- c("Sensitive", "Resistant")
table(trainGroup)

## Testing
testMat <- exprsChemo_TNBC$GSE32646
testGroup <- factor(do.call("c", ChemoSensitivity_ALL_TNBC[test]))
levels(testGroup) <- c("Sensitive", "Resistant")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/ChemoDataNew_GSE32646Out.rda")


##################################################################################
#####################
## All combined
allMat_TNBC <- do.call("cbind", exprsChemo_TNBC)
allGroup_TNBC <- unlist(ChemoSensitivity_ALL_TNBC)
allStudies_TNBC <- gsub("\\..+", "", names(allGroup_TNBC))

names(allGroup_TNBC) <- colnames(allMat_TNBC)
all(colnames(allMat_TNBC) == names(allGroup_TNBC))


####################################################################################
### Covariates for stratified sampling

#############################################################
### Grade

# Pheno1
#View(allPheno$GSE25055)
allPheno_TNBC$GSE25055$`grade:ch1` <- gsub("\\=.+", "", allPheno_TNBC$GSE25055$`grade:ch1`)
allPheno_TNBC$GSE25055$GRADE <- as.factor(allPheno_TNBC$GSE25055$`grade:ch1`)

# Pheno2
#View(allPheno$GSE25065)
allPheno_TNBC$GSE25065$GRADE <- as.factor(allPheno_TNBC$GSE25065$`grade:ch1`)

# Pheno3
#View(allPheno$GSE140494)
allPheno_TNBC$GSE140494$GRADE <- as.factor(allPheno_TNBC$GSE140494$`tumor grade:ch1`)

#Pheno4 >> No grade

# Pheno6
allPheno_TNBC$GSE20194$GRADE <- as.factor(allPheno_TNBC$GSE20194$`bmngrd:ch1`)
#allPheno_TNBC$GSE20194[allPheno_TNBC$GSE20194$GRADE == 'NA'] <- 'NA'
#levels(allPheno_TNBC$GSE20194$GRADE) <- c('2', '3', 'NA')

# Pheno7
#View(allPheno_TNBC$GSE20271)
allPheno_TNBC$GSE20271$GRADE <- as.factor(allPheno_TNBC$GSE20271$`bmn grade:ch1`)
levels(allPheno_TNBC$GSE20271$GRADE) <- c('1', '2', '3', 'NA')  
  
# Pheno8
#View(allPheno_TNBC$GSE32646)
allPheno_TNBC$GSE32646$GRADE <- as.factor(allPheno_TNBC$GSE32646$`histological grade:ch1`)

### Covariates of relevance select complete cases: GRADE
allGrade <- lapply(allPheno_TNBC, function(x) {
  i <- grep("GRADE", colnames(x))
  if (length(i) == 0) out <- factor(rep("NA", nrow(x)))
  else x <- factor(x[, i  ])
})


allGrade <- factor(unlist(allGrade))
table(allGrade)
################################################################################
########################
### T stage

# Pheno1
#View(allPheno_TNBC$GSE25055)
allPheno_TNBC$GSE25055$T_Stage <- as.factor(allPheno_TNBC$GSE25055$`clinical_t_stage:ch1`)

# Pheno2
#View(allPheno_TNBC$GSE25065)
allPheno_TNBC$GSE25065$T_Stage <- as.factor(allPheno_TNBC$GSE25065$`clinical_t_stage:ch1`)

# Pheno3
#View(allPheno_TNBC$GSE140494)
allPheno_TNBC$GSE140494$T_Stage <- as.factor(allPheno_TNBC$GSE140494$`ct stage:ch1`)

# Pheno4 > no T-stage

# Pheno6
allPheno_TNBC$GSE20194$T_Stage <- as.factor(allPheno_TNBC$GSE20194$`tbefore:ch1`)


# Pheno7
#View(allPheno_TNBC$GSE20271)
allPheno_TNBC$GSE20271$T_Stage <- as.factor(allPheno_TNBC$GSE20271$`prechemo t:ch1`)

# Pheno8
#View(allPheno_TNBC$GSE32646)
allPheno_TNBC$GSE32646$T_Stage <- as.factor(allPheno_TNBC$GSE32646$`clinical t stage:ch1`)

### Covariates of relevance select complete cases: RECURRENCE
allTstage<- lapply(allPheno_TNBC, function(x) {
  i <- grep("T_Stage", colnames(x))
  if (length(i) == 0) out <- factor(rep("NA", nrow(x)))
  else x <- factor(x[, i  ])
})

levels(allTstage$GSE25055) <- c("1", "2", "3", "4")
levels(allTstage$GSE25065) <- c("1", "2", "3", "4")
levels(allTstage$GSE140494) <- c("1", "1", "2", "2", "3", "4")
levels(allTstage$GSE103668)
levels(allTstage$GSE20194)
levels(allTstage$GSE20271) 
levels(allTstage$GSE32646) 


allTstage <- factor(unlist(allTstage))
table(allTstage)

################################################################################
########################
### N stage

# Pheno1
#View(allPheno_TNBC$GSE25055)
allPheno_TNBC$GSE25055$N_Stage <- as.factor(allPheno_TNBC$GSE25055$`clinical_nodal_status:ch1`)

# Pheno2
#View(allPheno_TNBC$GSE25065)
allPheno_TNBC$GSE25065$N_Stage <- as.factor(allPheno_TNBC$GSE25065$`clinical_nodal_status:ch1`)

# Pheno3
#View(allPheno_TNBC$GSE140494)
allPheno_TNBC$GSE140494$N_Stage <- as.factor(allPheno_TNBC$GSE140494$`cn stage:ch1`)

# Pheno4 > no N-stage

# Pheno6
allPheno_TNBC$GSE20194$N_Stage <- as.factor(allPheno_TNBC$GSE20194$`nbefore:ch1`)

# Pheno7
#View(allPheno_TNBC$GSE20271)
allPheno_TNBC$GSE20271$N_Stage <- as.factor(allPheno_TNBC$GSE20271$`prechemo n:ch1`)

# Pheno8 > no N-stage

### Covariates of relevance select complete cases: RECURRENCE
allNstage<- lapply(allPheno_TNBC, function(x) {
  i <- grep("N_Stage", colnames(x))
  if (length(i) == 0) out <- factor(rep("NA", nrow(x)))
  else x <- factor(x[, i  ])
})
# 
levels(allNstage$GSE25055) <- c("0", "1", "2", "3")
levels(allNstage$GSE25065) <- c("0", "1", "2", "3")
levels(allNstage$GSE140494) <- c("0", "1")
levels(allNstage$GSE103668)
levels(allNstage$GSE20194)
levels(allNstage$GSE20271) <- c("0", "1", "2", "3")
levels(allNstage$GSE32646) 


allNstage <- factor(unlist(allNstage))
table(allNstage)


#################################################################################
###############################
### Age
allPheno_TNBC$GSE25055$AGE <- as.numeric(allPheno_TNBC$GSE25055$`age_years:ch1`)
allPheno_TNBC$GSE25065$AGE <- as.numeric(allPheno_TNBC$GSE25065$`age_years:ch1`)
allPheno_TNBC$GSE140494$AGE <- as.numeric(allPheno_TNBC$GSE140494$`age:ch1`)
allPheno_TNBC$GSE20194$AGE <- as.numeric(allPheno_TNBC$GSE20194$`age:ch1`)
allPheno_TNBC$GSE20271$AGE <- as.numeric(allPheno_TNBC$GSE20271$`age:ch1`)
allPheno_TNBC$GSE32646$AGE <- as.numeric(allPheno_TNBC$GSE32646$`age:ch1`)

### Covariates of relevance select complete cases: AGE
allAGE <- lapply(allPheno_TNBC, function(x) {
  i <- grep("AGE", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})

allAGE <- unlist(allAGE)
summary(allAGE)

#########################################################################

### Assemble in one data.frame and turn numeric
covs <- data.frame(STUDIES=allStudies_TNBC,
                   GRADE=allGrade, 
                   T_Stage=allTstage,
                   N_Stage=allNstage, 
                   AGE=allAGE)

### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))) )

###########################################################################
###SAMPLING

### Balanced stratification
set.seed(333)
trainingOrTesting <- balancedstratification(
  covs[ , , drop=FALSE], strata=1*(allGroup_TNBC == "Sensitive"),
  pik=inclusionprobabilities(1:nrow(covs), nrow(covs) * 0.25),
  comment=TRUE, method=1)

### Show
apply(covs[, -ncol(covs),drop=FALSE], 2, table, allGroup_TNBC, trainingOrTesting)

### Subset Training
mixTrainMat <- allMat_TNBC[ , trainingOrTesting == 0]
mixTrainGroup <- allGroup_TNBC[ trainingOrTesting == 0]
mixTrainStudy <- allStudies_TNBC[ trainingOrTesting == 0]
usedTrainGrade <- as.character(allGrade[trainingOrTesting == 0])
usedTrainTstage <- as.character(allTstage[trainingOrTesting == 0])
usedTrainNstage <- as.character(allNstage[trainingOrTesting == 0])
usedTrainAge <- allAGE[trainingOrTesting == 0]

### Subset Testing
mixTestMat <- allMat_TNBC[ , trainingOrTesting == 1]
mixTestGroup <- allGroup_TNBC[ trainingOrTesting == 1]
mixTestStudy <- allStudies_TNBC[ trainingOrTesting == 1]
usedTestGrade <- as.character(allGrade[trainingOrTesting == 1])
usedTestTstage <- as.character(allTstage[trainingOrTesting == 1])
usedTestNstage <- as.character(allNstage[trainingOrTesting == 1])
usedTestAge <- allAGE[trainingOrTesting == 1]

table(mixTrainGroup)
table(mixTestGroup)

###########################################################################
### Save
save(exprsChemo_TNBC, 
     mixTrainMat, mixTrainGroup, mixTrainStudy,
     mixTestMat, mixTestGroup, mixTestStudy,
     file="./Objs/ChemoDataNew.rda")

