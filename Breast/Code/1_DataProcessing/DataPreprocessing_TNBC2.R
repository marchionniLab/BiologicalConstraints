#############################################################################
#############################################################################
## Mohamed Omar
## December 2, 2019
## Goal: Getting and organizing the data for the next step
############################################################################
rm(list = ls())

setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/BreastChemo")


library(GEOquery)
library(Biobase)
library(sampling)
library(limma)
library(genefilter)
library(edgeR)
library(survJamda)
library(caret)
library(sampling)

### Getting the data 
# dataset1 <- getGEO("GSE25055", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset1 <- dataset1$GSE25055_series_matrix.txt.gz
# 
# dataset2 <- getGEO("GSE25065", GSEMatrix = TRUE, AnnotGPL = TRUE)
# dataset2 <- dataset2$GSE25065_series_matrix.txt.gz

#save(dataset1, dataset2, file = "./Data/datasets.rda")

load("./Data/datasets.rda")
###################################################

## Get the phenotype
Pheno1 <- pData(dataset1)
Pheno2 <- pData(dataset2)

## Get the expression matrices
Expr1 <- exprs(dataset1)
Expr2 <- exprs(dataset2)

## Get the feature data
FeatData1 <- fData(dataset1)
FeatData2 <- fData(dataset2)

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
X1 <- Expr1
ffun <- pOverA(p = 1, A = 50)
Filt1 <- genefilter(2^X1, ffun)
Expr1 <- Expr1[Filt1, ]
dim(Expr1)

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
X2 <- Expr2
Filt2 <- genefilter(2^X2, ffun)
Expr2 <- Expr2[Filt2, ]
dim(Expr2)

#Expr2 <- t(scale(t(Expr2), center = TRUE, scale = TRUE))

###########################################
### Modify the phenotypes

#Keep only triple-negative Breast cancer TNBC
table(Pheno1$`er_status_ihc:ch1`)
Pheno1 <- Pheno1[Pheno1$`er_status_ihc:ch1` == "N", ]

table(Pheno1$`pr_status_ihc:ch1`)
Pheno1 <- Pheno1[Pheno1$`pr_status_ihc:ch1` == "N", ]

table(Pheno1$`her2_status:ch1`)
Pheno1 <- Pheno1[Pheno1$`her2_status:ch1` == "N", ]

# Convert to factor 
table(Pheno1$`pathologic_response_pcr_rd:ch1`)
Pheno1 <- Pheno1[!(Pheno1$`pathologic_response_pcr_rd:ch1` == "NA"), ]
Pheno1$ChemoResponse <- Pheno1$`pathologic_response_pcr_rd:ch1`
Pheno1$ChemoResponse[Pheno1$ChemoResponse == "RD"] <- "Resistant"
Pheno1$ChemoResponse[Pheno1$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno1$ChemoResponse)
Pheno1$ChemoResponse <- factor(Pheno1$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno1$ChemoResponse)

## Clean pheno1
colnames(Pheno1)
keep <- c("characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4",                   
          "characteristics_ch1.5", "characteristics_ch1.7",                    
          "characteristics_ch1.8", "characteristics_ch1.9", "characteristics_ch1.10", 
          "characteristics_ch1.13", "characteristics_ch1.14", "characteristics_ch1.15", 
          "characteristics_ch1.16", "characteristics_ch1.17",                  
          "characteristics_ch1.19", "characteristics_ch1.20", "ChemoResponse")

Pheno1 <- Pheno1[, keep]
colnames(Pheno1) <- c("Age_Years", "ER_Status", "PR_Status", "HER2_Status", "cT_Stage", 
                      "cN_Stage", "AJCC_Stage", "Grade", "DRFS_Event", "DRFS_Year", "ESR1_Status",
                      "ERBB2_Status", "SET_Class", "GGI_Class", "PAM50_Class", "ChemoResponse")


Expr1 <- Expr1[, colnames(Expr1)%in%rownames(Pheno1)]
all(rownames(Pheno1) == colnames(Expr1))

#####################
## Pheno2

# Keep only TNBC
table(Pheno2$`er_status_ihc:ch1`)
Pheno2 <- Pheno2[Pheno2$`er_status_ihc:ch1` == "N", ]

table(Pheno2$`pr_status_ihc:ch1`)
Pheno2 <- Pheno2[Pheno2$`pr_status_ihc:ch1` == "N", ]

table(Pheno2$`her2_status:ch1`)
Pheno2 <- Pheno2[Pheno2$`her2_status:ch1` == "N", ]

# 
table(Pheno2$`pathologic_response_pcr_rd:ch1`)
Pheno2 <- Pheno2[!(Pheno2$`pathologic_response_pcr_rd:ch1` == "NA"), ]
Pheno2$ChemoResponse <- Pheno2$`pathologic_response_pcr_rd:ch1`
Pheno2$ChemoResponse[Pheno2$ChemoResponse == "RD"] <- "Resistant"
Pheno2$ChemoResponse[Pheno2$ChemoResponse == "pCR"] <- "Sensitive"
table(Pheno2$ChemoResponse)
Pheno2$ChemoResponse <- factor(Pheno2$ChemoResponse, levels = c("Sensitive", "Resistant"))
levels(Pheno2$ChemoResponse)

## Clean pheno2
colnames(Pheno2)
keep <- c("characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.4",                   
          "characteristics_ch1.5", "characteristics_ch1.6", "characteristics_ch1.7",                    
          "characteristics_ch1.8", "characteristics_ch1.9", "characteristics_ch1.12",
          "characteristics_ch1.13", "characteristics_ch1.15", 
          "characteristics_ch1.16", "characteristics_ch1.17",                  
          "characteristics_ch1.19", "characteristics_ch1.20", "ChemoResponse")

Pheno2 <- Pheno2[, keep]
colnames(Pheno2) <- c("Age_Years", "ER_Status", "PR_Status", "HER2_Status", "cT_Stage", 
                      "cN_Stage", "AJCC_Stage", "Grade", "DRFS_Event", "DRFS_Year", "ESR1_Status",
                      "ERBB2_Status", "SET_Class", "GGI_Class",  "PAM50_Class", "ChemoResponse")

Expr2 <- Expr2[, colnames(Expr2)%in%rownames(Pheno2)]
all(rownames(Pheno2) == colnames(Expr2))


##################################

CommonGns <- intersect(rownames(Expr1), rownames(Expr2))

Expr1 <- Expr1[CommonGns, ]

Expr2 <- Expr2[CommonGns, ]

UsedTrainMat <- Expr1
UsedTrainGroup <- Pheno1$ChemoResponse

UsedTestMat <- Expr2
UsedTestGroup <- Pheno2$ChemoResponse

## Save the processed data
save(UsedTrainMat, UsedTestMat, UsedTrainGroup, UsedTestGroup, Pheno1, Pheno2, file = "./Objs/ChemoData2.rda")





