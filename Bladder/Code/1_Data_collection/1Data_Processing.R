#############################################################################
#############################################################################
## Mohamed Omar
## 10/4/2019
## Goal: Getting and organizing the data for the next step (Creating the K-TSP classifier)
############################################################################
rm(list = ls())

setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Bladder")


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
#dataset1 <- getGEO("GSE57813", GSEMatrix = TRUE, AnnotGPL = TRUE)
#dataset1 <- dataset1$GSE57813_series_matrix.txt.gz

#dataset2 <- getGEO("GSE39282", GSEMatrix = TRUE, AnnotGPL = TRUE)
#dataset2 <- dataset2$`GSE39282-GPL10150_series_matrix.txt.gz`
#dataset3 <- dataset2$`GSE39282-GPL15793_series_matrix.txt.gz`

#dataset4 <- getGEO("GSE37817", GSEMatrix = TRUE, AnnotGPL = TRUE)
#dataset4 <- dataset4$`GSE37817-GPL6102_series_matrix.txt.gz`
#dataset5 <- dataset4$`GSE37817-GPL8490_series_matrix.txt.gz`

#dataset6 <- getGEO("GSE19915", GSEMatrix = TRUE, AnnotGPL = TRUE)
#dataset6 <- dataset6$`GSE19915-GPL3883_series_matrix.txt.gz`
#dataset7 <- dataset6$`GSE19915-GPL4723_series_matrix.txt.gz`
#dataset8 <- dataset6$`GSE19915-GPL5186_series_matrix.txt.gz`

#save(dataset1, dataset2, file = "./Data/MyData.rda")


load("./Data/MyData.rda")

##################
## Getting old data

### Old objects
load("./Data/Old/allExprsData.rda")
load("./Data/Old/allPhenoData.rda")

### New from GEO
load("./Data/allGEOexprs.rda")
load("./Data/allGEOpheno.rda")

## Expression
expr1 <- exprs(dataset1)
#expr2 <- exprs(dataset2)

## Feature data
featData1 <- fData(dataset1)
#featData2 <- fData(dataset2)


## Annotation
rownames(expr1) <- featData1$Symbol
summary(is.na(rownames(expr1)))
rownames(expr1) <- gsub("-","", rownames(expr1))
rownames(expr1) <- gsub("_","",rownames(expr1))
sel <- which(apply(expr1, 1, function(x) all(is.finite(x)) ))
expr1 <- expr1[sel, ]
expr1 <- expr1[!is.na(rownames(expr1)),]
dim(expr1)


plot(density(expr1))

X1 <- expr1
ffun <- filterfun(pOverA(p = 0.5, A = 100))
filt1 <- genefilter(2^X1,ffun)
expr1 <- expr1[filt1,]

expr1 <- t(scale(t(expr1), center = TRUE, scale = TRUE))

#############################

# rownames(expr2) <- featData2$`Gene symbol`
# expr2 <- expr2[!(rownames(expr2) == ""), ]
# rownames(expr2) <- gsub("-","", rownames(expr2))
# rownames(expr2) <- gsub("_","",rownames(expr2))
# 
# plot(density(expr2))
# 
# X2 <- expr2
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt2 <- genefilter(2^X2,ffun)
# expr2 <- expr2[filt2,]
# dim(expr2)
# 
# expr2 <- t(scale(t(expr2), center = TRUE, scale = TRUE))


################################################

## Load expr3 raw read counts
expr3_raw <- read.delim(file = "./Data/rnaseq_gc19_extNc_chr_gtf_total_accepted_hits.htseqCounts", stringsAsFactors = FALSE)
head(expr3_raw$Gene, 20)
# remove rRNA
expr3_raw <- expr3_raw[-grep("rRNA", expr3_raw$Gene), ]
# remove tRNA
expr3_raw <- expr3_raw[-grep("tRNA", expr3_raw$Gene), ]

# Modify gene names (keep only gene IDs)
expr3_raw$Gene <- gsub(".+\\_", "", expr3_raw$Gene)

# Annotate genes (Convert IDs to gene symbols)
library(org.Hs.eg.db)
library(annotate)

tmp <- gsub("\\..*","", expr3_raw$Gene)
expr3_raw$GeneSymbol<- mapIds(org.Hs.eg.db,
                            keys=tmp,
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")

summary(!is.na(expr3_raw$GeneSymbol))
expr3_raw <- expr3_raw[!is.na(expr3_raw$GeneSymbol),]
expr3_raw$GeneSymbol <- gsub("-", "", expr3_raw$GeneSymbol)
expr3_raw$GeneSymbol <- gsub("_", "", expr3_raw$GeneSymbol)


expr3_raw <- expr3_raw[!duplicated(expr3_raw$GeneSymbol),]
rownames(expr3_raw) <- expr3_raw$GeneSymbol

expr3_raw$Gene <- NULL
expr3_raw$GeneSymbol <- NULL


## Filter expr3 by keeping only the genes with cpm > 1 (raw counts > 15) in at least 50% of samples
mycpm <- cpm(expr3_raw)
#plot(expr3_raw[,1],mycpm[,1],xlim=c(0,40),ylim=c(0,3))
#abline(v=15,col=2)
#abline(h=1,col=4)

thresh <- mycpm > 1
keep <- rowSums(thresh) >= 238
table(keep)

expr3_raw <- expr3_raw[keep,]
dim(expr3_raw)

## Visualize density before/after filtering
#plot(density(log2(as.matrix(expr3_raw))))
#plot(density(log2(as.matrix(expr3))))


expr3 <- DGEList(expr3_raw)
#barplot(expr3$samples$lib.size, names.arg = colnames(expr3), las=2)


expr3 <- calcNormFactors(expr3, method = c("TMM"))
#expr3$samples

plotMD(expr3,column=2)
abline(h=0,col="grey")

expr3 <- cpm(expr3, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
#boxplot(expr3[,1:100])

expr3 <- t(scale(t(expr3), center = TRUE, scale = TRUE))
###############################################

### Phenotype
pheno1 <- pData(dataset1)

#pheno2 <- pData(dataset2)

pheno3 <- read.csv(file = "./Data/E-MTAB-4321.sdrf.txt", sep = "\t")

### Modify the pheno

pheno1$PROGRESSION <- as.factor(pheno1$`liver sample group:ch1`)
levels(pheno1$PROGRESSION) <- c("NoProgression", "Progression")

# pheno2/pheno3 ## Not good (GSE39282)

# pheno2$PROGRESSION <- as.factor(pheno2$`progression:ch1`)
# levels(pheno2$PROGRESSION) <- c("NoProgression", "Progression")
# # Remove normal samples
# pheno2 <- pheno2[-c(1:6), ]
# # Modify expr2
# expr2 <- expr2[, colnames(expr2) %in% rownames(pheno2)]
# all(rownames(pheno2) == colnames(expr2))
# 


## Remove Ta from pheno3
pheno3 <- pheno3[which(!grepl("Ta", pheno3$Factor.Value.disease.staging.)), ]
pheno3$PROGRESSION <- as.factor(pheno3$Factor.Value.progression.to.T2..)
levels(pheno3$PROGRESSION) <- c("NoProgression", "Progression")
rownames(pheno3) <- pheno3$Source.Name
# Modify expr3
expr3 <- expr3[,colnames(expr3) %in% rownames(pheno3)]
all(rownames(pheno3) == colnames(expr3))

## Clean Pheno for Luigi
# keep <- c("Source.Name", "Characteristics.organism.", "Characteristics.disease.", "Characteristics.tumor.grading.", "Characteristics.disease.staging.",
#           "Characteristics.sex.", "Characteristics.age.", "Characteristics.tumor.size.", "Characteristics.growth.pattern.", "Characteristics.BCG.treatment.", 
#            "Characteristics.cystectomy.", "Characteristics.progression.free.survival.", "PROGRESSION")

# pheno3 <- pheno3[, keep]
# colnames(pheno3) <- c("SampleName", "Organism", "Disease", "Grade", "Stage", "Sex", "Age", "Size_cm", "Growth_Pattern", "BCG_TTT", "Cystectomy", "PFS", "Progression")
# 
# pheno3 <- apply(pheno3, 2, as.character)
# rownames(pheno3) <- pheno3[, "SampleName"]
# 
# ###writing the data to a file
# write.table(pheno3, file="./Objs/ForLuigi/EMTAB4321/Pheno.txt", sep="\t", row.names=T, col.names=T)
# 
# #save expr3
# write.table(expr3, file="./Objs/ForLuigi/EMTAB4321/Expr.txt", sep="\t", row.names=T, col.names=T)
############

allpheno1 <- list(pheno1, pheno3)
names(allpheno1) <- c("GSE57813", "E-MTAB-4321")

allExpr1 <- list(expr1, expr3)
names(allExpr1) <- c("GSE57813", "E-MTAB-4321")


### Filter phenotype information for the required samples
groupProgression1 <- mapply(x=allpheno1, FUN=function(x) {
  x <- x[,"PROGRESSION"]
  out <- factor(x, levels=c("NoProgression", "Progression"))
  out
})
###################################################################################
## Process old data

### Combine expression and phenotypes (OLD DATA)
allExprs2 <- c(allExprsData, allGEOexprs)
allPheno2 <- c(allPhenoData, allGEOpheno)


plot(density(allExprs2$GSE13507))
dim(allExprs2$GSE13507)
X3 <- allExprs2$GSE13507
ffun <- filterfun(pOverA(p = 0.5, A = 100))
filt3 <- genefilter(2^X3,ffun)
allExprs2$GSE13507 <- allExprs2$GSE13507[filt3,]
#

plot(density(allExprs2$pmid15930337))
dim(allExprs2$pmid15930337)
X4 <- allExprs2$pmid15930337
ffun <- filterfun(pOverA(p = 0.5, A = 100))
filt4 <- genefilter(2^X4,ffun)
allExprs2$pmid15930337 <- allExprs2$pmid15930337[filt4,]


# Z-score transformation
allExprs2$GSE13507 <- t(scale(t(allExprs2$GSE13507), center = TRUE, scale = TRUE))
allExprs2$pmid15930337 <- t(scale(t(allExprs2$pmid15930337), center = TRUE, scale = TRUE))
allExprs2$GSE32894 <- allExprs2$GSE32894/sd(allExprs2$GSE32894)

### Prepare objects
allPheno2 <- lapply(allPheno2, function(x) {
    colnames(x) <- toupper(colnames(x))
    colnames(x) <- gsub("\\.$", "", gsub("\\.\\.", ".", colnames(x)))
    colnames(x) <- gsub("\\.", "_", colnames(x))
    x
})

### Retain only datasets with survival in caucasians
keep <- c("GSE13507", "pmid15930337", "GSE32894")

### Drop unwanted datasets and format survival time
allPheno2 <- allPheno2[keep]
allExprs2 <- allExprs2[keep]

### Process information for consistency
allPheno2 <- lapply(allPheno2, function(x) {
    x[ , grep("_MO", colnames(x)) ] <- sapply(x[ , grep("_MO", colnames(x)) ], as.numeric)
    x[ , grep("DATEFACT", colnames(x)) ] <- sapply(x[ , grep("DATEFACT", colnames(x)) ], as.numeric)
    x[ , sapply(x, class) == "factor" ] <- sapply(x[ , sapply(x, class) == "factor"], as.character)
    x
})

### Make IDs as rownames
rownames(allPheno2$GSE13507) <- allPheno2$GSE13507$GSM
rownames(allPheno2$pmid15930337) <- allPheno2$pmid15930337$IDS
rownames(allPheno2$GSE32894) <- allPheno2$GSE3289$GEO_ACCESSION

###########################################################################
### Keep only relevant samples
keepSamples <- sapply(allPheno2, function(x) {
  if (length(grep("STAGE$", colnames(x))) > 0) {
    out <- x[ , grep("STAGE$", colnames(x))] %in% c("T1", "STAGE:T1N0M0") &
      !is.na(x[ , grep("STAGE$", colnames(x)) ]) &
      !is.na(x[ , grep("PROGR", colnames(x)) ])
  } else {
    out <- rep(TRUE, nrow(x))
  }
  out
})


### Filter phenotype information for the required samples
groupProgression2 <- mapply(x=allPheno2, y=keepSamples, FUN=function(x, y) {
  nms <- rownames(x)[y]
  x <- x[y, grep("PROGRESSION", colnames(x))]
  names(x) <- nms
  x[grep("YES", x)] <- "YES"
  x[x %in% c("1", "NO", "NonProg")] <- "NoProgression"
  x[x %in% c("2", "YES", "Prog")] <- "Progression"
  out <- factor(x, levels=c("NoProgression", "Progression"))
  out
})

### Modify allExprs2 to include only the relevant samples
allExprs2$GSE13507 <- allExprs2$GSE13507[,colnames(allExprs2$GSE13507)%in% names(groupProgression2$GSE13507)]
allExprs2$pmid15930337 <- allExprs2$pmid15930337[, colnames(allExprs2$pmid15930337)%in% names(groupProgression2$pmid15930337)]
allExprs2$GSE32894 <- allExprs2$GSE32894[, colnames(allExprs2$GSE32894) %in% names(groupProgression2$GSE32894)]
#######################

## For Luigi/Don

# View(allPheno2$GSE13507)
# pheno_GSE13507 <- allPheno2$GSE13507
# expr_GSE13507 <- allExprs2$GSE13507
# 
# 
# pheno_GSE13507 <- pheno_GSE13507[!is.na(pheno_GSE13507$PROGRESSION), ]
# pheno_GSE13507 <- pheno_GSE13507[(pheno_GSE13507$STAGE == "T1"), ]
# pheno_GSE13507$PROGRESSION <- as.factor(pheno_GSE13507$PROGRESSION)
# table(pheno_GSE13507$PROGRESSION)
# levels(pheno_GSE13507$PROGRESSION) <- c("NoProgression", "Progression")
# 
# # Clean
# colnames(pheno_GSE13507)
# keep <- c("GSM", "TYPE", "SEX", "AGE", "INVASIVENESS", "INTRAVESICAL_THERAPY", 
#           "SYSTEMIC_CHEMO", "STAGE", "GRADE", "RECURRENCE", "PROGRESSION", "OVERALL_SURVIVAL", 
#           "SURVIVALMONTH", "METS", "NODES", "STAGE_SIMPLE")
# pheno_GSE13507 <- pheno_GSE13507[, keep]
# pheno_GSE13507$NODES <- NULL
# colnames(pheno_GSE13507) <- c("SampleName", "Type", "Sex", "Age", "Invasiveness", "Intravesical_Therapy", "Systemic_Chemo", "Stage", "Grade", "Recurrence",
#                               "Progression", "OS", "Survival_Month", "Mets", "Stage")
# 
# all(rownames(pheno_GSE13507) == colnames(expr_GSE13507))
# 
# # SAve for Luigi
# pheno_GSE13507 <- apply(pheno_GSE13507, 2, as.character)
# rownames(pheno_GSE13507) <- pheno_GSE13507[, "SampleName"]
# 
# table(pheno_GSE13507[,"Progression"])
# ###writing the data to a file
# write.table(pheno_GSE13507, file="./Objs/ForLuigi/GSE13507/Pheno.txt", sep="\t", row.names=T, col.names=T)
# 
# #save expr3
# write.table(expr_GSE13507, file="./Objs/ForLuigi/GSE13507/Expr.txt", sep="\t", row.names=T, col.names=T)

######################
# View(allPheno2$GSE32894)
# pheno_GSE32894 <- allPheno2$GSE32894
# expr_GSE32894 <- allExprs2$GSE32894
# 
# 
# pheno_GSE32894 <- pheno_GSE32894[(pheno_GSE32894$TUMOR_STAGE == "T1"), ]
# pheno_GSE32894$PROGRESSION <- as.factor(pheno_GSE32894$PROGRESSION_TO_HIGHER_STAGE_OR_RADIOTHERAPY)
# table(pheno_GSE32894$PROGRESSION)
# levels(pheno_GSE32894$PROGRESSION) <- c("NoProgression", "Progression", "Progression", "Progression", "Progression", "Progression")
# pheno_GSE32894 <- pheno_GSE32894[!is.na(pheno_GSE32894$PROGRESSION), ]
# 
# # Clean
# colnames(pheno_GSE32894)
# keep <- c("ROW_NAMES", "GENDER", "AGE", "TUMOR_STAGE", "TUMOR_GRADE", "TIME_TO_DOD_MONTHS",
#           "DOD_EVENT", "PROGRESSION")
# pheno_GSE32894 <- pheno_GSE32894[, keep]
# colnames(pheno_GSE32894) <- c("SampleName", "Sex", "Age", "Stage", "Grade", "DOD_Month", "DOD_Event", "Progression")
# 
# all(rownames(pheno_GSE32894) == colnames(expr_GSE32894))
# 
# # SAve for Luigi
# pheno_GSE32894 <- apply(pheno_GSE32894, 2, as.character)
# rownames(pheno_GSE32894) <- pheno_GSE32894[, "SampleName"]
# 
# table(pheno_GSE32894[,"Progression"])
# ###writing the data to a file
# write.table(pheno_GSE32894, file="./Objs/ForLuigi/GSE32894/Pheno.txt", sep="\t", row.names=T, col.names=T)
# 
# #save expr3
# write.table(expr_GSE32894, file="./Objs/ForLuigi/GSE32894/Expr.txt", sep="\t", row.names=T, col.names=T)

#########################
###################################################################
### Combine old and new
allPheno <- c(allPheno2, allpheno1)
allExpr <- c(allExprs2, allExpr1)
groupProgression_ALL <- c(groupProgression2, groupProgression1)



### Keep good expression data and phenotypes
allExpr <- allExpr[ names(groupProgression_ALL) ]

### Find commom subset of genes
commonGenes <- Reduce("intersect", lapply(allExpr, rownames))


### Filter expression for the required samples
exprsProgression <- mapply(x=allExpr, FUN=function(x, gns) {
  x <- x[ gns ,]
}, MoreArgs=list(gns=commonGenes))

### Check
all(names(exprsProgression) == names(groupProgression_ALL))

### Check order
all(rownames(allPheno$GSE57813) == colnames(allExpr$GSE57813))
all(rownames(allPheno$`GSE37817-GPL6102`) == colnames(allExpr$`GSE37817-GPL6102`))
all(rownames(allPheno$`E-MTAB-4321`) == colnames(allExpr$`E-MTAB-4321`))


###################################################################
#############
## Cross-study validation
# Leave GSE13507 out
train <- c("pmid15930337","GSE32894", "GSE57813", "E-MTAB-4321")
test <- c("GSE13507")

## Training
trainMat <- do.call("cbind", exprsProgression[train])
trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
levels(trainGroup) <- c("NoProgression", "Progression")
table(trainGroup)

## Testing
testMat <- exprsProgression$GSE13507
testGroup <- factor(do.call("c", groupProgression_ALL[test]))
levels(testGroup) <- c("NoProgression", "Progression")
table(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/progressionDataGood_GSE13507Out.rda")

#############
# Leave pmid15930337 out
train <- c("GSE13507","GSE32894", "GSE57813", "E-MTAB-4321")
test <- c("pmid15930337")

## Training
trainMat <- do.call("cbind", exprsProgression[train])
trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
levels(trainGroup) <- c("NoProgression", "Progression")
table(trainGroup)
length(trainGroup)

## Testing
testMat <- exprsProgression$pmid15930337
testGroup <- factor(do.call("c", groupProgression_ALL[test]))
levels(testGroup) <- c("NoProgression", "Progression")
table(testGroup)
length(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/progressionDataGood_pmidOut.rda")

#############
# Leave GSE32894 out
train <- c("GSE13507","pmid15930337", "GSE57813", "E-MTAB-4321")
test <- c("GSE32894")

## Training
trainMat <- do.call("cbind", exprsProgression[train])
trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
levels(trainGroup) <- c("NoProgression", "Progression")
table(trainGroup)
length(trainGroup)

## Testing
testMat <- exprsProgression$GSE32894
testGroup <- factor(do.call("c", groupProgression_ALL[test]))
levels(testGroup) <- c("NoProgression", "Progression")
table(testGroup)
length(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/progressionDataGood_GSE32894Out.rda")

#############
# Leave GSE57813 out
train <- c("GSE13507","pmid15930337", "GSE32894", "E-MTAB-4321")
test <- c("GSE57813")

## Training
trainMat <- do.call("cbind", exprsProgression[train])
trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
levels(trainGroup) <- c("NoProgression", "Progression")
table(trainGroup)
length(trainGroup)

## Testing
testMat <- exprsProgression$GSE57813
testGroup <- factor(do.call("c", groupProgression_ALL[test]))
levels(testGroup) <- c("NoProgression", "Progression")
table(testGroup)
length(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/progressionDataGood_GSE57813Out.rda")

#############
# Leave EMTAB-4321 out
train <- c("GSE13507","pmid15930337", "GSE32894", "GSE57813")
test <- c("E-MTAB-4321")

## Training
trainMat <- do.call("cbind", exprsProgression[train])
trainGroup <- factor(do.call("c", groupProgression_ALL[train]))
levels(trainGroup) <- c("NoProgression", "Progression")
table(trainGroup)
length(trainGroup)

## Testing
testMat <- exprsProgression$`E-MTAB-4321`
testGroup <- factor(do.call("c", groupProgression_ALL[test]))
levels(testGroup) <- c("NoProgression", "Progression")
table(testGroup)
length(testGroup)

# Save for cross-study validation
save(trainMat, trainGroup, testMat, testGroup, file = "./Objs/progressionDataGood_EMTABOut.rda")

##################################################################
#####################
## All combined
allMat <- do.call("cbind", exprsProgression)
allGroup <- unlist(groupProgression_ALL)
allStudies <- gsub("\\..+", "", names(allGroup))

names(allGroup) <- colnames(allMat)
all(colnames(allMat) == names(allGroup))


### Prepare to extract additional info
myPheno1 <- mapply(function(x, y) x[y,], allPheno[1:3], keepSamples)
myPheno <- c(myPheno1, allpheno1)

#############################################################
### Grade

## allpheno$GSE57813$Grade <- allpheno$GSE57813$.. ##No grade
#myPheno$`GSE37817-GPL6102`$GRADE <- as.factor(myPheno$`GSE37817-GPL6102`$`grade:ch1`)
myPheno$`E-MTAB-4321`$GRADE <- as.factor(myPheno$`E-MTAB-4321`$Factor.Value.tumor.grading.)



### Covariates of relevance select complete cases: GRADE
allGrade <- lapply(myPheno, function(x) {
  i <- grep("GRADE", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})


levels(allGrade$GSE13507) <- c("low", "high")
levels(allGrade$pmid15930337) <- c("low", "high")
levels(allGrade$GSE32894) <- c("low", "high", "high")
#levels(allGrade$GSE57813) <- c("high", "low")
levels(allGrade$`E-MTAB-4321`) <- c("high", "low")


allGrade <- factor(unlist(allGrade))
levels(allGrade)[levels(allGrade) == "GX"] <- NA
levels(allGrade)[levels(allGrade) == ""] <- NA

allGrade <- ordered(allGrade, levels=c("low", "high"))

################################################################################
########################
### Recurrence

#myPheno$`GSE37817-GPL6102`$RECURRENCE <- as.factor(myPheno$`GSE37817-GPL6102`$`recurrence:ch1`)
#levels(myPheno$`GSE37817-GPL6102`$RECURRENCE) <- c("NoRecurrence", "Recurrence")

### Covariates of relevance select complete cases: RECURRENCE
# allRecurrence <- lapply(myPheno, function(x) {
#   i <- grep("RECURRENCE", colnames(x))
#   if (length(i) == 0) out <- factor(rep("", nrow(x)))
#   else x <- factor(x[, i  ])
# })
# 
# levels(allRecurrence$GSE13507) <- c("Recurrence", "NoRecurrence") # Value 1=Rec, value 2=NoRec
# 
# allRecurrence <- factor(unlist(allRecurrence))
# levels(allRecurrence)[levels(allRecurrence) == ""] <- NA

###############################################################################
############################
### Intravesical therapy

#myPheno$`GSE37817-GPL6102`$INTRAVESICAL <- as.factor(myPheno$`GSE37817-GPL6102`$`intravesical therapy:ch1`)
myPheno$`E-MTAB-4321`$INTRAVESICAL <- as.factor(myPheno$`E-MTAB-4321`$Factor.Value.BCG.treatment.)


### Covariates of relevance select complete cases: INTRAVESCICAL
allRX <- lapply(myPheno, function(x) {
  i <- grep("INTRAVESICAL", colnames(x))
  if (length(i) == 0) out <- factor(rep("", nrow(x)))
  else x <- factor(x[, i  ])
})

levels(allRX$GSE13507) <- c("YES", "NO") # Value 1=YES, value 2=NO
#levels(allRX$`GSE37817-GPL6102`) <- c("NO", "YES")
levels(allRX$`E-MTAB-4321`) <- c("NO", "YES")

allRX <- factor(unlist(allRX))
levels(allRX)[levels(allRX) == ""] <- NA

#################################################################################
###############################
### Age
#myPheno$`GSE37817-GPL6102`$AGE <- as.character(myPheno$`GSE37817-GPL6102`$`AGE:ch1`)
myPheno$`E-MTAB-4321`$AGE <- as.character(myPheno$`E-MTAB-4321`$Factor.Value.age.)

### Covariates of relevance select complete cases: AGE
allAGE <- lapply(myPheno, function(x) {
  i <- grep("^AGE$", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- as.numeric(x[, i  ])
})
allAGE <- unlist(allAGE)

##################################################################################
################################
### Sex
#myPheno$`GSE37817-GPL6102`$GENDER <- as.factor(myPheno$`GSE37817-GPL6102`$`SEX:ch1`)
myPheno$`E-MTAB-4321`$GENDER <- as.factor(myPheno$`E-MTAB-4321`$Factor.Value.sex.)
myPheno$GSE13507$GENDER <- myPheno$GSE13507$SEX

### Covariates of relevance select complete cases: SEX
allGENDER <- lapply(myPheno, function(x) {
  i <- grep("GENDER", colnames(x))
  if (length(i) == 0) out <- rep(NA, nrow(x))
  else x <- factor(x[, i  ])
})
levels(allGENDER$`E-MTAB-4321`) <- c("F", "M")
allGENDER <- factor(unlist(allGENDER))


#########################################################################

### Assemble in one data.frame and turn numeric
covs <- data.frame(STUDIES=allStudies,
                   GRADE=allGrade, RX=allRX,
                   GENDER=allGENDER, AGE=allAGE)

### Prepare vocs for sampling
covs <- sapply(covs , function(x) as.numeric(factor(paste(x))) )


###########################################################################
###SAMPLING

### Balanced stratification
set.seed(333)
trainingOrTesting <- balancedstratification(
  covs[ , , drop=FALSE], strata=1*(allGroup == "Progression"),
  pik=inclusionprobabilities(1:nrow(covs), 90),
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
save(exprsProgression, 
     mixTrainMat, mixTrainGroup, mixTrainStudy,
     mixTestMat, mixTestGroup, mixTestStudy,
     file="./Objs/progressionDataGood2.rda")

#########################################################################
#########################################################################
########################################################################









