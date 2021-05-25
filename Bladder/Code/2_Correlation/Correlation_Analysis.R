############################################################################
############################################################################
### August 1, 2019
### Luigi Marchionni
### Analyze expression FANTOM5 data for the Zoo project
### Collaboration with Michiel DeHoon
### Identify conserved genes across evolution using ICOR


############################################################################
### Clean working space ans set working directory
############################################################################

### Clean
rm(list=ls())

### Set working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Bladder")

### Set global options
options(width=130)


############################################################################
### Load libraries
library(Biobase)
library(MergeMaid)
library(parallel)
library(mclust)
library(limma)
library(dplyr)
library(edgeR)
library(genefilter)
###########################################################################
### Useful functions
###########################################################################

### Functions
source("/Users/mohamedomar/Documents/Research/Tutorials/plotEMdistr.R")

### Panel correlations
panel.cor <- function(x, y, digits = 3, prefix = "Correlation:\n", cex.cor=1)  {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}

### scatterPlot  with regression line
panel.scatterFit <- function(x, y) {
  points(x, y, pch=16)
  abline(lsfit(x,y), col="blue", lty=2)
}

### Geometric mean
geoMean <- function(x, ...) exp(mean(log(x), na.rm=TRUE, ...))

### Nice histogram
niceHist <- function(x, probability = TRUE, nclass=25,
                     xlab="Overall Integrative Correlation", ylab="Density", ...) {
  hist(x, probability = probability, col="salmon",
       ylab=ylab, nclass=nclass, xlab=xlab, ...)
  lines(density(x), lwd=2, col="blue")
  abline(v=median(x), lwd=2, lty=2, col="blue")
}


###########################################################################
### Read/Load data
###########################################################################
load("./Data/MyData.rda")

##################
## Getting old data

### Old objects
load("./Data/Old/allExprsData.rda")
load("./Data/Old/allPhenoData.rda")

### New from GEO
load("./Data/allGEOexprs.rda")
load("./Data/allGEOpheno.rda")

## Get the expression matrices (of the new datasets)
expr1 <- exprs(dataset1)
#expr2 <- exprs(dataset2)
expr3_raw <- read.delim(file = "./Data/rnaseq_gc19_extNc_chr_gtf_total_accepted_hits.htseqCounts", stringsAsFactors = FALSE)


## Feature data
featData1 <- fData(dataset1)
#featData2 <- fData(dataset2)

#############################################################################

## Annotation of expr1
rownames(expr1) <- featData1$Symbol
summary(is.na(rownames(expr1)))
rownames(expr1) <- gsub("-","", rownames(expr1))
rownames(expr1) <- gsub("_","",rownames(expr1))
sel <- which(apply(expr1, 1, function(x) all(is.finite(x)) ))
expr1 <- expr1[sel, ]
expr1 <- expr1[!is.na(rownames(expr1)),]
dim(expr1)

# Group by the duplicate gene names (take the mean)
# expr1 <- as.data.frame(expr1)
# expr1 %>% 
#   group_by(rownames(expr1)) %>% 
#   summarise_all(funs(median))

X1 <- expr1
ffun <- filterfun(pOverA(p = 0.5, A = 100))
filt1 <- genefilter(2^X1,ffun)
expr1 <- expr1[filt1,]

#expr1 <- t(scale(t(expr1), center = TRUE, scale = TRUE))

###############################
## Annotation of expr2
# rownames(expr2) <- featData2$`Gene symbol`
# expr2 <- expr2[!(rownames(expr2) == ""), ]
# rownames(expr2) <- gsub("-","", rownames(expr2))
# rownames(expr2) <- gsub("_","",rownames(expr2))
# 
# ## group the duplicate gene symbols (take the mean)
# # expr2 <- as.data.frame(expr2)
# # expr2 %>% 
# #   group_by(rownames(expr2)) %>% 
# #   summarise_all(funs(median))
# 
# X2 <- expr2
# ffun <- filterfun(pOverA(p = 0.5, A = 100))
# filt2 <- genefilter(2^X2,ffun)
# expr2 <- expr2[filt2,]
# # 
# expr2 <- t(scale(t(expr2), center = TRUE, scale = TRUE))
# ###############################
# ## Annotation of expr3

## Load expr3 (RNA-seq) raw read counts
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
# plot(expr3_raw[,1],mycpm[,1],xlim=c(0,40),ylim=c(0,3))
# abline(v=15,col=2)
# abline(h=1,col=4)
# 
thresh <- mycpm > 1
keep <- rowSums(thresh) >= 238
table(keep)
# 
expr3_raw <- expr3_raw[keep,]
dim(expr3_raw)
# 
# ## Visualize density before/after filtering
# plot(density(log2(as.matrix(expr3_raw))))
# plot(density(log2(as.matrix(expr3))))

# Convert to DGE object
expr3 <- DGEList(expr3_raw)
barplot(expr3$samples$lib.size, names.arg = colnames(expr3), las=2)

# Calculate the normalization factors
expr3 <- calcNormFactors(expr3, method = c("TMM"))
#expr3$samples

plotMD(expr3,column=2)
abline(h=0,col="grey")

# Take the CPM accounting for the lib.sizes
expr3 <- cpm(expr3, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
#boxplot(expr3[,1:100])

# Group the duplicate gene symbols (take the mean)
# expr3 <- as.data.frame(expr3)
# expr3 %>% 
#   group_by(rownames(expr3)) %>% 
#   summarise_all(funs(median))

#expr3 <- t(scale(t(expr3), center = TRUE, scale = TRUE))

##############################
## Collect all the expr of the new datasets
allExpr_New <- list(expr1, expr3)
names(allExpr_New) <- c("GSE57813", "E-MTAB-4321")

###################################################################################
#########################################################################

## Process old data

### Combine the expression of the OLD DATA
allExprs_Old <- c(allExprsData, allGEOexprs)

# Keep the relevant datasets
keep <- c("GSE13507", "pmid15930337", "GSE32894")
allExprs_Old <- allExprs_Old[keep]

####################
## Group the duplicate gene symbols (take the mean)
# allExprs_Old$GSE13507 <- as.data.frame(allExprs_Old$GSE13507)
# allExprs_Old$GSE13507 %>% 
#   group_by(rownames(allExprs_Old$GSE13507)) %>% 
#   summarise_all(funs(median))

dim(allExprs_Old$GSE13507)
X3 <- allExprs_Old$GSE13507
ffun <- filterfun(pOverA(p = 0.5, A = 100))
filt3 <- genefilter(2^X3,ffun)
allExprs_Old$GSE13507 <- allExprs_Old$GSE13507[filt3,]

#allExprs_Old$GSE13507 <- t(scale(t(allExprs_Old$GSE13507), center = TRUE, scale = TRUE))

########
# allExprs_Old$pmid15930337 <- as.data.frame(allExprs_Old$pmid15930337)
# allExprs_Old$pmid15930337 %>% 
#   group_by(rownames(allExprs_Old$pmid15930337)) %>% 
#   summarise_all(funs(median))


dim(allExprs_Old$pmid15930337)
X4 <- allExprs_Old$pmid15930337
ffun <- filterfun(pOverA(p = 0.5, A = 100))
filt4 <- genefilter(2^X4,ffun)
allExprs_Old$pmid15930337 <- allExprs_Old$pmid15930337[filt4,]

#allExprs_Old$pmid15930337 <- t(scale(t(allExprs_Old$pmid15930337), center = TRUE, scale = TRUE))

########
# allExprs_Old$GSE32894 <- as.data.frame(allExprs_Old$GSE32894)
# allExprs_Old$GSE32894 %>% 
#   group_by(rownames(allExprs_Old$GSE32894)) %>% 
#   summarise_all(funs(median))

#allExprs_Old$GSE32894 <- allExprs_Old$GSE32894/sd(allExprs_Old$GSE32894)

##############################################################
## Combine the old and the new expression matrices in a single list
allExpr <- c(allExprs_Old, allExpr_New)


### Quantile normalize
allExpr <- lapply(allExpr, normalizeBetweenArrays)

### Check dimensions
sapply(allExpr, dim)


###########################################################################
### Prepapre object for MergeMaid analysis
###########################################################################

###########################################################################
### Uncollapsed data
### Get matrices out of the list (necessary to work with MergeMaid)
expr1 <- as.matrix(allExpr$GSE57813)
#expr2 <- as.matrix(allExpr$`GSE37817-GPL6102`)
expr3 <- as.matrix(allExpr$`E-MTAB-4321`)
expr4 <- as.matrix(allExpr$GSE13507)
expr5 <- as.matrix(allExpr$pmid15930337)
expr6 <- as.matrix(allExpr$GSE32894)

expr1 <- expr1[!duplicated(rownames(expr1)), ]
#expr2 <- expr2[!duplicated(rownames(expr2)), ]
expr3 <- expr3[!duplicated(rownames(expr3)), ]
expr4 <- expr4[!duplicated(rownames(expr4)), ]
expr5 <- expr5[!duplicated(rownames(expr5)), ]
expr6 <- expr6[!duplicated(rownames(expr6)), ]

### Merge data and clean object names
mgdExpr <- mergeExprs(expr1, expr3, expr4, expr5, expr6)



###########################################################################
### Now calculate the correlations
###########################################################################

### Calculate integrative correlation: original values
iCorP <- intCor(mgdExpr, method="pearson", exact = FALSE)

###########################################################################
### Explore reproducibility using plots
###########################################################################

###########################################################################
### ALL DATA

### Open device
png("./Figs/Correlation/iCorDensPearson.png", width=2000, height=2000, res=200)
### Set graphic layout
par(oma=rep(0.5,4), mar=c(2.75,2.75,0.75,0.25), mgp=c(1.75, 0.75, 0))
### Compute and plot distributions
set.seed(333)
iCorNullP <- intcorDens(mgdExpr, cex.legend=1.1, lwd=2)
### Close device
dev.off() 

### Open device
png("./Figs/Correlation/iCorHistPearson.png", width=2000, height=2000, res=200)
### Set graphic layout
par(oma=rep(0.25,4), mar=c(2.75,2.75,0.75,0.25), mgp=c(1.75, 0.75, 0))
### Compute and plot distributions
### hist(iCorP, col="salmon")
niceHist(integrative.cors(iCorP), main="Integrative Correlation")
### Close device
dev.off()


###########################################################################

###########################################################################
### EM-algorithm to split genes based on iCOR
###########################################################################
set.seed(333)

null <- unlist(iCorNullP)
iCorNullP_Modified <- na.omit(sample(null, length(integrative.cors(iCorP)), replace = TRUE))
### Create list of integrative correlations omit NAs
myDat <- data.frame("Corr_All"=na.omit(integrative.cors(iCorP)),
                    "Null_Distribution_All"= iCorNullP_Modified)
                    

### ### Use Mclust with 3 normals using mclapply from the "parallel" package
datEM <- mclapply(myDat, Mclust, mc.cores=detectCores(), G=1:2)


###########################################################################
### Plot
png("./Figs/Correlation/iCorPearsonEMclass.png", height=3500, width=3500, res=350)

### Set layout
nc <- floor(sqrt(ncol(myDat)))
nr <- ceiling((ncol(myDat))/nc)
layout(matrix(1:(nc*nr), ncol=nc, nrow= nr, byrow=FALSE))
### Set par
par(list(mar=c(3, 3, 4, 1), mgp=c(1.5, 0.5, 0)))
### Plot the distributions
splitPoint <- mapply(exprDat=myDat, emOut=datEM, name=gsub("_mRNA", "", colnames(myDat)),
                     FUN=plotEMdistr, MoreArgs=list(quant = c(0.1, 0.9)),
                     forXlab="Integrative Correlation", forTitle="Integrative Correlation")
### Close screen and device
dev.off()


##########################################################################
###########################################################################
### Reproducible genes: Pearson
iCorPint <- integrative.cors(iCorP)
### More then selected quantile
thr <- 0.15
#thr <- quantile(iCorPint, 0.66)

table(iCorPint > thr)
iCorPint_T <- iCorPint[iCorPint > thr]

RGenes <- names(iCorPint_T)
save(iCorPint, RGenes, file = "./Objs/Correlation/RGenes.rda")


###########################################################################
###########################################################################
### Find split point for each pair-wise correlation in comparison to null

myDatWithNull <- data.frame(
  ## Join real and null data to run EM-classification for ALL
  Corr_All = c(myDat$Corr_All, myDat$Null_Distribution_All))


### ### Use Mclust with 3 normals using mclapply from the "parallel" package
datEM <- mclapply(myDatWithNull, Mclust, mc.cores=detectCores(), G=1:2)


###########################################################################
### Plot
png("figs/iCorPearsonEMclassNULLvsREAL.png", height=2500, width=3500, res=300)
### Set layout
nc <- floor(sqrt(ncol(myDat)))
nr <- ceiling((ncol(myDat))/nc)
layout(matrix(1:(nc*nr), ncol=nc, nrow= nr, byrow=TRUE))
### Set par
par(list(mar=c(3, 3, 2.5, 1), mgp=c(1.5, 0.5, 0)))
### Plot the distributions
splitPointVSnull <- mapply(exprDat=myDatWithNull, emOut=datEM, name=gsub("_mRNA", "", colnames(myDat)),
                           FUN=plotEMdistr, MoreArgs=list(quant = c(0.1, 0.9),  upper=2.1),
                           forXlab="Integrative Correlation", forTitle="Integrative Correlation")
### Close screen and device
dev.off()


###########################################################################
### Save
save(eDat, eDatMed, eDatMean, eDatQuant, mgdExpr, mgdExprMed, mgdExprMean, mgdExprQuant,
     iCorP, iCorMedP, iCorMeanP,  iCorQuantP,  iCorNullP, iCorNullMedP, iCorNullMeanP, iCorNullQuantP,
     splitPoint, splitPointVSnull, file="objs/iCorData.rda")


###########################################################################
### Session information
sessionInfo()

### Clean working space and quit
rm(list=ls())
q("no")
