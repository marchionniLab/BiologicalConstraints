

rm(list = ls())

library(tidyverse)
library(data.table)
library(limma)

#######################################
## Load the data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#######################################
# List of stromal genes
StromalGns_up <- c('AEBP1', 'ANTXR1', 'BGN', 'C1QA', 'C1QB', 'C1QC', 'C1R', 'C1S', 'CDH11', 'COL1A1', 'COL3A1', 'FBLN5', 'FCGR2B', 'HLA-DRB1', 'HLA-DRB1', 'LUM', 'MOXD1', 'PRELP', 'RNASE1')

StromalGns_up <- StromalGns_up[StromalGns_up %in% rownames(usedTrainMat)]

#######
## For each gene, get the genes negatively correleted with it:

# init an empty list
mylist <- list()

CorList <- lapply(StromalGns_up, function(x){
  cr <- apply(usedTrainMat, 1, FUN = cor, usedTrainMat[x, ])
  cr_ord <- order(cr, decreasing = F)
  neg_cr <- cr_ord[1:500]
  neg_cr
  tmp <- c(rownames(usedTrainMat)[neg_cr])
  mylist[[x]] <- tmp
})

names(CorList) <- StromalGns_up

CorMat <- do.call("rbind", CorList)

CorMat <- CorMat %>%
  as.data.frame() %>%
  rownames_to_column(var = 'stromal_up') %>%
  pivot_longer(values_to = 'stromal_down', cols = -stromal_up) %>%
  select(-name)
      

StromalPairs <- as.matrix(CorMat)
StromalPairs2 <- expand_grid(CorMat[, 1], CorMat[, 2])
colnames(StromalPairs) <- c("proStroma","antiStroma")
colnames(StromalPairs2) <- c("proStroma","antiStroma")

save(StromalPairs, file = "./objs/StromaPairs.rda")
save(StromalPairs2, file = "./objs/StromaPairs2.rda")










