######################################################################
#### Clear List and set options
rm(list = ls())


#########################################################################
#### Load 
load("./Data/geneSetCollectionSet.MSigDBv6.1.rda")


#########################################################################
#### Package
library(GSEABase)


#########################################################################
#### Retrieve the Lists Needed
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('DIABETES', names(x), value=TRUE))
nms <- nms[sapply(nms,  length) > 0]


#########################################################################
#### Select Lists and Extract Genes From All Lists
GeneList <- mapply(x=nms, y=geneSetCollectionSet.msigdb[names(nms)],  function(x, y) {
  geneIds(y[x])
})

### Remove a level
GeneList <- unlist(GeneList, recursive = FALSE)
names(GeneList) <- gsub(".+\\.",  "",  names(GeneList))

######################################################
#### Pro Diabetes
ProDiabetes <- GeneList[c("GSE9006_TYPE_1_DIABETES_AT_DX_VS_1MONTH_POST_DX_PBMC_UP",
                        'GSE9006_TYPE_1_DIABETES_AT_DX_VS_4MONTH_POST_DX_PBMC_UP'
                        
                        
)]


#### Generate Consensus
allGns <- unique(unlist(ProDiabetes))
ProDiabetes <- sapply(ProDiabetes, function(x,y) y %in% x, y=allGns)
rownames(ProDiabetes) <- allGns

#### Summary
table(rowSums(ProDiabetes))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPI <- names(which(rowSums(ProDiabetes) > 0))
length(ConsensusPI)


######################################################
#### Anti diabetes
AntiDiabetes <- GeneList[c("GSE9006_TYPE_1_DIABETES_AT_DX_VS_1MONTH_POST_DX_PBMC_DN",
                         "GSE9006_TYPE_1_DIABETES_AT_DX_VS_4MONTH_POST_DX_PBMC_DN"
              )]


#### Generate Consensus
allGns <- unique(unlist(AntiDiabetes))
AntiDiabetes <- sapply(AntiDiabetes, function(x,y) y %in% x, y=allGns)
rownames(AntiDiabetes) <- allGns

#### Summary
table(rowSums(AntiDiabetes))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAI <- names(which(rowSums(AntiDiabetes) > 0))
length(ConsensusAI)


#####################################################################
#### Combine
DiabetesList<- list(ConsensusPI, ConsensusAI)

#### Make Mutually Exclusive
allGns <- unique(unlist(DiabetesList))
DiabetesList<- sapply(DiabetesList, function(x,y) y %in% x, y=allGns)
rownames(DiabetesList) <- allGns

#### Summary
table(rowSums(DiabetesList))

#### Get Genes
DiabetesList <- DiabetesList[ rowSums(DiabetesList) == 1 , ]

#### Make as  List
DiabetesList <- apply(DiabetesList, 2, function(x) names(x[x]))


######################################################################
#### Make Pairs
DiabetesPairs <- expand.grid(DiabetesList)

### Rename
colnames(DiabetesPairs) <- c("ProDiabetes","AntiDiabetes")

### Make into a matrix
DiabetesPairs <- as.matrix(DiabetesPairs)

######################################################################
####Save Results
save(DiabetesPairs,file="./Objs/DiabetesPairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
