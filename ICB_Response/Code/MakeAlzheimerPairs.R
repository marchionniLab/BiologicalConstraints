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
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('ALZHEIMER', names(x), value=TRUE))
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
ProAlzheimer <- GeneList[c("BLALOCK_ALZHEIMERS_DISEASE_UP",
                          'WU_ALZHEIMER_DISEASE_UP'
                          
                          
)]


#### Generate Consensus
allGns <- unique(unlist(ProAlzheimer))
ProAlzheimer <- sapply(ProAlzheimer, function(x,y) y %in% x, y=allGns)
rownames(ProAlzheimer) <- allGns

#### Summary
table(rowSums(ProAlzheimer))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPI <- names(which(rowSums(ProAlzheimer) > 0))
length(ConsensusPI)


######################################################
#### Anti diabetes
AntiAlzheimer <- GeneList[c("BLALOCK_ALZHEIMERS_DISEASE_DN",
                           "WU_ALZHEIMER_DISEASE_DN"
)]


#### Generate Consensus
allGns <- unique(unlist(AntiAlzheimer))
AntiAlzheimer <- sapply(AntiAlzheimer, function(x,y) y %in% x, y=allGns)
rownames(AntiAlzheimer) <- allGns

#### Summary
table(rowSums(AntiAlzheimer))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAI <- names(which(rowSums(AntiAlzheimer) > 0))
length(ConsensusAI)


#####################################################################
#### Combine
AlzheimerList<- list(ConsensusPI, ConsensusAI)

#### Make Mutually Exclusive
allGns <- unique(unlist(AlzheimerList))
AlzheimerList<- sapply(AlzheimerList, function(x,y) y %in% x, y=allGns)
rownames(AlzheimerList) <- allGns

#### Summary
table(rowSums(AlzheimerList))

#### Get Genes
AlzheimerList <- AlzheimerList[ rowSums(AlzheimerList) == 1 , ]

#### Make as  List
AlzheimerList <- apply(AlzheimerList, 2, function(x) names(x[x]))


######################################################################
#### Make Pairs
AlzheimerPairs <- expand.grid(AlzheimerList)

### Rename
colnames(AlzheimerPairs) <- c("ProAlzheimer","AntiAlzheimer")

### Make into a matrix
AlzheimerPairs <- as.matrix(AlzheimerPairs)

######################################################################
####Save Results
save(AlzheimerPairs,file="./Objs/AlzheimerPairs.rda")

