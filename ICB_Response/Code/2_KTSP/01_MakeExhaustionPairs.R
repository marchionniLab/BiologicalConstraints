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
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('EXHAUS', names(x), value=TRUE))
nms <- nms[sapply(nms,  length) > 0]

#########################################################################
#### Select Lists and Extract Genes From All Lists
GeneList <- mapply(x=nms, y=geneSetCollectionSet.msigdb[names(nms)],  function(x, y) {
  geneIds(y[x])
}, SIMPLIFY = F)

### Remove a level
GeneList <- unlist(GeneList, recursive = FALSE)
names(GeneList) <- gsub(".+\\.",  "",  names(GeneList))
######################################################
#### Pro Notch
ProExhaustion <- GeneList[c('GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN',
                        'GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN',
                        'GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP'
                        )]


#### Generate Consensus
allGns <- unique(unlist(ProExhaustion))
ProExhaustion <- sapply(ProExhaustion, function(x,y) y %in% x, y=allGns)
rownames(ProExhaustion) <- allGns

#### Summary
table(rowSums(ProExhaustion))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPI <- names(which(rowSums(ProExhaustion) > 0))
length(ConsensusPI)


######################################################
#### Anti Notch
AntiExhausion <- GeneList[c("GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP", 
                         "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
                         "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN"
                         )]


#### Generate Consensus
allGns <- unique(unlist(AntiExhausion))
AntiExhausion <- sapply(AntiExhausion, function(x,y) y %in% x, y=allGns)
rownames(AntiExhausion) <- allGns

#### Summary
table(rowSums(AntiExhausion))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAI <- names(which(rowSums(AntiExhausion) > 0))
length(ConsensusAI)


#####################################################################
#### Combine
ExhausionList<- list(ConsensusPI, ConsensusAI)

#### Make Mutually Exclusive
allGns <- unique(unlist(ExhausionList))
ExhausionList<- sapply(ExhausionList, function(x,y) y %in% x, y=allGns)
rownames(ExhausionList) <- allGns

#### Summary
table(rowSums(ExhausionList))

#### Get Genes
ExhausionList <- ExhausionList[ rowSums(ExhausionList) == 1 , ]

#### Make as  List
ExhausionList <- apply(ExhausionList, 2, function(x) names(x[x]))


######################################################################
####Make Pairs
Pairs <- expand.grid(ExhausionList)

### Rename
colnames(Pairs) <- c("ProExhaustion","AntiExhausion")

### Make into a matrix
ExhausionPairs <- as.matrix(Pairs)

# ImmunePairs[,1] <- gsub("\\-", "", ImmunePairs[,1])
# ImmunePairs[,2] <- gsub("\\-", "", ImmunePairs[,2])

######################################################################
####Save Results
save(ExhausionList, file="./Objs/ExhausionListConsensus.v6.1.rda")
save(ExhausionPairs,file="./Objs/ExhausionPairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
