######################################################################
#### Clear List and set options
rm(list = ls())


#########################################################################
#### Load 
load("./Data/geneSetCollectionSet.MSigDBv6.1.rda")


#########################################################################
#### Package
require(GSEABase)


#########################################################################
#### Retrieve the Lists Needed
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('NOTCH', names(x), value=TRUE))
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
#### Pro Notch
ProNotch <- GeneList[c("NGUYEN_NOTCH1_TARGETS_UP",
                      "VILIMAS_NOTCH1_TARGETS_UP",
                      "REACTOME_PRE_NOTCH_TRANSCRIPTION_AND_TRANSLATION",
                      "REACTOME_ACTIVATED_NOTCH1_TRANSMITS_SIGNAL_TO_THE_NUCLEUS",
                      "REACTOME_SIGNALING_BY_NOTCH4",
                      "REACTOME_SIGNALING_BY_NOTCH2",
                      "REACTOME_PRE_NOTCH_EXPRESSION_AND_PROCESSING",
                      "REACTOME_NOTCH1_INTRACELLULAR_DOMAIN_REGULATES_TRANSCRIPTION",
                      "REACTOME_PRE_NOTCH_PROCESSING_IN_GOLGI",
                      "REACTOME_SIGNALING_BY_NOTCH1",
                      "REACTOME_SIGNALING_BY_NOTCH3",
                      "REACTOME_SIGNALING_BY_NOTCH",
                      "GO_POSITIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY"
)]
                                      

#### Generate Consensus
allGns <- unique(unlist(ProNotch))
ProNotch <- sapply(ProNotch, function(x,y) y %in% x, y=allGns)
rownames(ProNotch) <- allGns

#### Summary
table(rowSums(ProNotch))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPN <- names(which(rowSums(ProNotch) > 0))
length(ConsensusPN)


######################################################
#### Anti Notch
AntiNotch <- GeneList[c("NGUYEN_NOTCH1_TARGETS_DN", 
                       "VILIMAS_NOTCH1_TARGETS_DN",
                       "GO_NEGATIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY"
)]


#### Generate Consensus
allGns <- unique(unlist(AntiNotch))
AntiNotch <- sapply(AntiNotch, function(x,y) y %in% x, y=allGns)
rownames(AntiNotch) <- allGns

#### Summary
table(rowSums(AntiNotch))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAN <- names(which(rowSums(AntiNotch) > 0))
length(ConsensusAN)


#####################################################################
#### Combine
NotchList<- list(ConsensusPN, ConsensusAN)

#### Make Mutually Exclusive
allGns <- unique(unlist(NotchList))
NotchList<- sapply(NotchList, function(x,y) y %in% x, y=allGns)
rownames(NotchList) <- allGns

#### Summary
table(rowSums(NotchList))

#### Get Genes
NotchList <- NotchList[ rowSums(NotchList) == 1 , ]

#### Make as  List
NotchList <- apply(NotchList, 2, function(x) names(x[x]))


######################################################################
####Make Pairs
Pairs <- expand.grid(NotchList)

### Rename
colnames(Pairs) <- c("proNotch","antiNotch")

### Make into a matrix
NotchPairs <- as.matrix(Pairs)

NotchPairs[,1] <- gsub("\\-", "", NotchPairs[,1])
NotchPairs[,2] <- gsub("\\-", "", NotchPairs[,2])

######################################################################
####Save Results
save(NotchList, file="./Objs/NotchListConsensus.v6.1.rda")
save(NotchPairs,file="./Objs/NotchPairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
