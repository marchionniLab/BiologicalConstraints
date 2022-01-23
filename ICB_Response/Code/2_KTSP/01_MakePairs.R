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
ProImmune <- GeneList[c(
)]


#### Generate Consensus
allGns <- unique(unlist(ProImmune))
ProImmune <- sapply(ProImmune, function(x,y) y %in% x, y=allGns)
rownames(ProImmune) <- allGns

#### Summary
table(rowSums(ProImmune))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPI <- names(which(rowSums(ProImmune) > 0))
length(ConsensusPI)


######################################################
#### Anti Notch
AntiImmune <- GeneList[c("NGUYEN_NOTCH1_TARGETS_DN", 
                        "VILIMAS_NOTCH1_TARGETS_DN",
                        "GO_NEGATIVE_REGULATION_OF_NOTCH_SIGNALING_PATHWAY"
)]


#### Generate Consensus
allGns <- unique(unlist(AntiImmune))
AntiImmune <- sapply(AntiImmune, function(x,y) y %in% x, y=allGns)
rownames(AntiImmune) <- allGns

#### Summary
table(rowSums(AntiImmune))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAI <- names(which(rowSums(AntiImmune) > 0))
length(ConsensusAI)


#####################################################################
#### Combine
ImmuneList<- list(ConsensusPI, ConsensusAI)

#### Make Mutually Exclusive
allGns <- unique(unlist(ImmuneList))
ImmuneList<- sapply(ImmuneList, function(x,y) y %in% x, y=allGns)
rownames(ImmuneList) <- allGns

#### Summary
table(rowSums(ImmuneList))

#### Get Genes
ImmuneList <- ImmuneList[ rowSums(ImmuneList) == 1 , ]

#### Make as  List
ImmuneList <- apply(ImmuneList, 2, function(x) names(x[x]))


######################################################################
####Make Pairs
Pairs <- expand.grid(ImmuneList)

### Rename
colnames(Pairs) <- c("ProImmune","AntiImmune")

### Make into a matrix
NotchPairs <- as.matrix(Pairs)

NotchPairs[,1] <- gsub("\\-", "", NotchPairs[,1])
NotchPairs[,2] <- gsub("\\-", "", NotchPairs[,2])

######################################################################
####Save Results
save(ImmuneList, file="./Objs/ImmuneListConsensus.v6.1.rda")
save(NotchPairs,file="./Objs/ImmunePairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
