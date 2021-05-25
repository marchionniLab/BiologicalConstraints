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
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('MYC', names(x), value=TRUE))
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
ProMYC <- GeneList[c(#"LEE_LIVER_CANCER_MYC_UP",
                     #"LEE_LIVER_CANCER_MYC_TGFA_UP",
                     #"LEE_LIVER_CANCER_MYC_E2F1_UP",
                     #"COLLER_MYC_TARGETS_UP",
                     #"YU_MYC_TARGETS_UP",
                     "BILD_MYC_ONCOGENIC_SIGNATURE",
                     #"ELLWOOD_MYC_TARGETS_UP",
                     "PID_MYC_ACTIV_PATHWAY",
                     "HALLMARK_MYC_TARGETS_V2"
)]


#### Generate Consensus
allGns <- unique(unlist(ProMYC))
ProMYC <- sapply(ProMYC, function(x,y) y %in% x, y=allGns)
rownames(ProMYC) <- allGns

#### Summary
table(rowSums(ProMYC))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPM <- names(which(rowSums(ProMYC) > 0))
length(ConsensusPM)


######################################################
#### Anti Notch
AntiMYC <- GeneList[c(#"LEE_LIVER_CANCER_MYC_E2F1_DN", 
                      #"LEE_LIVER_CANCER_MYC_DN",
                      #"COLLER_MYC_TARGETS_DN",
                      #"YU_MYC_TARGETS_DN",
                      #"ELLWOOD_MYC_TARGETS_DN",
                      "PID_MYC_REPRESS_PATHWAY"
)]


#### Generate Consensus
allGns <- unique(unlist(AntiMYC))
AntiMYC <- sapply(AntiMYC, function(x,y) y %in% x, y=allGns)
rownames(AntiMYC) <- allGns

#### Summary
table(rowSums(AntiMYC))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAM <- names(which(rowSums(AntiMYC) > 0))
length(ConsensusAM)


#####################################################################
#### Combine
MYCList<- list(ConsensusPM, ConsensusAM)

#### Make Mutually Exclusive
allGns <- unique(unlist(MYCList))
MYCList<- sapply(MYCList, function(x,y) y %in% x, y=allGns)
rownames(MYCList) <- allGns

#### Summary
table(rowSums(MYCList))

#### Get Genes
MYCList <- MYCList[ rowSums(MYCList) == 1 , ]

#### Make as  List
MYCList <- apply(MYCList, 2, function(x) names(x[x]))


######################################################################
####Make Pairs
Pairs <- expand.grid(MYCList)

### Rename
colnames(Pairs) <- c("proMYC","antiMYC")

### Make into a matrix
MycPairs <- as.matrix(Pairs)

MycPairs[,1] <- gsub("\\-", "", MycPairs[,1])
MycPairs[,2] <- gsub("\\-", "", MycPairs[,2])

######################################################################
####Save Results
save(MYCList, file="./Objs/MycListConsensus.v6.1.rda")
save(MycPairs,file="./Objs/MycPairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
