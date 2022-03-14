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
#### Pro alzheimer
ProAlzheimer <- GeneList[c('WU_ALZHEIMER_DISEASE_UP'

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
#### Anti alzheimer
AntiAlzheimer <- GeneList[c(
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

##################################
## subsetting

# subset the TF-miRNA targets
# load("../../Genes/allTSPs.rda")
# filt <- which(AlzheimerPairs[,1] %in% myTSPs[,1] & AlzheimerPairs[,2] %in% myTSPs[,2])
# AlzheimerPairs <- AlzheimerPairs[-filt, ]
# 
# #####
# # subset the NOTCH-MYC targets
# Notch = load("../Breast/Objs/NotchPairs.rda")
# MYC = load("../Breast/Objs/MycPairs.rda")
# myTSPs <- rbind(NotchPairs, MycPairs)
# filt <- which(AlzheimerPairs[,1] %in% myTSPs[,1] & AlzheimerPairs[,2] %in% myTSPs[,2])
# AlzheimerPairs <- AlzheimerPairs[-filt, ]
# 
# ######
# # subset the adhesion, activation and O2 response pairs
# Genes1 <- read.delim("../Prostate/objs/GO_Adhesion.txt")
# Genes1 <- as.matrix(Genes1)
# Genes1 <- Genes1[-1,]
# 
# Genes2 <- read.delim("../Prostate/objs/GO_Activation.txt")
# Genes2 <- as.matrix(Genes2)
# Genes2 <- Genes2[-1,]
# 
# Genes3 <- read.delim("../Prostate/objs/GO_O2Response.txt")
# Genes3 <- as.matrix(Genes3)
# Genes3 <- Genes3[-1,]
# 
# Genes <- c(Genes1,Genes2, Genes3)
# Genes <- Genes[!duplicated(Genes)]
# 
# myTSPs <- t(combn(Genes,2))
# 
# filt <- which(AlzheimerPairs[,1] %in% myTSPs[,1] & AlzheimerPairs[,2] %in% myTSPs[,2])
# AlzheimerPairs <- AlzheimerPairs[-filt, ]

#####
# subset the T-cell activation pairs
# load("./Objs/ImmunePairs.rda")
# filt <- which(AlzheimerPairs[,1] %in% ImmunePairs[,1] & AlzheimerPairs[,2] %in% ImmunePairs[,2])
# AlzheimerPairs <- AlzheimerPairs[-filt, ]

###############################################
### Rename
colnames(AlzheimerPairs) <- c("ProAlzheimer","AntiAlzheimer")

### Make into a matrix
AlzheimerPairs <- as.matrix(AlzheimerPairs)

write.csv(AlzheimerPairs, file = '/Users/mohamedomar/Dropbox (MechPred)/MechPred/USER/Mohamed/MechanisticModels/BiologicalConstraints/GenePairs/Alzheimer_pairs.csv')
######################################################################
####Save Results
save(AlzheimerPairs,file="./Objs/AlzheimerPairs.rda")

