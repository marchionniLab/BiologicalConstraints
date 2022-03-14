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
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('INFECTION', names(x), value=TRUE))
nms <- nms[sapply(nms,  length) > 0]
nms

#########################################################################
#### Select Lists and Extract Genes From All Lists
GeneList <- mapply(x=nms, y=geneSetCollectionSet.msigdb[names(nms)],  function(x, y) {
  geneIds(y[x])
})

### Remove a level
GeneList <- unlist(GeneList, recursive = FALSE)
names(GeneList) <- gsub(".+\\.",  "",  names(GeneList))

######################################################
#### Pro Viral Infection
ProVirInfection <- GeneList[c("AKL_HTLV1_INFECTION_UP",
                              "MARSHALL_VIRAL_INFECTION_RESPONSE_UP",
                              "DORN_ADENOVIRUS_INFECTION_48HR_UP"
                     
)]


#### Generate Consensus
allGns <- unique(unlist(ProVirInfection))
ProVirInfection <- sapply(ProVirInfection, function(x,y) y %in% x, y=allGns)
rownames(ProVirInfection) <- allGns

#### Summary
table(rowSums(ProVirInfection))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPI <- names(which(rowSums(ProVirInfection) > 0))
length(ConsensusPI)


######################################################
#### Anti viral infection
AntiVirInfection <- GeneList[c("AKL_HTLV1_INFECTION_DN",
                               "MARSHALL_VIRAL_INFECTION_RESPONSE_DN",
                               "DORN_ADENOVIRUS_INFECTION_48HR_DN"
                               
                               
)]


#### Generate Consensus
allGns <- unique(unlist(AntiVirInfection))
AntiVirInfection <- sapply(AntiVirInfection, function(x,y) y %in% x, y=allGns)
rownames(AntiVirInfection) <- allGns

#### Summary
table(rowSums(AntiVirInfection))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAI <- names(which(rowSums(AntiVirInfection) > 0))
length(ConsensusAI)


#####################################################################
#### Combine
InfectionList<- list(ConsensusPI, ConsensusAI)

#### Make Mutually Exclusive
allGns <- unique(unlist(InfectionList))
InfectionList<- sapply(InfectionList, function(x,y) y %in% x, y=allGns)
rownames(InfectionList) <- allGns

#### Summary
table(rowSums(InfectionList))

#### Get Genes
InfectionList <- InfectionList[ rowSums(InfectionList) == 1 , ]

#### Make as  List
InfectionList <- apply(InfectionList, 2, function(x) names(x[x]))


######################################################################
#### Make Pairs
VirInfectionPairs <- expand.grid(InfectionList)

##################################
## subsetting

# subset the TF-miRNA targets
# load("../../Genes/allTSPs.rda")
# filt <- which(VirInfectionPairs[,1] %in% myTSPs[,1] & VirInfectionPairs[,2] %in% myTSPs[,2])
# VirInfectionPairs <- VirInfectionPairs[-filt, ]

#####
# subset the NOTCH-MYC targets
# Notch = load("../Breast/Objs/NotchPairs.rda")
# MYC = load("../Breast/Objs/MycPairs.rda")
# myTSPs <- rbind(NotchPairs, MycPairs)
# filt <- which(VirInfectionPairs[,1] %in% myTSPs[,1] & VirInfectionPairs[,2] %in% myTSPs[,2])
# VirInfectionPairs <- VirInfectionPairs[-filt, ]
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
# filt <- which(VirInfectionPairs[,1] %in% myTSPs[,1] & VirInfectionPairs[,2] %in% myTSPs[,2])
# VirInfectionPairs <- VirInfectionPairs[-filt, ]
# 
# #####
# # subset the T-cell activation pairs
# load("./Objs/ImmunePairs.rda")
# filt <- which(VirInfectionPairs[,1] %in% ImmunePairs[,1] & VirInfectionPairs[,2] %in% ImmunePairs[,2])
# VirInfectionPairs <- VirInfectionPairs[-filt, ]
# 
### Rename
colnames(VirInfectionPairs) <- c("ProInfection","AntiInfection")
 
### Make into a matrix
VirInfectionPairs <- as.matrix(VirInfectionPairs)

write.csv(VirInfectionPairs, file = '/Users/mohamedomar/Dropbox (MechPred)/MechPred/USER/Mohamed/MechanisticModels/BiologicalConstraints/GenePairs/Infection_pairs.csv')

######################################################################
####Save Results
save(VirInfectionPairs,file="./Objs/VirInfectionPairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
