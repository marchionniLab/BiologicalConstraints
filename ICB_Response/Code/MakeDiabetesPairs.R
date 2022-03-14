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
ConsensusPI <- names(which(rowSums(ProDiabetes) > 1))
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
ConsensusAI <- names(which(rowSums(AntiDiabetes) > 1))
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

##################################
## subsetting

# subset the TF-miRNA targets
#load("../../Genes/allTSPs.rda")
#filt <- which(DiabetesPairs[,1] %in% myTSPs[,1] & DiabetesPairs[,2] %in% myTSPs[,2])
#DiabetesPairs <- DiabetesPairs[-filt, ]

# #####
# # subset the NOTCH-MYC targets
# Notch = load("../Breast/Objs/NotchPairs.rda")
# MYC = load("../Breast/Objs/MycPairs.rda")
# myTSPs <- rbind(NotchPairs, MycPairs)
# filt <- which(DiabetesPairs[,1] %in% myTSPs[,1] & DiabetesPairs[,2] %in% myTSPs[,2])
# DiabetesPairs <- DiabetesPairs[-filt, ]
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
# filt <- which(DiabetesPairs[,1] %in% myTSPs[,1] & DiabetesPairs[,2] %in% myTSPs[,2])
# DiabetesPairs <- DiabetesPairs[-filt, ]

#####
# subset the T-cell activation pairs
# load("./Objs/ImmunePairs.rda")
# filt <- which(DiabetesPairs[,1] %in% ImmunePairs[,1] & DiabetesPairs[,2] %in% ImmunePairs[,2])
# DiabetesPairs <- DiabetesPairs[-filt, ]

### Rename
colnames(DiabetesPairs) <- c("ProDiabetes","AntiDiabetes")

### Make into a matrix
DiabetesPairs <- as.matrix(DiabetesPairs)

write.csv(DiabetesPairs, file = '/Users/mohamedomar/Dropbox (MechPred)/MechPred/USER/Mohamed/MechanisticModels/BiologicalConstraints/GenePairs/Diabetes_pairs.csv')

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
