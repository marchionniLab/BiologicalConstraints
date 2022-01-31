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
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('T_CELL', names(x), value=TRUE))
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
ProImmune <- GeneList[c('POSITIVE_REGULATION_OF_T_CELL_ACTIVATION',
                        'GO_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE',
                        'GO_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_CYTOTOXICITY',
                        'GO_POSITIVE_REGULATION_OF_T_CELL_CYTOKINE_PRODUCTION',
                        'GO_POSITIVE_REGULATION_OF_ACTIVATED_T_CELL_PROLIFERATION',
                        'GO_POSITIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNITY',
                        'GO_POSITIVE_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION',
                        'GO_POSITIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION',
                        #######
                        "GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_UP", 
                        "GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_UP",
                        "GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_DN"
            
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
AntiImmune <- GeneList[c("GO_NEGATIVE_REGULATION_OF_T_CELL_RECEPTOR_SIGNALING_PATHWAY", 
                        "GO_REGULATION_OF_T_CELL_CYTOKINE_PRODUCTION",
                        "GO_NEGATIVE_REGULATION_OF_T_CELL_MEDIATED_IMMUNITY",
                        'GO_NEGATIVE_REGULATION_OF_ALPHA_BETA_T_CELL_ACTIVATION',
                        'GO_NEGATIVE_REGULATION_OF_CD4_POSITIVE_ALPHA_BETA_T_CELL_ACTIVATION',
                        #####
                        'GSE9650_NAIVE_VS_EXHAUSTED_CD8_TCELL_DN',
                        'GSE9650_EFFECTOR_VS_EXHAUSTED_CD8_TCELL_DN',
                        'GSE9650_EXHAUSTED_VS_MEMORY_CD8_TCELL_UP'
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
ImmunePairs <- as.matrix(Pairs)

# ImmunePairs[,1] <- gsub("\\-", "", ImmunePairs[,1])
# ImmunePairs[,2] <- gsub("\\-", "", ImmunePairs[,2])

######################################################################
####Save Results
save(ImmuneList, file="./Objs/ImmuneListConsensus.v6.1.rda")
save(ImmunePairs,file="./Objs/ImmunePairs2.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
