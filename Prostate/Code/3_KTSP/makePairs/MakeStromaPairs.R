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
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('STROMA|FIBROBLAST|ENDOTHELIAL', names(x), value=TRUE))
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
#### ProLipogenesis
ProStroma <- GeneList[c("SUNG_METASTASIS_STROMA_UP",
                        #"CROONQUIST_STROMAL_STIMULATION_UP",
                        #"CROONQUIST_NRAS_VS_STROMAL_STIMULATION_DN",
                        "DURAND_STROMA_S_UP"
                        #"MISHRA_CARCINOMA_ASSOCIATED_FIBROBLAST_UP",
                        #"GO_POSITIVE_REGULATION_OF_FIBROBLAST_PROLIFERATION",
                        #"LU_TUMOR_ENDOTHELIAL_MARKERS_UP",
                        #"GO_POSITIVE_REGULATION_OF_ENDOTHELIAL_CELL_DIFFERENTIATION"
                              )]
                                          

#### Generate Consensus
allGns <- unique(unlist(ProStroma))
ProStroma <- sapply(ProStroma, function(x,y) y %in% x, y=allGns)
rownames(ProStroma) <- allGns

#### Summary
table(rowSums(ProStroma))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPS <- names(which(rowSums(ProStroma) > 0))
length(ConsensusPS)


######################################################
#### AntiAngiogenesis
AntiStroma <- GeneList[c("SUNG_METASTASIS_STROMA_DN",
                         #"CROONQUIST_STROMAL_STIMULATION_DN",
                         #"CROONQUIST_NRAS_VS_STROMAL_STIMULATION_UP",
                         "DURAND_STROMA_NS_UP"
                         #"MISHRA_CARCINOMA_ASSOCIATED_FIBROBLAST_DN",
                         #"GO_NEGATIVE_REGULATION_OF_FIBROBLAST_PROLIFERATION",
                         #"LU_TUMOR_ENDOTHELIAL_MARKERS_DN",
                         #"GO_NEGATIVE_REGULATION_OF_ENDOTHELIAL_CELL_PROLIFERATION"
                       )]


#### Generate Consensus
allGns <- unique(unlist(AntiStroma))
AntiStroma <- sapply(AntiStroma, function(x,y) y %in% x, y=allGns)
rownames(AntiStroma) <- allGns

#### Summary
table(rowSums(AntiStroma))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAS <- names(which(rowSums(AntiStroma) > 0))
length(ConsensusAS)


#####################################################################
#### Combine
StromaList<- list(ConsensusPS, ConsensusAS)

#### Make Mutually Exclusive
allGns <- unique(unlist(StromaList))
StromaList<- sapply(StromaList, function(x,y) y %in% x, y=allGns)
rownames(StromaList) <- allGns

#### Summary
table(rowSums(StromaList))

#### Get Genes
StromaList <- StromaList[ rowSums(StromaList) == 1 , ]

#### Make as  List
StromaList <- apply(StromaList, 2, function(x) names(x[x]))


######################################################################
####Make Pairs
Pairs <- expand.grid(StromaList)

### Rename
colnames(Pairs) <- c("proStroma","antiStroma")

### Make into a matrix
StromaPairs <- as.matrix(Pairs)

StromaPairs[,1] <- gsub("\\-", "", StromaPairs[,1])
StromaPairs[,2] <- gsub("\\-", "", StromaPairs[,2])

######################################################################
####Save Results
save(StromaList, file="./Objs/StromaListConsensus.v6.1.rda")
save(StromaPairs,file="objs/StromaPairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
