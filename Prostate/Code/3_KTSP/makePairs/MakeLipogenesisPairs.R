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
nms <- lapply(geneSetCollectionSet.msigdb, function(x) grep('FATTY|LIPID', names(x), value=TRUE))
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
ProLipogenesis <- GeneList[c("KEGG_BIOSYNTHESIS_OF_UNSATURATED_FATTY_ACIDS",  
                    "KEGG_FATTY_ACID_METABOLISM",
                    "REACTOME_FATTY_ACYL_COA_BIOSYNTHESIS",
                    "REACTOME_FATTY_ACID_TRIACYLGLYCEROL_AND_KETONE_BODY_METABOLISM",  
                    "REACTOME_SYNTHESIS_OF_VERY_LONG_CHAIN_FATTY_ACYL_COAS",
                    "GO_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS",
                    "GO_REGULATION_OF_FATTY_ACID_METABOLIC_PROCESS",
                    "GO_POSITIVE_REGULATION_OF_FATTY_ACID_BIOSYNTHETIC_PROCESS",
                    "GO_FATTY_ACYL_COA_METABOLIC_PROCESS",
                    #"GO_FATTY_ACID_HOMEOSTASIS",
                    "GO_SHORT_CHAIN_FATTY_ACID_METABOLIC_PROCESS",
                    "GO_VERY_LONG_CHAIN_FATTY_ACID_METABOLIC_PROCESS",
                    "GO_FATTY_ACID_BIOSYNTHETIC_PROCESS",
                    "GO_FATTY_ACID_METABOLIC_PROCESS",
                    "GO_REGULATION_OF_FATTY_ACID_BIOSYNTHETIC_PROCESS",  
                    "GO_POSITIVE_REGULATION_OF_FATTY_ACID_METABOLIC_PROCESS",
                    "GO_UNSATURATED_FATTY_ACID_METABOLIC_PROCESS", 
                    "GO_FATTY_ACID_DERIVATIVE_BIOSYNTHETIC_PROCESS",
                    "GO_UNSATURATED_FATTY_ACID_BIOSYNTHETIC_PROCESS", 
                    "GO_FATTY_ACID_ELONGATION", 
                    #"GO_LONG_CHAIN_FATTY_ACID_BINDING", 
                    #"GO_FATTY_ACYL_COA_BINDING", 
                    #"GO_FATTY_ACID_BINDING", 
                    #"GO_LONG_CHAIN_FATTY_ACID_COA_LIGASE_ACTIVITY", 
                    #"GO_FATTY_ACID_LIGASE_ACTIVITY",
                    #"GO_REGULATION_OF_FATTY_ACID_TRANSPORT",
                    #"GO_FATTY_ACID_DERIVATIVE_TRANSPORT",
                    #"GO_LONG_CHAIN_FATTY_ACID_TRANSPORT",
                    #"GO_FATTY_ACID_TRANSPORT",
                    "HALLMARK_FATTY_ACID_METABOLISM", 
                    "FATTY_ACID_METABOLIC_PROCESS", 
                    "KEGG_GLYCEROLIPID_METABOLISM", 
                    "KEGG_GLYCEROPHOSPHOLIPID_METABOLISM", 
                    "KEGG_ETHER_LIPID_METABOLISM", 
                    "KEGG_SPHINGOLIPID_METABOLISM", 
                    "KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_LACTO_AND_NEOLACTO_SERIES", 
                    "KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_GLOBO_SERIES", 
                    "KEGG_GLYCOSPHINGOLIPID_BIOSYNTHESIS_GANGLIO_SERIES", 
                    "REACTOME_SPHINGOLIPID_DE_NOVO_BIOSYNTHESIS", 
                    "REACTOME_GLYCOSPHINGOLIPID_METABOLISM", 
                    "REACTOME_PHOSPHOLIPID_METABOLISM", 
                    "REACTOME_GLYCEROPHOSPHOLIPID_BIOSYNTHESIS", 
                    #"REACTOME_HDL_MEDIATED_LIPID_TRANSPORT", 
                    "REACTOME_PEROXISOMAL_LIPID_METABOLISM", 
                    "REACTOME_SPHINGOLIPID_METABOLISM", 
                    "REACTOME_METABOLISM_OF_LIPIDS_AND_LIPOPROTEINS", 
                    #"REACTOME_LIPID_DIGESTION_MOBILIZATION_AND_TRANSPORT", 
                    #"REACTOME_CHYLOMICRON_MEDIATED_LIPID_TRANSPORT", 
                    #"GO_LIPID_MODIFICATION", 
                    #"GO_REGULATION_OF_LIPID_STORAGE", 
                    #"GO_CELLULAR_RESPONSE_TO_LIPID", 
                    "GO_NEGATIVE_REGULATION_OF_LIPID_METABOLIC_PROCESS", 
                    #"GO_POSITIVE_REGULATION_OF_LIPID_KINASE_ACTIVITY", 
                    #"GO_SPHINGOLIPID_MEDIATED_SIGNALING_PATHWAY",
                    #"GO_PROTEIN_LIPID_COMPLEX_ASSEMBLY", 
                    "GO_REGULATION_OF_MEMBRANE_LIPID_METABOLIC_PROCESS", 
                    #"GO_REGULATION_OF_MEMBRANE_LIPID_DISTRIBUTION", 
                    #"GO_POSITIVE_REGULATION_OF_LIPID_TRANSPORT", 
                    "GO_REGULATION_OF_PHOSPHOLIPID_METABOLIC_PROCESS", 
                    "GO_PHOSPHOLIPID_BIOSYNTHETIC_PROCESS", 
                    #"GO_LIPID_HOMEOSTASIS", 
                    "GO_GLYCEROLIPID_METABOLIC_PROCESS", 
                    "GO_REGULATION_OF_LIPID_METABOLIC_PROCESS", 
                    "GO_SPHINGOLIPID_METABOLIC_PROCESS", 
                    "GO_POSITIVE_REGULATION_OF_PHOSPHOLIPID_METABOLIC_PROCESS", 
                    #"GO_REGULATION_OF_LIPID_KINASE_ACTIVITY", 
                    #"GO_PROTEIN_LIPID_COMPLEX_SUBUNIT_ORGANIZATION", 
                    #"GO_LIPID_TRANSLOCATION", 
                    #"GO_REGULATION_OF_LIPID_TRANSPORT", 
                    "GO_GLYCEROLIPID_BIOSYNTHETIC_PROCESS", 
                    "GO_GLYCOSPHINGOLIPID_METABOLIC_PROCESS", 
                    "GO_GLYCEROPHOSPHOLIPID_METABOLIC_PROCESS", 
                    #"GO_RESPONSE_TO_LIPID", 
                    "GO_NEUTRAL_LIPID_BIOSYNTHETIC_PROCESS", 
                    "GO_LIPID_METABOLIC_PROCESS", 
                    "GO_NEUTRAL_LIPID_METABOLIC_PROCESS", 
                    #"GO_CELLULAR_RESPONSE_TO_FATTY_ACID", 
                    "GO_LIPID_BIOSYNTHETIC_PROCESS", 
                    "GO_REGULATION_OF_PHOSPHOLIPID_BIOSYNTHETIC_PROCESS", 
                    "GO_SPHINGOLIPID_BIOSYNTHETIC_PROCESS", 
                    "GO_NEGATIVE_REGULATION_OF_LIPID_CATABOLIC_PROCESS", 
                    "GO_POSITIVE_REGULATION_OF_LIPID_METABOLIC_PROCESS", 
                    #"GO_INTRACELLULAR_LIPID_TRANSPORT", 
                    "GO_POSITIVE_REGULATION_OF_LIPID_BIOSYNTHETIC_PROCESS",
                    #"GO_LIPID_STORAGE", 
                    "GO_MEMBRANE_LIPID_METABOLIC_PROCESS", 
                    "GO_UNSATURATED_FATTY_ACID_METABOLIC_PROCESS", 
                    #"GO_POSITIVE_REGULATION_OF_LIPID_STORAGE", 
                    "GO_MEMBRANE_LIPID_BIOSYNTHETIC_PROCESS", 
                    "GO_OLIGOSACCHARIDE_LIPID_INTERMEDIATE_BIOSYNTHETIC_PROCESS")]
                    #"GO_PROTEIN_LIPID_COMPLEX_BINDING", 
                    #"GO_GLYCOLIPID_BINDING", 
                    #"GO_BIOACTIVE_LIPID_RECEPTOR_ACTIVITY", 
                    #"GO_LIPID_BINDING")]                                             

#### Generate Consensus
allGns <- unique(unlist(ProLipogenesis))
ProLipogenesis <- sapply(ProLipogenesis, function(x,y) y %in% x, y=allGns)
rownames(ProLipogenesis) <- allGns

#### Summary
table(rowSums(ProLipogenesis))

#### Get the Consensus Genes from ProAngiogenesis
ConsensusPL <- names(which(rowSums(ProLipogenesis) > 0))
length(ConsensusPL)


######################################################
#### AntiAngiogenesis
AntiLipogenesis <- GeneList[c(#"FATTY_ACID_BETA_OXIDATION", 
                               #"FATTY_ACID_OXIDATION",
                     #"REACTOME_ACTIVATED_AMPK_STIMULATES_FATTY_ACID_OXIDATION_IN_MUSCLE",
                     #"REACTOME_MITOCHONDRIAL_FATTY_ACID_BETA_OXIDATION",
                     "GO_FATTY_ACID_CATABOLIC_PROCESS",
                     #"GO_FATTY_ACID_BETA_OXIDATION",
                     #"GO_FATTY_ACID_BETA_OXIDATION_USING_ACYL_COA_OXIDASE",
                     #"GO_REGULATION_OF_FATTY_ACID_OXIDATION",
                     #"GO_POSITIVE_REGULATION_OF_FATTY_ACID_OXIDATION",
                     #"GO_REGULATION_OF_FATTY_ACID_BETA_OXIDATION", 
                     "LIPID_CATABOLIC_PROCESS", 
                     "CELLULAR_LIPID_CATABOLIC_PROCESS", 
                     "GO_NEGATIVE_REGULATION_OF_LIPID_METABOLIC_PROCESS", 
                     "GO_MEMBRANE_LIPID_CATABOLIC_PROCESS", 
                     "GO_GLYCOLIPID_CATABOLIC_PROCESS", 
                     "GO_NEUTRAL_LIPID_CATABOLIC_PROCESS",
                     "GO_CELLULAR_LIPID_CATABOLIC_PROCESS", 
                     "GO_POSITIVE_REGULATION_OF_LIPID_CATABOLIC_PROCESS", 
                     "GO_NEGATIVE_REGULATION_OF_LIPID_STORAGE", 
                     "GO_LIPID_CATABOLIC_PROCESS", 
                     "GO_GLYCEROLIPID_CATABOLIC_PROCESS", 
                     "GO_REGULATION_OF_LIPID_CATABOLIC_PROCESS")]


#### Generate Consensus
allGns <- unique(unlist(AntiLipogenesis))
AntiLipogenesis <- sapply(AntiLipogenesis, function(x,y) y %in% x, y=allGns)
rownames(AntiLipogenesis) <- allGns

#### Summary
table(rowSums(AntiLipogenesis))

#### Get the Consensus Genes from AntiAgenesis
ConsensusAL <- names(which(rowSums(AntiLipogenesis) > 0))
length(ConsensusAL)


#####################################################################
#### Combine
lipoList<- list(ConsensusPL, ConsensusAL)

#### Make Mutually Exclusive
allGns <- unique(unlist(lipoList))
lipoList<- sapply(lipoList, function(x,y) y %in% x, y=allGns)
rownames(lipoList) <- allGns

#### Summary
table(rowSums(lipoList))

#### Get Genes
lipoList <- lipoList[ rowSums(lipoList) == 1 , ]

#### Make as  List
lipoList <- apply(lipoList, 2, function(x) names(x[x]))


######################################################################
####Make Pairs
Pairs <- expand.grid(lipoList)

### Rename
colnames(Pairs) <- c("proLipogenesis","antiLipogenesis")

### Make into a matrix
lipogenesisPairs <- as.matrix(Pairs)

lipogenesisPairs[,1] <- gsub("\\-", "", lipogenesisPairs[,1])
lipogenesisPairs[,2] <- gsub("\\-", "", lipogenesisPairs[,2])

######################################################################
####Save Results
save(lipoList, file="./Objs/LipoListConsensus.v6.1.rda")
save(lipogenesisPairs,file="objs/LipogenesisPairs.rda")


#####################################################################
### Session information
date()
sessionInfo()

#####################################################################
### Quitting in a clean way
rm(list=ls())
q("no")
