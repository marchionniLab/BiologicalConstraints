
## Venn diagrams of the mechanistic gene sets
# 1. TF_MiRNA genes.
# 2. Immune surveillance genes.
# 3. MiRNA genes.
# 4. Encode genes.
# 5. Pathway Commons genes.

#######################################################
rm(list = ls())

library(VennDiagram)
library(RColorBrewer)

## Load the Gene sets
load("./Objs/KTSP/KeepGns_TF_MiR.rda")
load("./Objs/KTSP/KeepGns_MiRNA.rda")
load("./Objs/KTSP/KeepGns_ImmuneSurv.rda")
load("./Objs/KTSP/KeepGns_ENC.rda")
load("./Objs/KTSP/KeepGns_PC.rda")


# Creat a list of different gene sets
MechGeneList <- list()
MechGeneList$TF_MiRNA <- keepGns_TF_MiR
MechGeneList$ImmuneSurv <- keepGns_ImmuneSurv
MechGeneList$MiRNA <- keepGns_MiRNA
MechGeneList$ENC <- keepGns_ENC
MechGeneList$PC <- keepGns_PC

## The venn diagram
# Prepare a palette of 5 colors with R colorbrewer:
myCol <- brewer.pal(5, "Pastel2")

venn.diagram(
  x = MechGeneList,
  filename = 'venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 2000 , 
  width = 2000 , 
  resolution = 300,

  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .5,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.8,
  cat.fontface = "bold",
  cat.dist = c(0.055, 0.055, 0.055, 0.055, 0.055),
  cat.fontfamily = "sans", 
  margin = 0.05, 
  main = "Different gene sets used in the mechanistic K-TSPs classifiers"
  )


