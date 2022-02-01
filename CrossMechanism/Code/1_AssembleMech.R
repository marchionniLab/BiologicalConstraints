
rm(list = ls())
gc()

## Assemble the different mechanisms

############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(enrichR)
library(mltools)
library(xtable)
library(boot)
library(patchwork)

## -------------

list.mech = list()

## -------------
## TF-MiR (Bladder) mechanism

l1 = load("../../../Genes/allTSPs.rda")

list.mech[["TF_MIR"]] = myTSPs

rm(myTSPs)

## -------------
## NOTCH and MYC (Breast) mechanism

l2 = load("../../Breast/Objs/NotchPairs.rda")
l3 = load("../../Breast/Objs/MycPairs.rda")

myTSPs <- rbind(NotchPairs, MycPairs)
colnames(myTSPs) <- c("BadGene", "GoodGene")

list.mech[["NOTCH_MYC"]] = myTSPs

rm(myTSPs)

## -------------
## GO (Prostate) mechanism


Genes1 <- read.delim("../../Prostate/objs/GO_Adhesion.txt")
Genes1 <- as.matrix(Genes1)
Genes1 <- Genes1[-1,]

Genes2 <- read.delim("../../Prostate/objs/GO_Activation.txt")
Genes2 <- as.matrix(Genes2)
Genes2 <- Genes2[-1,]

Genes3 <- read.delim("../../Prostate/objs/GO_O2Response.txt")
Genes3 <- as.matrix(Genes3)
Genes3 <- Genes3[-1,]

Genes <- c(Genes1,Genes2, Genes3)
Genes <- Genes[!duplicated(Genes)]

myTSPs <- t(combn(Genes,2))

list.mech[["GO"]] = myTSPs

## -------------
## Alzheimer mechanism

l4 = load("../Objs/CancerUnrelatedMechanisms/AlzheimerPairs.rda")

list.mech[["Alzheimer"]] = AlzheimerPairs

rm(list=l4)

## -------------
## Diabetes mechanism

l5 = load("../Objs/CancerUnrelatedMechanisms/DiabetesPairs.rda")

list.mech[["Diabetes"]] = DiabetesPairs

rm(list=l5)

## -------------
## Viral infection mechanism

l6 = load("../Objs/CancerUnrelatedMechanisms/VirInfectionPairs.rda")

list.mech[["Viral"]] = VirInfectionPairs

rm(list=l6)

## -------------

save(list.mech, file="../Objs/list.mech.rda")





