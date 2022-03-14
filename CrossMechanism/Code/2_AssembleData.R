
rm(list = ls())
gc()

## Assemble the datasets

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

list.data = list()

## -------------
## Bladder

l1 = load("../Bladder/Objs/progressionDataGood2.rda")
l2 = load("../Bladder/Objs/Correlation/RGenes.rda")

usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

list.data[["Bladder"]] = list(
  trainMat = usedTrainMat,
  testMat = usedTestMat,
  trainGroup = usedTrainGroup,
  testGroup = usedTestGroup
)

rm(list=c(l1, l2))

## -------------
## Breast

l3 = load("../Breast/Objs/ChemoDataNew.rda")

usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

list.data[["Breast"]] = list(
  trainMat = usedTrainMat,
  testMat = usedTestMat,
  trainGroup = usedTrainGroup,
  testGroup = usedTestGroup
)

rm(list=l3)

## -------------
## Prostate

l4 = load("../Prostate/Objs/MetastasisDataGood.rda")
l5 = load("../Prostate/Objs/Correlation/RGenes.rda")

usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

list.data[["Prostate"]] = list(
  trainMat = usedTrainMat,
  testMat = usedTestMat,
  trainGroup = usedTrainGroup,
  testGroup = usedTestGroup
)

rm(list = c(l4,l5))

## -------------

save(list.data, file="../CrossMechanism/Objs/list.data.rda")

