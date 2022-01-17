###################################################################################
### Mohamed Omar
### April 1, 2020
### ### Goal: Creating the restricted ktsp classifier.
### invasion genes
#################################################################################

###### 
# Clean Work space
rm(list = ls())
# Set work directory
setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

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

###########################################################################
### Load expression and phenotype data
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

############################################################################
## Load the selected genes
Genes <- load("/Users/mohamedomar/Desktop/Genes/MECHANISM/CancerHallmarks/Invasion/objs/emtPairs.v6.1.rda")

myTSPs <- emtPairs

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(usedTrainMat))
#keepGns_TF_MiR <- keepGns
#save(keepGns_TF_MiR, file = "./Objs/KTSP/KeepGns_TF_MiR.rda")

#usedTrainMat <- usedTrainMat[keepGns, ]
#usedTestMat <- usedTestMat[keepGns, ]

### Associated groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

#print(xtable(myTSPs, type = "latex"), file = "./Objs/KTSP/Restricted_Pairs.tex")
#write.csv(myTSPs, file = "./Objs/KTSP/Restricted_Pairs.csv")

###########################################################################
### TRAINING using restricted pairs
###########################################################################

### Set Feature number and max k
ktsp <- c(3:25) # 23 
featNo <- nrow(usedTrainMat)

### Train a classifier using default filtering function based on Wilcoxon
set.seed(333)

ktspPredictorRes <- SWAP.Train.KTSP(
  usedTrainMat, usedTrainGroup, krange=ktsp,
  FilterFunc = SWAP.Filter.Wilcoxon, featureNo=featNo, RestrictedPairs = myTSPs)

ktspPredictorRes

#save(ktspPredictorRes, file = "./Objs/KTSP/MechKTSP94Pairs.rda")

#Mechanistic_KTSP <- cbind(ktspPredictorRes$TSPs, ktspPredictorRes$score)
#colnames(Mechanistic_KTSP) <- c("gene1", "gene2", "score")

#write.csv(Mechanistic_KTSP, file = "./Objs/KTSP/MechanisticKTSP.csv")
#print(xtable(Mechanistic_KTSP, type = "latex"), file = "./Objs/KTSP/Mechanistic.tex")

## Save the mechanistic Pairs
#MechanisticKTSP_Pairs <- c("COX7B>EZH2", "RAD21>ARHGEF11", "GTF2B>SF3B4", "TP53>MTMR4", "BTG2>CHD4", "MZF1>PITX1", "NCOR1>TRIP13")
#save(MechanisticKTSP_Pairs, file = "./Objs/KTSP/MechanisticKTSP_Pairs.rda")

### Check how many TSP
#xx <- SWAP.Filter.Wilcoxon(usedTrainGroup,usedTrainMat,featureNo=featNo)
#dim(myTSPs[myTSPs[,1] %in% xx & myTSPs[,2] %in% xx ,])
#choose(78,2)

############################################################################
### Compute the sum and find the best threshold: All training samples
ktspStatsTrainRes <- SWAP.KTSP.Statistics(inputMat = usedTrainMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTrainRes$statistics)

### Threshold
thr <- coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "best")["threshold"]
thr

### Print ROC curve local maximas
coords(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"), "local maximas")

### Plot Curve: note that you must reorder the levels!!!
### ("good" goes first, "bad" goes second, the opposite of confusionMatrix)
roc(usedTrainGroup, ktspStatsTrainRes$statistics, plot = F, print.thres=thr, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main="Mechanistic KTSP performance in the training data")

### Get predictions based on best threshold from ROC curve
usedTrainPredictionRes <- SWAP.KTSP.Classify(usedTrainMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the TRAINING set
confusionMatrix(usedTrainPredictionRes, usedTrainGroup, positive = "Mets", mode = "everything")

MCC_Mechanistic_Train <- mltools::mcc(pred = usedTrainPredictionRes, actuals = usedTrainGroup)
MCC_Mechanistic_Train

#########################################################################
#########################################################################
### Testing

## Compute the sum and find the best threshold
ktspStatsTestRes <- SWAP.KTSP.Statistics(inputMat = usedTestMat, classifier = ktspPredictorRes, CombineFunc = sum)
summary(ktspStatsTestRes$statistics)

## Plot curve
roc(usedTestGroup, ktspStatsTestRes$statistics, plot = F, print.auc=TRUE, print.auc.col="black", levels = c("No_Mets", "Mets"), direction = "<", col="blue", lwd=2, grid=TRUE, main= "Mechanistic KTSP using TF_MiR Gns")

### Get predictions based on best threshold from ROC curve
usedTestPredictionRes <- SWAP.KTSP.Classify(usedTestMat, ktspPredictorRes, DecisionFunc = function(x) sum(x) > thr)

### Resubstitution performance in the Test set
confusionMatrix(usedTestPredictionRes, usedTestGroup, positive = "Mets", mode = "everything")

MCC_Mechanistic_Test <- mltools::mcc(pred = usedTestPredictionRes, actuals = usedTestGroup)
MCC_Mechanistic_Test

## OR (calculate accuracy, sens, spec and AUC)
# SWAP.GetKTSP.PredictionStats(usedTestPredictionRes, usedTestGroup, decision_values = ktspStatsTestRes$statistics)
# 
# ### basic plots of various metrics
# library(precrec)
# precrec_obj <- evalmod(scores = ktspStatsTestRes$statistics, labels = usedTestGroup, mode="basic")
# autoplot(precrec_obj)   


############################################################################
###########################################################################
### Plot genes in the training set
## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i , function(i, x) as.character(unlist(x[i, 1:2])), x=ktspPredictorRes$TSPs)

## Assemble
dfTspTrain <- lapply(tsp, function(i,x,g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair= paste("TSP:", paste(i, collapse = "-")), out)
}, x=usedTrainMat, g=usedTrainGroup)

names(dfTspTrain) <- rownames(ktspPredictorRes$TSPs)

# Change the names of elements inside each element in dfTspTrain (For Plotting)  
for(i in seq_along(dfTspTrain)) names(dfTspTrain[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")



## Reduce
datTrain <- Reduce("rbind", dfTspTrain)

## Rename columns
#colnames(datTrain)[colnames(datTrain) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
# png("./Figs/KTSP/mechanistic.trainKTSPexprs.png", width = 3000, height = 1500, res = 200)
# bxplt <- ggplot(na.omit(datTrain), aes(x=Gene, y=Expression, fill=Group)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) +
#   facet_wrap(~pair, scales = "free", nrow = 2) +
#   theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size = 17.5), strip.text.x = element_text(face = "bold", size = 11))
# bxplt
# dev.off()

#####################
## Scatter plots
png(filename = "./Figs/KTSP/MechanisticKTSP_Train_ScatterPlot_Growth.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datTrain), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 12),
        axis.title=element_text(face="bold", size = 12),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=17.5),
        legend.title = element_text(face="bold", size=17.5),
        strip.text.x = element_text(face="bold", size=11))
sctplt
dev.off()


########################################################################
#######################################################################
### Plot genes in the testing set

## Which TSPs
i <- 1:nrow(ktspPredictorRes$TSPs)

## Assemble in a data frame
tsp <- lapply(i, function(i, x) as.character(unlist(x[i,1:2])), x=ktspPredictorRes$TSPs)

## Assemble
dfTspTest <- lapply(tsp, function(i, x, g){
  out <- data.frame(t(x[i, ]), Group=g)
  out <- cbind(pair=paste("TSP:", paste(i, collapse = "-")), out)
}, x=usedTestMat, g=usedTestGroup)


names(dfTspTest) <- rownames(ktspPredictorRes$TSPs)

# Change the names of elements inside each element in dfTspTrain 
for(i in seq_along(dfTspTest)) names(dfTspTest[[i]]) <-  c("pair", "Gene1", "Gene2", "Group")


## Reduce
datTest <- Reduce("rbind", dfTspTest)

## Rename columns
#colnames(datTest)[colnames(datTest) %in% c("variable", "value")] <- c("Gene", "Expression")

########################
## Make paired boxplot
# png("./Figs/KTSP/mechanistic.testKTSPexprs.png", width = 3000, height = 1500, res = 200)
# bxplt <- ggplot(na.omit(datTest), aes(x=Gene, y=Expression, fill=Group)) +
#   geom_boxplot(outlier.shape = NA) + 
#   geom_point(aes(group=Group), position = position_jitterdodge(), size=0.75, color=rgb(0.2,0.2,0.2,0.5)) + 
#   facet_wrap(~pair, scales = "free", nrow = 2) +
#   theme(axis.text = element_text(face = "bold", size = 12), axis.title = element_text(face = "bold", size = 12), legend.position = "bottom", legend.text = element_text(face = "bold", size = 17.5), legend.title = element_text(face = "bold", size= 17.5), strip.text.x = element_text(face = "bold", size = 11))
# bxplt
# dev.off()

#####################
## Scatter plots
png(filename = "./Figs/KTSP/MechanisticKTSP_Test_ScatterPlot_Growth.png", width=3000, height=1500, res=300)
### Prepare ggplot
sctplt <- ggplot(na.omit(datTest), aes(x=Gene1, y=Gene2, color=Group)) +
  geom_point(size=0.75) + geom_abline(intercept=0, slope=1, alpha=0.5) +
  facet_wrap(~ pair, scales="free", nrow=2) +
  theme(axis.text=element_text(face="bold", size = 12),
        axis.title=element_text(face="bold", size = 12),
        legend.position="bottom",
        legend.text=element_text(face="bold", size=17.5),
        legend.title = element_text(face="bold", size=17.5),
        strip.text.x = element_text(face="bold", size=11))
sctplt
dev.off()

#########################################################################
######################################################

## GGPLOT to compare the two classifiers
### Prepare the legend
forLegend_KTSP <- apply(rbind(
  ci(roc(usedTrainGroup, ktspStatsTrainRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTrainGroup, ktspStatsTrainUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<")),
  ci(roc(usedTestGroup, ktspStatsTestUnRes$statistics, levels = c("No_Mets", "Mets"), direction = "<"))
),  1, function(x) {
  x <- format(round(x, digits=2), nsmall=2)
  paste("AUC: ", x[[2]], ";", "95% CI: ", x[[1]], "-", x[[3]])
})


#################################################################
### ROC curves Using ggplot2

### Training
datTrn_KTSP <- melt(data.frame(
  ## Training Group
  Training= usedTrainGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Training = ktspStatsTrainUnRes$statistics,
  ## Mechanistic KTSP SUM training
  Mechanistic.Training= ktspStatsTrainRes$statistics))
### Change Colnames
colnames(datTrn_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")


### Testing
datTst_KTSP <- melt(data.frame(
  ## Testing group
  Testing= usedTestGroup,
  ## Agnostic KTSP SUM: the lowest mus be for disease status
  Agnostic.Testing=ktspStatsTestUnRes$statistics,
  ## Mechanistic KTSP SUM training
  Mechanistic.Testing=ktspStatsTestRes$statistics))
### Change Colnames
colnames(datTst_KTSP) <- c("Status", "KTSP_type", "KTSP_sum")

### Combine
dat_KTSP <- rbind(datTrn_KTSP, datTst_KTSP)

####
### Replace levels
levels(dat_KTSP$KTSP_type) <- gsub("\\.", "-", levels(dat_KTSP$KTSP_type))
levels(dat_KTSP$KTSP_type) <- paste(levels(dat_KTSP$KTSP_type), forLegend_KTSP[c(3,1,4,2)])

#################################################################
### Plot Curve
png("./Figs/KTSP/CompareAUCggplot_Invasion.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic (Invasion Pairs) vs agnostic K-TSP"
legendTitle <- paste("Mechanistic (", nrow(ktspPredictorRes$TSPs), " pairs)",
                     " Agnostic (", nrow(ktspPredictorUnRes$TSPs), " pairs)",  sep="")
### Plot
basicplot_KTSP_Invasion <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=KTSP_type,
                                              linetype = KTSP_type)) +
  geom_roc(increasing = F, cutoffs.at = seq(1,20,2)) +
  style_roc(theme = theme_grey) + ggtitle(plotTitle) +
  theme(plot.title = element_text(face="bold", size=16, hjust = 0.5),
        axis.text=element_text(face="plain", size = 11),
        axis.title=element_text(face="bold", size = 13),
        legend.justification=c(1,0),  legend.position=c(1,0),
        legend.background=element_rect(fill="lightblue1"),
        legend.text=element_text(face="plain", size = 10),
        legend.title = element_text(face="bold", size=12)) +
  scale_color_manual(legendTitle, values=rep(myCol, 2)) +
  scale_linetype_manual(legendTitle, values=rep(c("solid", "dotted"), each=2)) +
  guides(colour = guide_legend(override.aes = list(size=3)))
### Plot
basicplot_KTSP_Invasion
### Close device
dev.off()

save(basicplot_KTSP_Invasion, file = "./Objs/KTSP/BasicPlot_KTSP_Invasion.rda")


##############################################################
