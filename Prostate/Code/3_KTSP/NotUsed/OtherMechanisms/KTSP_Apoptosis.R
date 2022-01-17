
rm(list = ls()) 

setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")


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
Genes <- load("/Users/mohamedomar/Desktop/Genes/MECHANISM/CancerHallmarks/Apoptosis/objs/apoptosisPairs.v6.1.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

### Common genes
myTSPs <- apoptosisPairs
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
ktsp <- c(3:25) #24
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
###############################################################
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
dat_KTSP$Status <- as.numeric(dat_KTSP$Status)-1

####
### Replace levels
levels(dat_KTSP$KTSP_type) <- gsub("\\.", "-", levels(dat_KTSP$KTSP_type))
levels(dat_KTSP$KTSP_type) <- paste(levels(dat_KTSP$KTSP_type), forLegend_KTSP[c(3,1,4,2)])

#################################################################
### Plot Curve
png("./Figs/KTSP/CompareAUCggplot_Apoptosis.png",
    width=3000, height=3000, res=360)
### Color
myCol <- brewer.pal(3, "Dark2")[c(2,1)]
### Plot and legend titles
plotTitle <- "AUC mechanistic (Apoptosis pairs) vs agnostic K-TSP"
legendTitle <- paste("Mechanistic (", nrow(ktspPredictorRes$TSPs), " pairs)",
                     " Agnostic (", nrow(ktspPredictorUnRes$TSPs), " pairs)",  sep="")
### Plot
basicplot_KTSP_Apoptosis <- ggplot(dat_KTSP, aes(d=Status, m=KTSP_sum, color=KTSP_type,
                                       linetype = KTSP_type)) +
  geom_roc(cutoffs.at = seq(1,20,2)) +
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
basicplot_KTSP_Apoptosis
### Close device
dev.off()

save(basicplot_KTSP_Apoptosis, file = "./Objs/KTSP/BasicPlot_KTSP_Apoptosis.rda")


##############################################################
#rmarkdown::render(input = "./Code/TF_MiR/10KTSP_TF_MiR.R", output_dir = "./HTML", output_format = "html_document")