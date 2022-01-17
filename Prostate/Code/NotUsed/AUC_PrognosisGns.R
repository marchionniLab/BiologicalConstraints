###############################################################################
### Mohamed Omar
### 01/07/2019
### Goal : Creating the restricted K-TSP classifier
### Using : Prognosis genes
### AUC_Sensitivity plot to choose the best K
#################################################################################

rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

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

###################################################################
## Load Metastasis Data
load("./Objs/MetastasisData.rda")


Prognosis_genes1 <- read.csv("/Users/mohamedomar/Desktop/Prostate_Prognosis.csv")
Prognosis_genes1 <- apply(Prognosis_genes1, 2, as.character)

Prognosis_genes2 <- read.csv("/Users/mohamedomar/Desktop/Prostate_Prognosis2.csv")
Prognosis_genes2 <- apply(Prognosis_genes2, 2, as.character)


myTSPs <- rbind(Prognosis_genes1, Prognosis_genes2)
### Common genes
keepGns <- intersect(as.vector(myTSPs), rownames(mixTrainMat))

## Normalization between arrays
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
boxplot(usedTrainMat, outline =FALSE)

usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")
boxplot(usedTestMat, outline = FALSE)

#
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

### For the TSP
myTSPs <- myTSPs[myTSPs[,1] %in% keepGns & myTSPs[,2] %in% keepGns , ]

### Set Feature number and max k
featN <- nrow(usedTrainMat) #nrow(usedTrainMat)
ktsp <- c(1:50)

### Train a classifier with the desired maximum number of TSPs
allTSPs <- SWAP.Train.KTSP(usedTrainMat, usedTrainGroup, krange=ktsp,
                           FilterFunc = SWAP.Filter.Wilcoxon,
                           featureNo=featN,
                           RestrictedPairs=myTSPs)
allTSPs

### Extract statistics for increasing number k of TSPs
statsTSPs <- lapply(2:length(allTSPs$score), function(i, ...) {
  ktsp <- list(name=paste(i, "TSPs", sep=""), TSPs=allTSPs$TSPs[1:i,],
               score=allTSPs$score[1:i], labels=allTSPs$labels)
  SWAP.KTSP.Statistics(usedTrainMat, ktsp, sum)$statistics
}, allTSPs, usedTrainMat)

###########################################################################
### Extract statistics: AUC (In the testing set)
###########################################################################

### Compute AUC for plotting
zz <- 1:length(statsTSPs)
str(zz)

pp <- sapply(statsTSPs[zz], function(x, ...) {
  outAUC <- ci(roc(usedTrainGroup, x, levels=c("No_Mets", "Mets"),
                   direction=">",))[2]
  outOther <- coords(roc(usedTrainGroup, x, levels=c("No_Mets", "Mets"),
                         direction=">"), "best")
  if ( is.matrix(outOther) ) outOther <- outOther[ , which.max(outOther["sensitivity",])]
  c(AUC=outAUC, outOther)
})


############################################################################
### Plot
png("./Figs/Prognosis_genes/mechanistic.AUCbyKTSP.png", width = 3000, height = 1500, res = 200)
## Prepare pallette and layout
myCol <- c(brewer.pal(3, "Dark2")[c(3)], rgb(0.11, 0.62, 0.47, 0.65), rgb(0.85, 0.37, 0.1, 0.65))
par(margin= c(5,5,3,5), mgp=c(3,0.75,0))
## Prepeare coordinates
xCoord <- 2:(length(pp["AUC",])+1)
yRange <- range(pp["AUC", ]) + c(-0.01,0.01)
## Plot
plot(xCoord, pp["AUC",], las=1, cex.lab=1.5, cex.main=2, xaxp=c(0,100,100), yaxp=c(0,1,100), ylim=yRange, xlab="Number of combined TSPs", ylab="AUC",main="AUC versus Number of TSPs")
## Add grid
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "grey90")
abline(v=xCoord, h=axTicks(1,c(0,1,100)), col="white", lty=3)
## Add lines and points
points(xCoord, pp["AUC", ], col=myCol[1], type = "b", lwd=4, pch= 16)
text(xCoord, pp["AUC", ], xCoord, col = rgb(0.1,0.1,0.1,0.85), offset = c(1.5), pos = 3, cex = 1.25)
## Add line points on a new scale
par(new = T)
par(new = T)
plot(xCoord, rep(range(pp[c("specificity", "sensitivity"),]),
                 length.out = length(xCoord)), col=NA, axes=F, xlab=NA, ylab=NA)
points(xCoord, pp["sensitivity",], col=myCol[2], type = "b", lwd=4, pch=16)
points(xCoord, pp["specificity",], col=myCol[3], type = "b", lwd=4, pch=16)
## Label axis
axis(side = 4, las=1, yaxp=c(0,1,50), crt=90)
mtext(side = 4, line = 3, "Sensitivity and 1-Specificity", cex = 1.5)
## Add legend
legend("bottom", legend = c("AUC", "Sensitivity", "1-Specificity"), col = myCol, lwd = 20, cex=1.75, bg="grey90")
dev.off()

###############################################################################
### Plot
png("./Figs/Prognosis_genes/mechanistic.AUCbyKTSP.fixed.png", width = 3000, height = 1500, res = 200)
## Prepare pallette and layout
myCol <- c(brewer.pal(3, "Dark2")[c(3)], rgb(0.11, 0.62, 0.47, 0.65), rgb(0.85, 0.37, 0.1, 0.65))
par(margin= c(5,5,3,5), mgp=c(3,0.75,0))
## Prepare coordinates
xCoord <- 2:(length(pp["AUC", ])+1)
yRange <- range(pp[-grep("threshold", rownames(pp)), ]) + c(-0.01, 0.01)
## Plot
plot(xCoord, pp["AUC",], las=1, cex.lab=1.5, cex.main=2, xaxp=c(0,100,100), yaxp=c(0,1,50), ylim=yRange, xlab="Number of combined TSPs", ylab="AUC", main="AUC versus Number of TSPs")
## Add grid
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
abline(v=xCoord, h=axTicks(1, c(0, 1, 50)), col="white", lty=3)
### Add lines and points
points(xCoord, pp["AUC",], col=myCol[1], type="b", lwd=4, pch=16)
text(xCoord, pp["AUC",], xCoord, col=rgb(0.1,0.1,0.1,0.85),
     offset=c(1.5), pos=3, cex=1.25)
points(xCoord, pp["sensitivity",], col=myCol[2],type="b", lwd=4, pch=16)
points(xCoord, pp["specificity",], col=myCol[3],type="b", lwd=4, pch=16)
### Add Legend
legend("bottom", legend=c("AUC", "Sensitivity", "1-Specificity"),
       col=myCol, lwd=20, , cex=1.75, bg="grey90")
dev.off()

#######################################
# Create an HTML file
#rmarkdown::render(input = "./Code/3Restricted_KTSP_AUC_TF_MiR.R", output_dir = "./HTML", output_format = "html_document")

