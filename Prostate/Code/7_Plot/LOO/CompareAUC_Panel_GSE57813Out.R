################################################
## Mohamed Omar
## October 20, 2019
## graph panel of the ROC curves of the 4 classifiers
## GSE32894 out
#############################################
rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Bladder")

library(ggplot2)
library("ggpubr")

#####################
## Load the 4 graphs
load("./Objs/KTSP/BasicPlot_KTSP_GSE57813Out.rda")
load("./Objs/SVM/BasicPlot_SVM_GSE57813Out.rda")
load("./Objs/RF/BasicPlot_RF_GSE57813Out.rda")
load("./Objs/XGB/BasicPlot_XGB_GSE57813Out.rda")


basicplot_KTSP_GSE57813Out$labels$title <- "KTSP (cross-study validation (GSE57813))"
basicplot_SVM_GSE57813Out$labels$title <- "SVM (cross-study validation (GSE57813))"
basicplot_RF_GSE57813Out$labels$title <- "RF (cross-study validation (GSE57813))"
basicplot_XGB_GSE57813Out$labels$title <- "XGB (cross-study validation (GSE57813))"

## Plot the panel
png("./Figs/CompareAUC_Panel_GSE57813Out.png", width = 6000, height = 4000, res = 360)
figure <- ggarrange(basicplot_KTSP_GSE57813Out, basicplot_SVM_GSE57813Out, basicplot_RF_GSE57813Out, basicplot_XGB_GSE57813Out, 
                    font.label = list(size = 20),
                    ncol = 2, nrow = 2)
figure
dev.off()
