################################################
## Mohamed Omar
## October 20, 2019
## graph panel of the ROC curves of the 4 classifiers
## pmid out
#############################################
rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Bladder")

library(ggplot2)
library("ggpubr")

#####################
## Load the 4 graphs
load("./Objs/KTSP/BasicPlot_KTSP_pmidOut.rda")
load("./Objs/SVM/BasicPlot_SVM_pmidOut.rda")
load("./Objs/RF/BasicPlot_RF_pmidOut.rda")
load("./Objs/XGB/BasicPlot_XGB_pmidOut.rda")


basicplot_KTSP_pmidOut$labels$title <- "KTSP (cross-study validation (pmid15930337))"
basicplot_SVM_pmidOut$labels$title <- "SVM (cross-study validation (pmid15930337))"
basicplot_RF_pmidOut$labels$title <- "RF (cross-study validation (pmid15930337))"
basicplot_XGB_pmidOut$labels$title <- "XGB (cross-study validation (pmid15930337))"

## Plot the panel
png("./Figs/CompareAUC_Panel_pmidOut.png", width = 6000, height = 4000, res = 360)
figure <- ggarrange(basicplot_KTSP_pmidOut, basicplot_SVM_pmidOut, basicplot_RF_pmidOut, basicplot_XGB_pmidOut, 
                    font.label = list(size = 20),
                    ncol = 2, nrow = 2)
figure
dev.off()
