################################################
## Mohamed Omar
## October 20, 2019
## graph panel of the ROC curves of the 4 classifiers
## EMTAB-4321 out
#############################################
rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Bladder")

library(ggplot2)
library("ggpubr")

#####################
## Load the 4 graphs
load("./Objs/KTSP/BasicPlot_KTSP_EMTABOut.rda")
load("./Objs/SVM/BasicPlot_SVM_EMTABOut.rda")
load("./Objs/RF/BasicPlot_RF_EMTABOut.rda")
load("./Objs/XGB/BasicPlot_XGB_EMTABOut.rda")


basicplot_KTSP_EMTABOut$labels$title <- "KTSP (cross-study validation (EMTAB-4321))"
basicplot_SVM_EMTABOut$labels$title <- "SVM (cross-study validation (EMTAB-4321))"
basicplot_RF_EMTABOut$labels$title <- "RF (cross-study validation (EMTAB-4321))"
basicplot_XGB_EMTABOut$labels$title <- "XGB (cross-study validation (EMTAB-4321))"

## Plot the panel
png("./Figs/CompareAUC_Panel_EMTABOut.png", width = 6000, height = 4000, res = 360)
figure <- ggarrange(basicplot_KTSP_EMTABOut, basicplot_SVM_EMTABOut, basicplot_RF_EMTABOut, basicplot_XGB_EMTABOut, 
                    font.label = list(size = 20),
                    ncol = 2, nrow = 2)
figure
dev.off()
