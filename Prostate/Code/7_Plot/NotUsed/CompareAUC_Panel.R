################################################
## Mohamed Omar
## October 20, 2019
## graph panel of the ROC curves of the 4 classifiers

#############################################
rm(list = ls())

setwd("/Volumes/Macintosh/Dropbox (MechPred)/MechPred/User/Mohamed/MechanisticModels/Prostate")

library(ggplot2)
library("ggpubr")

#####################
## Load the 4 graphs
load("./Objs/KTSP/BasicPlot_KTSP_Combined.rda")
load("./Objs/SVM/BasicPlot_SVM.rda")
load("./Objs/RF/BasicPlot_RF_Pairs.rda")
load("./Objs/XGB/BasicPlot_XGB_Pairs.rda")


basicplot_KTSP_Combined$labels$title <- "KTSP"
basicplot_SVM$labels$title <- "SVM"
basicplot_RF_Pairs$labels$title <- "RF"
basicplot_XGB_Pairs$labels$title <- "XGB"

## Plot the panel
png("./Figs/CompareAUC_Panel_Combined.png", width = 6000, height = 4000, res = 360)
figure <- ggarrange(basicplot_KTSP_Combined, basicplot_SVM, basicplot_RF_Pairs, basicplot_XGB_Pairs, 
                   font.label = list(size = 20),
                    ncol = 2, nrow = 2)
figure
dev.off()

################