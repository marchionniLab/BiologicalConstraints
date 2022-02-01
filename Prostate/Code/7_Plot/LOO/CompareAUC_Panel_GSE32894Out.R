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
load("./Objs/KTSP/BasicPlot_KTSP_GSE32894Out.rda")
load("./Objs/SVM/BasicPlot_SVM_GSE32894Out.rda")
load("./Objs/RF/BasicPlot_RF_GSE32894Out.rda")
load("./Objs/XGB/BasicPlot_XGB_GSE32894Out.rda")


basicplot_KTSP_GSE32894Out$labels$title <- "KTSP (cross-study validation (GSE32894))"
basicplot_SVM_GSE32894Out$labels$title <- "SVM (cross-study validation (GSE32894))"
basicplot_RF_GSE32894Out$labels$title <- "RF (cross-study validation (GSE32894))"
basicplot_XGB_GSE32894Out$labels$title <- "XGB (cross-study validation (GSE32894))"

## Plot the panel
png("./Figs/CompareAUC_Panel_GSE32894Out.png", width = 6000, height = 4000, res = 360)
figure <- ggarrange(basicplot_KTSP_GSE32894Out, basicplot_SVM_GSE32894Out, basicplot_RF_GSE32894Out, basicplot_XGB_GSE32894Out, 
                    font.label = list(size = 20),
                    ncol = 2, nrow = 2)
figure
dev.off()
