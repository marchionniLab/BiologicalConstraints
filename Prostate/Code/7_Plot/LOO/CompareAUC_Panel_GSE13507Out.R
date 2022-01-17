################################################
## Mohamed Omar
## October 20, 2019
## graph panel of the ROC curves of the 4 classifiers
## GSE13507 out
#############################################
rm(list = ls())

setwd("/Users/mohamedomar/Documents/Research/Projects/Bladder")

library(ggplot2)
library("ggpubr")

#####################
## Load the 4 graphs
load("./Objs/KTSP/BasicPlot_KTSP_GSE13507Out.rda")
load("./Objs/SVM/BasicPlot_SVM_GSE13507Out.rda")
load("./Objs/RF/BasicPlot_RF_GSE13507Out.rda")
load("./Objs/XGB/BasicPlot_XGB_GSE13507Out.rda")


basicplot_KTSP_GSE13507Out$labels$title <- "KTSP (cross-study validation (GSE13507))"
basicplot_SVM_GSE13507Out$labels$title <- "SVM (cross-study validation (GSE13507))"
basicplot_RF_GSE13507Out$labels$title <- "RF (cross-study validation (GSE13507))"
basicplot_XGB_GSE13507Out$labels$title <- "XGB (cross-study validation (GSE13507))"

## Plot the panel
png("./Figs/CompareAUC_Panel_GSE13507Out.png", width = 6000, height = 4000, res = 360)
figure <- ggarrange(basicplot_KTSP_GSE13507Out, basicplot_SVM_GSE13507Out, basicplot_RF_GSE13507Out, basicplot_XGB_GSE13507Out, 
                    font.label = list(size = 20),
                    ncol = 2, nrow = 2)
figure
dev.off()
