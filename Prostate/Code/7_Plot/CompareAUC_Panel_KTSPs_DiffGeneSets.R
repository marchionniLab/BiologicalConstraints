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
## Load the 6 graphs
load("./Objs/KTSP/BasicPlot_KTSP_Combined.rda")
load("./Objs/KTSP/BasicPlot_KTSP_TF_MiRNA.rda")
load("./Objs/KTSP/BasicPlot_KTSP_TFs.rda")
load("./Objs/KTSP/BasicPlot_KTSP_ImmSurv.rda")
load("./Objs/KTSP/BasicPlot_KTSP_Angio.rda")
load("./Objs/KTSP/BasicPlot_KTSP_Apoptosis.rda")


basicplot_KTSP_Combined$labels$title <- "Cell adhesion, activation and O2 response"
basicplot_KTSP_TF_MiRNA$labels$title <- "Feedforward loops"
basicplot_KTSP_TFs$labels$title <- "gene regulatory network"
basicplot_KTSP_ImmSurv$labels$title <- "Immune Surveillance"
basicplot_KTSP_Angio$labels$title <- "Angiogenesis"
basicplot_KTSP_Apoptosis$labels$title <- "Apoptosis"


## Plot the panel
png("./Figs/CompareAUC_Panel_DifferentKTSPs.png", width = 6000, height = 4000, res = 360)
figure <- ggarrange(basicplot_KTSP_Combined, 
                    basicplot_KTSP_TF_MiRNA, 
                    basicplot_KTSP_TFs, 
                    basicplot_KTSP_ImmSurv, 
                    basicplot_KTSP_Angio,
                    basicplot_KTSP_Apoptosis,
                    font.label = list(size = 10),
                    ncol = 3, nrow = 2)
figure
dev.off()

################