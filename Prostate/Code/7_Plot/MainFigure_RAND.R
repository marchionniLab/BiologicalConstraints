rm(list = ls())

library(ggplot2)
library(viridis)
library(ggridges)
library(ggsci)

## Load model comparisons: mechanistic vs agnostic (as genes) 
load("./Objs/KTSP/ModelCompare_KTSP_RAND.rda")
load("./Objs/RF/ModelCompare_RF_RAND.rda")
load("./Objs/SVM/ModelCompare_SVM_50_RAND.rda")
load("./Objs/XGB/ModelCompare_XGB_RAND.rda")


## Bind the 4 together in one data frame
AllModelCompare_Prostate_AgnIndGenes <- rbind(ModelCompare_KTSP, ModelCompare_RF, ModelCompare_SVM, ModelCompare_XGB)
AllModelCompare_Prostate_AgnIndGenes$data_type <- factor(AllModelCompare_Prostate_AgnIndGenes$data_type, levels = c("Training", "Testing"))
AllModelCompare_Prostate_AgnIndGenes$modelType <- factor(AllModelCompare_Prostate_AgnIndGenes$modelType, levels = c("Mechanistic", "Random_genes"))
levels(AllModelCompare_Prostate_AgnIndGenes$modelType) <- c("Mechanistic", "Random genes")
############################################################################

## Load model comparisons: mechanistic vs agnostic (as pairs) 
# load("./Objs/KTSP/ModelCompare_KTSP.rda") # K-TSPs is the same
# load("./Objs/RF/ModelCompare_RF_AgnPairs_new.rda")
# load("./Objs/SVM/ModelCompare_SVM_AgnPairs_new.rda")
# load("./Objs/XGB/ModelCompare_XGB_AgnPairs_new.rda")
# 
# ## Bind the 4 together in one data frame
# AllModelCompare_Prostate_AgnPairs <- rbind(ModelCompare_KTSP, ModelCompare_RF, ModelCompare_SVM, ModelCompare_XGB)
# AllModelCompare_Prostate_AgnPairs$data_type <- factor(AllModelCompare_Prostate_AgnPairs$data_type, levels = c("Training", "Testing"))
# AllModelCompare_Prostate_AgnPairs$modelType <- factor(AllModelCompare_Prostate_AgnPairs$modelType, levels = c("Agnostic_DEGs", "Agnostic_Pairs", "Mechanistic"))
# levels(AllModelCompare_Prostate_AgnPairs$modelType) <- c("Agnostic (top 74 DEGs)", "Agnostic (37 pairs)", "Mechanistic (37 FFLs)")
# 
# 
# ###########################
# # Bind the two big data frames together
# AllModelCompare_Prostate <- rbind(AllModelCompare_Prostate_AgnIndGenes, AllModelCompare_Prostate_AgnPairs)
# AllModelCompare_Prostate$modelType <- factor(AllModelCompare_Prostate$modelType, levels = c("Mechanistic (37 FFLs)", "Agnostic (top 74 DEGs)", "Agnostic (37 pairs)"))
# 
# sel <- which(AllModelCompare_Prostate$algorithm == "KTSP" & AllModelCompare_Prostate$modelType == "Agnostic (37 pairs)")
# AllModelCompare_Prostate$modelType[sel] <- "Agnostic (top 74 DEGs)"

############################################################################
## Plot
My_Theme = theme(
  axis.title.x = element_text(size = 10),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 8),
  strip.text.x = element_text(size = 10),
  plot.title = element_text(size=11, face = "bold", hjust = 0.5)
)


## Density plot
png(filename = "./Figs/Prostate_BS_AllModels_Density_RAND.png", width = 2500, height = 1500, res = 200)
BS_AUC_ModelCompare <- ggplot(AllModelCompare_Prostate_AgnIndGenes, aes(x = AUC, y = modelType, fill = modelType, alpha = data_type, height = ..ndensity..)) + 
  geom_density_ridges(stat = "density", bw = 0.8, adjust= 0.01, scale=1.2) +
  scale_x_continuous(limits = c(0.45, 1.02), breaks = seq(0.5, 1, by = 0.1), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title="Predicting prostate cancer metastasis: density of the AUC distribution of the mechanistic models vs models trained on random genes") + 
  ylab("Model Type") +
  My_Theme +
  #scale_fill_manual(values = c("red3", "#78B7C5", "#3B9AB2")) +
  #scale_fill_jco() +
  scale_alpha_manual(values = c(0.4, 1))+
  scale_fill_viridis(discrete = TRUE)+
  coord_cartesian(clip = "off") +
  facet_wrap(~algorithm, dir = "v", ncol = 1, scales = "free_y") + theme(axis.text = element_text(face = "bold"), 
                                                                         panel.background = element_rect(fill = "gray93"), 
                                                                         plot.background = element_rect(fill = "white")) +labs(fill = "model type", alpha = "data type")

BS_AUC_ModelCompare
dev.off()

