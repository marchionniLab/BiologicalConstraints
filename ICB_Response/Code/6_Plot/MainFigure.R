rm(list = ls())

library(ggplot2)
library(viridis)
library(ggridges)
library(ggsci)
library(wesanderson)

## Load model comparisons: mechanistic vs agnostic (as genes) 
load("./Objs/KTSP/ModelCompare_KTSP_Pre.rda")
load("./Objs/RF/ModelCompare_RF_new_Pre.rda")
load("./Objs/SVM/ModelCompare_SVM_new_Pre.rda")
load("./Objs/XGB/ModelCompare_XGB_new_Pre.rda")


## Bind the 4 together in one data frame
AllModelCompare_ICB_AgnIndGenes <- rbind(ModelCompare_KTSP, ModelCompare_RF, ModelCompare_SVM, ModelCompare_XGB)
AllModelCompare_ICB_AgnIndGenes$data_type <- factor(AllModelCompare_ICB_AgnIndGenes$data_type, levels = c("Training", "Testing"))
AllModelCompare_ICB_AgnIndGenes$modelType <- factor(AllModelCompare_ICB_AgnIndGenes$modelType, levels = c("Agnostic_DEGs", "Mechanistic"))
levels(AllModelCompare_ICB_AgnIndGenes$modelType) <- c("Agnostic (top 50 DEGs)", "Mechanistic (27 pairs)")

############################################################################

## Load model comparisons: mechanistic vs agnostic (as pairs) 
load("./Objs/KTSP/ModelCompare_KTSP_Pre.rda") # K-TSPs is the same
load("./Objs/RF/ModelCompare_RF_AgnPairs_new_Pre.rda")
load("./Objs/SVM/ModelCompare_SVM_AgnPairs_new_Pre.rda")
load("./Objs/XGB/ModelCompare_XGB_AgnPairs_new_Pre.rda")

## Bind the 4 together in one data frame
AllModelCompare_ICB_AgnPairs <- rbind(ModelCompare_KTSP, ModelCompare_RF, ModelCompare_SVM, ModelCompare_XGB)
AllModelCompare_ICB_AgnPairs$data_type <- factor(AllModelCompare_ICB_AgnPairs$data_type, levels = c("Training", "Testing"))
AllModelCompare_ICB_AgnPairs$modelType <- factor(AllModelCompare_ICB_AgnPairs$modelType, levels = c("Agnostic_DEGs", "Agnostic_Pairs", "Mechanistic"))
levels(AllModelCompare_ICB_AgnPairs$modelType) <- c("Agnostic (top 50 DEGs)", "Agnostic (25 pairs)", "Mechanistic (27 pairs)")


###########################
# Bind the two big data frames together
AllModelCompare_ICB <- rbind(AllModelCompare_ICB_AgnIndGenes, AllModelCompare_ICB_AgnPairs)
AllModelCompare_ICB$modelType <- factor(AllModelCompare_ICB$modelType, levels = c("Mechanistic (27 pairs)", "Agnostic (top 50 DEGs)", "Agnostic (25 pairs)"))

#sel <- which(AllModelCompare_ICB$algorithm == "KTSP" & AllModelCompare_ICB$modelType == "Agnostic (25 pairs)")
#AllModelCompare_ICB$modelType[sel] <- "Agnostic (top 74 DEGs)"

############################################################################
## Plot
My_Theme = theme(
  axis.title.x = element_text(size = 10),
  axis.text.x = element_text(size = 8),
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 8),
  strip.text.x = element_text(size = 10),
  plot.title = element_text(size=12, face = "bold", hjust = 0.5)
)


## Density plot
png(filename = "./Figs/ICB_BS_AllModels_Density_new_Pre.png", width = 2500, height = 1500, res = 200)
BS_AUC_ModelCompare <- ggplot(AllModelCompare_ICB, aes(x = AUC, y = modelType, fill = modelType, alpha = data_type, height = ..ndensity..)) + 
  geom_density_ridges(stat = "density", bw = 0.8, adjust= 0.01, scale=1.2) +
  scale_x_continuous(limits = c(0.25, 1.02), breaks = seq(0.3, 1, by = 0.1), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title="Predicting the response to ICIs in melanoma patients: density of the AUC distribution of the agnostic and mechanistic models") + 
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


####################################
####################################
# for the N of pairs
load('./objs/ktsp/ModelCompare_KTSP_Npairs.rda')

ModelCompare_KTSP_Npairs$NofFeatAgn <- factor(ModelCompare_KTSP_Npairs$NofFeatAgn, levels = c('50 Genes', '100 Genes', '200 Genes', '500 Genes'))
table(ModelCompare_KTSP_Npairs$NofFeatAgn)

ModelCompare_KTSP_Npairs$modelType <- factor(ModelCompare_KTSP_Npairs$modelType, levels = c('Mechanistic', 'Agnostic'))
table(ModelCompare_KTSP_Npairs$modelType)

colnames(ModelCompare_KTSP_Npairs) <- c('N of output pairs', 'model type', 'N of training features')

png(filename = "./Figs/ICB_BS_AllModels_Npairs.png", width = 3000, height = 1500, res = 200)
BS_AUC_ModelCompare <- ggplot(ModelCompare_KTSP_Npairs, aes(x = `N of output pairs`, y = `N of training features`, fill = `model type`, height = ..ndensity..)) + 
  geom_density_ridges(stat = "density", scale=1.2, alpha = 0.8) +
  scale_x_continuous(limits = c(0, 55), breaks = seq(0, 50, by = 5), expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title="The number of pairs returned by the K-TSPs model in the melanoma response to ICIs case") + 
  ylab("Model Type") +
  My_Theme +
  #scale_fill_manual(values = c("red3", "#78B7C5", "#3B9AB2")) +
  #scale_fill_jco() +
  #scale_alpha_manual(values = c(0.4, 1))+
  scale_fill_viridis(discrete = TRUE)+
  coord_cartesian(clip = "off") +
  #facet_wrap(~algorithm, dir = "v", ncol = 1, scales = "free_y") + 
  theme(axis.text = element_text(face = "bold"), 
        panel.background = element_rect(fill = "gray93"), 
        plot.background = element_rect(fill = "white")) + 
  labs(fill = "model type")

BS_AUC_ModelCompare
dev.off()
