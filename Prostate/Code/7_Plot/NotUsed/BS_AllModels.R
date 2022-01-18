rm(list = ls())

library(ggplot2)

load("./Objs/KTSP/ModelCompare_KTSP.rda")
load("./Objs/RF/ModelCompare_RF.rda")
load("./Objs/SVM/ModelCompare_SVM.rda")
load("./Objs/XGB/ModelCompare_XGB.rda")

AllModelCompare_Prostate <- rbind(ModelCompare_KTSP, ModelCompare_RF, ModelCompare_SVM, ModelCompare_XGB)
save(AllModelCompare_Prostate, file = "./Objs/BS_AllModels_Prostate.rda")


#######
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 10),
  axis.text.x = element_text(size = 10),
  axis.title.y = element_text(size = 10),
  axis.text.y = element_text(size = 10),
  plot.title = element_text(size=10)
)

png(filename = "./Figs/BS_AllModels.png", width = 3000, height = 1500, res = 300)
BS_AUC_ModelCompare <- ggplot(AllModelCompare, aes(AUC, fill = modelType, linetype = data_type)) + 
  geom_density(alpha = 0.5) +
  #scale_x_continuous(limits = c(0.4, 1.1)) +
  scale_y_continuous(limits = c(0, 60)) +
  labs(title="AUC distribution of the agnostic and mechanistic models") + My_Theme+
  facet_grid(algorithm ~ .)
BS_AUC_ModelCompare
dev.off()

