
rm(list = ls())

library(ggplot2)
library(dplyr)
library(tidyverse)

## Load the stuff
load("./Objs/KTSP/AUCs_KTSP_StrSampling.rda")
load("./Objs/RF/AUCs_RF_StrSampling.rda")
load("./Objs/SVM/AUCs_SVM_StrSampling.rda")
load("./Objs/XGB/AUCs_XGB_StrSampling.rda")

AllStratified <- rbind(AUCs_KTSP_StrSampling, AUCs_RF_StrSampling, AUCs_SVM_StrSampling, AUCs_XGB_StrSampling)

#AllStratified <- AllStratified[order(AllStratified$model_type, decreasing = F), ]

## Change the names of the models
AllStratified[AllStratified == "Agnostic"] <- "Agnostic (100 pairs)"
AllStratified[AllStratified == "Mechanistic"] <- "Mechanistic (100 pairs)"
AllStratified[AllStratified == "Agnostic200"] <- "Agnostic (200 pairs)"
AllStratified[AllStratified == "Agnostic500"] <- "Agnostic (500 pairs)"
AllStratified[AllStratified == "Agnostic50Feat"] <- "Agnostic (50 features)"
AllStratified[AllStratified == "Agnostic100Feat"] <- "Agnostic (100 features)"
AllStratified[AllStratified == "Agnostic500Feat"] <- "Agnostic (500 features)"


AllStratified$data_type <- factor(AllStratified$data_type, levels = c("Training", "Testing"))
AllStratified$model_type <- factor(AllStratified$model_type, levels = c("Agnostic (100 pairs)", "Agnostic (200 pairs)", "Agnostic (500 pairs)", "Agnostic (50 features)", "Agnostic (100 features)", "Agnostic (500 features)", "Mechanistic (100 pairs)"))

levels_model <-  rev(as.character(unique(sort(AllStratified$model_type))))


###########
################################
# Black and White
png(filename = "./Figs/AUCForestPlot.png", width = 3000, height = 2000, res = 300)
AllStratified %>% 
  # this is to reorder the levels
  mutate(model_type = factor(model_type, levels = levels_model)) %>% 
  # this create a var that has only the type of model
  separate(model_type, c("type", "pairs"), sep = " ", remove = FALSE, extra = "merge") %>% 
  ggplot(aes(x = AUC, y = model_type, shape =  type, fill = type)) +
  geom_linerange(aes(xmin = CI_low, xmax = CI_High)) +
  geom_point(size = 2.5) +
  geom_vline(xintercept = 0.5, linetype = 2) + 
  facet_grid(algorithm ~ data_type, scales = "free") +
  scale_x_continuous(lim = c(0.5, 1)) +
  scale_shape_manual(values = c(21, 22)) +
  scale_fill_manual(values = c("white", "grey50")) +
  labs(y = NULL,
       x = "AUC (95% Confidence Interval)") +
  theme(legend.position = "none")
dev.off()



