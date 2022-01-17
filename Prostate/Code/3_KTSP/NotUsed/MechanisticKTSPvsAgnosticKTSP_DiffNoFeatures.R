

rm(list = ls())
#########################
# Load the performance metrics of the mechanistic KTSP and agnostic KTSP with different N of features and group them together
########################

load("./Objs/KTSP/MechanisticKTSP_Perf.rda")
colnames(MechanisticKTSP_Perf) <- paste0("MechanisticKTSP", colnames(MechanisticKTSP_Perf))

load("./Objs/KTSP/AgnosticKTSP_AllFeat_Perf.rda")
colnames(AgnosticKTSP_AllFeat_Perf) <- paste0("AgnosticKTSP_AllFeat", colnames(AgnosticKTSP_AllFeat_Perf))

load("./Objs/KTSP/AgnosticKTSP_50Feat_Perf.rda")
colnames(AgnosticKTSP_50Feat_Perf) <- paste0("AgnosticKTSP_50Feat", colnames(AgnosticKTSP_50Feat_Perf))

load("./Objs/KTSP/AgnosticKTSP_100Feat_Perf.rda")
colnames(AgnosticKTSP_100Feat_Perf) <- paste0("AgnosticKTSP_100Feat", colnames(AgnosticKTSP_100Feat_Perf))

load("./Objs/KTSP/AgnosticKTSP_200Feat_Perf.rda")
colnames(AgnosticKTSP_200Feat_Perf) <- paste0("AgnosticKTSP_200Feat", colnames(AgnosticKTSP_200Feat_Perf))

load("./Objs/KTSP/AgnosticKTSP_500Feat_Perf.rda")
colnames(AgnosticKTSP_500Feat_Perf) <- paste0("AgnosticKTSP_500Feat", colnames(AgnosticKTSP_500Feat_Perf))


All_KTSP_Performance <- cbind(MechanisticKTSP_Perf, AgnosticKTSP_AllFeat_Perf, AgnosticKTSP_50Feat_Perf, AgnosticKTSP_100Feat_Perf, AgnosticKTSP_200Feat_Perf, AgnosticKTSP_500Feat_Perf)

save(All_KTSP_Performance, file = "./Objs/KTSP/MechanisticKTSPvsAgnosticKTSP_diffNoFeatures")