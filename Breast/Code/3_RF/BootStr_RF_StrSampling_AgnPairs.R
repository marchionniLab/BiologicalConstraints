###################################################
## RF
# Mech

rm(list = ls())

## Load necessary packages
library(caret)
library(limma)
library(pROC)
library(genefilter)
library(DMwR)
library(randomForest)
library(patchwork)
library(boot)

# ## Load data
# load("./Objs/KTSP/TNBC_KTSP_STATs_Mechanistic_NotchAndMyc2.rda")
# load("./Objs/ChemoDataNew.rda")
# 
# 
# ### Quantile normalize
# usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
# usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")
# 
# ####
# usedTrainGroup <- mixTrainGroup
# usedTestGroup <- mixTestGroup
# 
# ######
# predictor_data_Train_Mech <- t(KTSP_STATs_Train_Mechanistic)
# predictor_data_Test_Mech <- t(KTSP_STATs_Test_Mechanistic)
# 
# DataMech_Train <- cbind(predictor_data_Train_Mech, usedTrainGroup)
# DataMech_Train <- as.data.frame(DataMech_Train)
# DataMech_Train$usedTrainGroup <- as.factor(DataMech_Train$usedTrainGroup)
# levels(DataMech_Train[, "usedTrainGroup"]) <- c("Sensitive", "Resistant")
# 
# names(DataMech_Train) <- make.names(names(DataMech_Train))
# 
# colnames(predictor_data_Train_Mech) <- make.names(colnames(predictor_data_Train_Mech))
# colnames(predictor_data_Test_Mech) <- make.names(colnames(predictor_data_Test_Mech))
# 
# # The function for bootstraping
# RF_Strap <- function(data, indices) {
#   d <- data[indices, ] # allows boot to select sample
#   # Select the minimum sample size
#   tmp <- as.vector(table(d$usedTrainGroup))
#   num_classes <- length(tmp)
#   min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
#   sampsizes <- rep(min_size,num_classes)
#   RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
#   Importance_Mech <- randomForest::importance(RF, scale=FALSE, type = 2)
#   Importance_Mech <- Importance_Mech[order(Importance_Mech, decreasing = TRUE), ]
#   Importance_Mech <- Importance_Mech[Importance_Mech > 0]
#   N_ImportanVariables <- length(Importance_Mech)
#   PhenoTrain <- d$usedTrainGroup
#   PredictorTrainData <- d
#   PredictorTrainData$usedTrainGroup <- NULL
#   train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
#   test_preds <- predict(RF, newdata = predictor_data_Test_Mech, type = "vote")
#   ROCTrainMech <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
#   ROCTestMech <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
#   return(c(ROCTrainMech$auc, ROCTestMech$auc, N_ImportanVariables))
# }
# 
# 
# set.seed(333)
# bootobjectMech <- boot(data= DataMech_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 
# 
# AUCs_RF_Mech <- bootobjectMech$t
# colnames(AUCs_RF_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

load('./objs/RF/RF_MechBootObject.rda')
####################################################################################
####################################################################################
## RF
# Agnostic 

# 25 pairs

load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_25.rda")
load("./Objs/ChemoDataNew.rda")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

####
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

######
predictor_data_Train_Agnostic <- t(KTSP_STATs_Train_Agnostic_25)
predictor_data_Test_Agnostic <- t(KTSP_STATs_Test_Agnostic_25)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("Sensitive", "Resistant")

names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_25 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15)

####################################################################################
####################################################################################
## RF
# Agnostic 

# 50 pairs

load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_50.rda")
load("./Objs/ChemoDataNew.rda")


### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#########
predictor_data_Train_Agnostic <- t(KTSP_STATs_Train_Agnostic_50)
predictor_data_Test_Agnostic <- t(KTSP_STATs_Test_Agnostic_50)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("Sensitive", "Resistant")

names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)  
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_50 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

####################################################################################
####################################################################################
## RF
# Agnostic 

# 100 pairs

load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_100.rda")
load("./Objs/ChemoDataNew.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

#########
predictor_data_Train_Agnostic <- t(KTSP_STATs_Train_Agnostic_100)
predictor_data_Test_Agnostic <- t(KTSP_STATs_Test_Agnostic_100)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("Sensitive", "Resistant")

names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)  
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_100 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

####################################################################################
####################################################################################
## RF
# Agnostic 

# 119 pairs
# 
# load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_119.rda")
# load("./Objs/ChemoDataNew.rda")
# 
# ### Quantile normalize
# usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
# usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")
# 
# usedTrainGroup <- mixTrainGroup
# usedTestGroup <- mixTestGroup
# 
# #########
# predictor_data_Train_Agnostic <- t(KTSP_STATs_Train_Agnostic_119)
# predictor_data_Test_Agnostic <- t(KTSP_STATs_Test_Agnostic_119)
# 
# DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
# DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
# DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
# levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("Sensitive", "Resistant")
# 
# names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))
# 
# colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
# colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))
# 
# # The function for bootstraping
# RF_Strap <- function(data, indices) {
#   d <- data[indices, ] # allows boot to select sample
#   # Select the minimum sample size
#   tmp <- as.vector(table(d$usedTrainGroup))
#   num_classes <- length(tmp)
#   min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
#   sampsizes <- rep(min_size,num_classes)
#   RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)  
#   Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
#   Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
#   Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
#   N_ImportanVariables <- length(Importance_Agnostic)
#   PhenoTrain <- d$usedTrainGroup
#   PredictorTrainData <- d
#   PredictorTrainData$usedTrainGroup <- NULL
#   train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
#   test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
#   ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
#   ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
#   return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
# }
# 
# 
# set.seed(333)
# bootobjectAgnostic_119 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 
# 
####################################################################################
####################################################################################
## RF
# Agnostic 

# 250 pairs

load("./Objs/KTSP/TNBC_KTSP_STATs_Agnostic_250.rda")
load("./Objs/ChemoDataNew.rda")

### Quantile normalize
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")

usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup

predictor_data_Train_Agnostic <- t(KTSP_STATs_Train_Agnostic_250)
predictor_data_Test_Agnostic <- t(KTSP_STATs_Test_Agnostic_250)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("Sensitive", "Resistant")

names(DataAgnostic_Train) <- make.names(names(DataAgnostic_Train))

colnames(predictor_data_Train_Agnostic) <- make.names(colnames(predictor_data_Train_Agnostic))
colnames(predictor_data_Test_Agnostic) <- make.names(colnames(predictor_data_Test_Agnostic))

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)  
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("Sensitive", "Resistant"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic_250 <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 15) 

####################################################################################
################################################################################

## Save all bootobjects
save(bootobjectMech, bootobjectAgnostic_25, bootobjectAgnostic_50, bootobjectAgnostic_100, bootobjectAgnostic_250, file = "./Objs/RF/RFBootObjects_NotchAndMyc_AgnPairs_new.rda")

## load
load("./Objs/RF/RFBootObjects_NotchAndMyc_AgnPairs_new.rda")

##################################################################################
################################################################################

## Work with Agnostic bootobject 25 vs mechanistic

AUCs_RF_Agnostic_25 <- bootobjectAgnostic_25$t
colnames(AUCs_RF_Agnostic_25) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


AUCs_RF_Mech <- bootobjectMech$t
colnames(AUCs_RF_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_25 <- AUCs_RF_Agnostic_25[, "AUC_Train"] - AUCs_RF_Agnostic_25[, "AUC_Test"]
quantile(DiffAgnostic_25, c(0.025, 0.975))

## Calculate the difference and CI of the difference between mechanistic training and mechanistic testing
DiffMech <- AUCs_RF_Mech[, "AUC_Train"] - AUCs_RF_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUC_Train <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Train"])
AgnosticAUC_Train_25 <- data.frame(AUC = AUCs_RF_Agnostic_25[, "AUC_Train"])

MechanisticAUC_Train$modelType <- "Mechanistic"
AgnosticAUC_Train_25$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_25 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_25)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUC_Test <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Test"])
AgnosticAUC_Test_25 <- data.frame(AUC = AUCs_RF_Agnostic_25[, "AUC_Test"])

MechanisticAUC_Test$modelType <- "Mechanistic"
AgnosticAUC_Test_25$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_25 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_25)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_25$data_type <- "Training"
ModelCompareAUC_Test_25$data_type <- "Testing"

ModelCompareAUC_Train_25$NofFeatAgn <- "25_Pairs"
ModelCompareAUC_Test_25$NofFeatAgn <- "25_Pairs"


#############################################################################
## Save for the main figure
# ModelCompare_RF <- rbind(ModelCompareAUC_Train_25, ModelCompareAUC_Test_25)
# ModelCompare_RF$algorithm <- "RF"
# save(ModelCompare_RF, file = "./Objs/RF/ModelCompare_RF_AgnPairs.rda")

#####################################################################
############################################################
#Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/RF/ModelCompareAUC_50_new.rda")

# Combine
ModelCompareAUC_Train_25_50 <- rbind(ModelCompareAUC_Train_25, ModelCompareAUC_Train_50)
ModelCompareAUC_Test_25_50 <- rbind(ModelCompareAUC_Test_25, ModelCompareAUC_Test_50)

##########
###################################################################################3
###################################################################################
## Work with Agnostic bootobject 50 vs mechanistic

AUCs_RF_Agnostic_50 <- bootobjectAgnostic_50$t
colnames(AUCs_RF_Agnostic_50) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_RF_Mech <- bootobjectMech$t
colnames(AUCs_RF_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_50 <- AUCs_RF_Agnostic_50[, "AUC_Train"] - AUCs_RF_Agnostic_50[, "AUC_Test"]
quantile(DiffAgnostic_50, c(0.025, 0.975))

## Calculate the difference and CI of the difference between mechanistic training and mechanistic testing
DiffMech <- AUCs_RF_Mech[, "AUC_Train"] - AUCs_RF_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUC_Train <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Train"])
AgnosticAUC_Train_50 <- data.frame(AUC = AUCs_RF_Agnostic_50[, "AUC_Train"])

MechanisticAUC_Train$modelType <- "Mechanistic"
AgnosticAUC_Train_50$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_50 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_50)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUC_Test <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Test"])
AgnosticAUC_Test_50 <- data.frame(AUC = AUCs_RF_Agnostic_50[, "AUC_Test"])

MechanisticAUC_Test$modelType <- "Mechanistic"
AgnosticAUC_Test_50$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_50 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_50)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_50$data_type <- "Training"
ModelCompareAUC_Test_50$data_type <- "Testing"

ModelCompareAUC_Train_50$NofFeatAgn <- "50_Pairs"
ModelCompareAUC_Test_50$NofFeatAgn <- "50_Pairs"

##############################################################################
## Save for the main figure
# ModelCompare_RF <- rbind(ModelCompareAUC_Train_50, ModelCompareAUC_Test_50)
# ModelCompare_RF$algorithm <- "RF"
# save(ModelCompare_RF, file = "./Objs/RF/ModelCompare_RF_AgnPairs.rda")

#####################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/RF/ModelCompareAUC_100_new.rda")

# Combine
ModelCompareAUC_Train_50_100 <- rbind(ModelCompareAUC_Train_50, ModelCompareAUC_Train_100)
ModelCompareAUC_Test_50_100 <- rbind(ModelCompareAUC_Test_50, ModelCompareAUC_Test_100)

###################################################################################3
###################################################################################
## Work with Agnostic bootobject 100 vs mechanistic

AUCs_RF_Agnostic_100 <- bootobjectAgnostic_100$t
colnames(AUCs_RF_Agnostic_100) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_100 <- AUCs_RF_Agnostic_100[, "AUC_Train"] - AUCs_RF_Agnostic_100[, "AUC_Test"]
quantile(DiffAgnostic_100, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_100 <- data.frame(AUC = AUCs_RF_Agnostic_100[, "AUC_Train"])

AgnosticAUC_Train_100$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_100 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_100)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_100 <- data.frame(AUC = AUCs_RF_Agnostic_100[, "AUC_Test"])

AgnosticAUC_Test_100$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_100 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_100)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_100$data_type <- "Training"
ModelCompareAUC_Test_100$data_type <- "Testing"

ModelCompareAUC_Train_100$NofFeatAgn <- "100_Pairs"
ModelCompareAUC_Test_100$NofFeatAgn <- "100_Pairs"

#####################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/RF/ModelCompareAUC_200_new.rda")

# Combine
ModelCompareAUC_Train_100_200 <- rbind(ModelCompareAUC_Train_100, ModelCompareAUC_Train_200)
ModelCompareAUC_Test_100_200 <- rbind(ModelCompareAUC_Test_100, ModelCompareAUC_Test_200)

###################################################################################3
###################################################################################
## Work with Agnostic bootobject 119 vs mechanistic
# 
# AUCs_RF_Agnostic_119 <- bootobjectAgnostic_119$t
# colnames(AUCs_RF_Agnostic_119) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")
# 
# 
# ## Calculate the difference and CI of the difference between agnostic training and agnostic testing
# DiffAgnostic_119 <- AUCs_RF_Agnostic_119[, "AUC_Train"] - AUCs_RF_Agnostic_119[, "AUC_Test"]
# quantile(DiffAgnostic_119, c(0.025, 0.975))
# 
# ## Plot the distributions of the AUCs from both methods in the training data
# AgnosticAUC_Train_119 <- data.frame(AUC = AUCs_RF_Agnostic_119[, "AUC_Train"])
# 
# AgnosticAUC_Train_119$modelType <- "Agnostic_Pairs"
# 
# ModelCompareAUC_Train_119 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_119)
# 
# ## Plot the distributions of the AUCs from both methods in the testing data
# AgnosticAUC_Test_119 <- data.frame(AUC = AUCs_RF_Agnostic_119[, "AUC_Test"])
# 
# AgnosticAUC_Test_119$modelType <- "Agnostic_Pairs"
# 
# ModelCompareAUC_Test_119 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_119)
# 
# ## Save the AUCs in the training and testing data
# ModelCompareAUC_Train_119$data_type <- "Training"
# ModelCompareAUC_Test_119$data_type <- "Testing"
# 
# ModelCompareAUC_Train_119$NofFeatAgn <- "119_Pairs"
# ModelCompareAUC_Test_119$NofFeatAgn <- "119_Pairs"

#####################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
# load("./Objs/RF/ModelCompareAUC_238.rda")
# 
# # Combine
# ModelCompareAUC_Train_119_238 <- rbind(ModelCompareAUC_Train_119, ModelCompareAUC_Train_238)
# ModelCompareAUC_Test_119_238 <- rbind(ModelCompareAUC_Test_119, ModelCompareAUC_Test_238)


###################################################################################3
###################################################################################
## Work with Agnostic bootobject 500 vs mechanistic

AUCs_RF_Agnostic_250 <- bootobjectAgnostic_250$t
colnames(AUCs_RF_Agnostic_250) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


## Calculate the difference and CI of the difference between agnostic training and agnostic testing
DiffAgnostic_250 <- AUCs_RF_Agnostic_250[, "AUC_Train"] - AUCs_RF_Agnostic_250[, "AUC_Test"]
quantile(DiffAgnostic_250, c(0.025, 0.975))

## Plot the distributions of the AUCs from both methods in the training data
AgnosticAUC_Train_250 <- data.frame(AUC = AUCs_RF_Agnostic_250[, "AUC_Train"])

AgnosticAUC_Train_250$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Train_250 <- rbind(MechanisticAUC_Train, AgnosticAUC_Train_250)

## Plot the distributions of the AUCs from both methods in the testing data
AgnosticAUC_Test_250 <- data.frame(AUC = AUCs_RF_Agnostic_250[, "AUC_Test"])

AgnosticAUC_Test_250$modelType <- "Agnostic_Pairs"

ModelCompareAUC_Test_250 <- rbind(MechanisticAUC_Test, AgnosticAUC_Test_250)

## Save the AUCs in the training and testing data
ModelCompareAUC_Train_250$data_type <- "Training"
ModelCompareAUC_Test_250$data_type <- "Testing"

ModelCompareAUC_Train_250$NofFeatAgn <- "250_Pairs"
ModelCompareAUC_Test_250$NofFeatAgn <- "250_Pairs"

#############################################################################
## Save for the main figure
ModelCompare_RF <- rbind(ModelCompareAUC_Train_250, ModelCompareAUC_Test_250)
ModelCompare_RF$algorithm <- "RF"
save(ModelCompare_RF, file = "./Objs/RF/ModelCompare_RF_AgnPairs_new.rda")

#####################################################################
############################################################
# Load the AUC comparisons from the indivdial genes and combine them with pairs
load("./Objs/RF/ModelCompareAUC_500_new.rda")

# Combine
ModelCompareAUC_Train_250_500 <- rbind(ModelCompareAUC_Train_250, ModelCompareAUC_Train_500)
ModelCompareAUC_Test_250_500 <- rbind(ModelCompareAUC_Test_250, ModelCompareAUC_Test_500)

####################################################################################
####################################################################################
####################################################################################
# Combine all together in one dataframe

ModelCompare_RF_DiffNoFeat <- rbind(ModelCompareAUC_Train_25_50,
                                    ModelCompareAUC_Test_25_50,
                                    ModelCompareAUC_Train_50_100,
                                    ModelCompareAUC_Test_50_100,
                                    ModelCompareAUC_Train_100_200,
                                    ModelCompareAUC_Test_100_200,
                                    ModelCompareAUC_Train_250_500,
                                    ModelCompareAUC_Test_250_500
)

save(ModelCompare_RF_DiffNoFeat, file = "./Objs/RF/ModelCompare_RF_DiffNoFeat_new.rda")

