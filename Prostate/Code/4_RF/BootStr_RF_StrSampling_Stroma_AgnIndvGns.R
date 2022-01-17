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

## Load data
load("./Objs/KTSP/KTSP_STATs_Mechanistic_Stroma.rda")
load("./Objs/MetastasisDataGood.rda")


usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


predictor_data_Train_Mech <- t(KTSP_STATs_Train_Mechanistic)
predictor_data_Test_Mech <- t(KTSP_STATs_Test_Mechanistic)

DataMech_Train <- cbind(predictor_data_Train_Mech, usedTrainGroup)
DataMech_Train <- as.data.frame(DataMech_Train)
DataMech_Train$usedTrainGroup <- as.factor(DataMech_Train$usedTrainGroup)
levels(DataMech_Train[, "usedTrainGroup"]) <- c("No_Mets", "Mets")

names(DataMech_Train) <- make.names(names(DataMech_Train))

colnames(predictor_data_Train_Mech) <- make.names(colnames(predictor_data_Train_Mech))
colnames(predictor_data_Test_Mech) <- make.names(colnames(predictor_data_Test_Mech))

# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  #RF <- randomForest(usedTrainGroup~., data=d, importance = TRUE, ntree = 500, mtry = 10, proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 5, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Mech <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Mech <- Importance_Mech[order(Importance_Mech, decreasing = TRUE), ]
  Importance_Mech <- Importance_Mech[Importance_Mech > 0]
  N_ImportanVariables <- length(Importance_Mech)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Mech, type = "vote")
  ROCTrainMech <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestMech <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainMech$auc, ROCTestMech$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectMech <- boot(data= DataMech_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 12) 

AUCs_RF_Mech <- bootobjectMech$t
colnames(AUCs_RF_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

###################################################
## RF
# Agnostic
load("./Objs/MetastasisDataGood.rda")
load("./Objs/Correlation/RGenes.rda")

### Normalization
usedTrainMat <- normalizeBetweenArrays(mixTrainMat, method = "quantile")[RGenes, ]
usedTestMat <- normalizeBetweenArrays(mixTestMat, method = "quantile")[RGenes, ]

# Groups
usedTrainGroup <- mixTrainGroup
usedTestGroup <- mixTestGroup


names(usedTrainGroup) <- colnames(usedTrainMat)
all(names(usedTrainGroup) == colnames(usedTrainMat))

names(usedTestGroup) <- colnames(usedTestMat)
all(names(usedTestGroup) ==colnames(usedTestMat))

#########
## Detect Top DE genes
TopDEgenes <- SWAP.Filter.Wilcoxon(phenoGroup = usedTrainGroup, inputMat = usedTrainMat, featureNo = nrow(KTSP_STATs_Test_Mechanistic)*2)

## Subset the expression matrix to the top DE genes only
usedTrainMat <- usedTrainMat[TopDEgenes, ]
usedTestMat <- usedTestMat[TopDEgenes, ]


########################


predictor_data_Train_Agnostic <- t(usedTrainMat)
predictor_data_Test_Agnostic <- t(usedTestMat)

DataAgnostic_Train <- cbind(predictor_data_Train_Agnostic, usedTrainGroup)
DataAgnostic_Train <- as.data.frame(DataAgnostic_Train)
DataAgnostic_Train$usedTrainGroup <- as.factor(DataAgnostic_Train$usedTrainGroup)
levels(DataAgnostic_Train[, "usedTrainGroup"]) <- c("No_Mets", "Mets")


# The function for bootstraping
RF_Strap <- function(data, indices) {
  d <- data[indices, ] # allows boot to select sample
  # Select the minimum sample size
  tmp <- as.vector(table(d$usedTrainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  #RF <- randomForest(usedTrainGroup~., data=d, importance = TRUE, ntree = 500, mtry = 10 ,proximity=TRUE, na.action = na.omit, sampsize = sampsizes)
  RF <- tuneRF(x = d[,!colnames(d) == "usedTrainGroup"], y = d$usedTrainGroup, mtryStart = 5, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, sampsize = sampsizes)
  Importance_Agnostic <- randomForest::importance(RF, scale=FALSE, type = 2)
  Importance_Agnostic <- Importance_Agnostic[order(Importance_Agnostic, decreasing = TRUE), ]
  Importance_Agnostic <- Importance_Agnostic[Importance_Agnostic > 0]
  N_ImportanVariables <- length(Importance_Agnostic)
  PhenoTrain <- d$usedTrainGroup
  PredictorTrainData <- d
  PredictorTrainData$usedTrainGroup <- NULL
  train_preds <- predict(RF, newdata = PredictorTrainData, type = "vote")
  test_preds <- predict(RF, newdata = predictor_data_Test_Agnostic, type = "vote")
  ROCTrainAgnostic <- roc(PhenoTrain, train_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  ROCTestAgnostic <- roc(usedTestGroup, test_preds[,2], plot = F, print.auc = TRUE, levels = c("No_Mets", "Mets"), direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  return(c(ROCTrainAgnostic$auc, ROCTestAgnostic$auc, N_ImportanVariables))
}


set.seed(333)
bootobjectAgnostic <- boot(data= DataAgnostic_Train, statistic= RF_Strap, R= 1000, parallel = "multicore", ncpus = 12) 

AUCs_RF_Agnostic <- bootobjectAgnostic$t
colnames(AUCs_RF_Agnostic) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")


save(bootobjectAgnostic, bootobjectMech, file = "./Objs/RF/RFBootObjects_Stroma_AgnIndvGns.rda")

load("./Objs/RF/RFBootObjects_Lipogenesis_AgnIndvGns.rda")

AUCs_RF_Agnostic <- bootobjectAgnostic$t
colnames(AUCs_RF_Agnostic) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

AUCs_RF_Mech <- bootobjectMech$t
colnames(AUCs_RF_Mech) <- c("AUC_Train", "AUC_Test", "N_ImportanVariables")

## Calculate the difference and CI of the difference in the training data
DiffAgnostic <- AUCs_RF_Agnostic[, "AUC_Train"] - AUCs_RF_Agnostic[, "AUC_Test"]
quantile(DiffAgnostic, c(0.025, 0.975))
#colnames(Diff) <- "Diff"

## Calculate the difference and CI of the difference between mechanistic training and mechanistic testing
DiffMech <- AUCs_RF_Mech[, "AUC_Train"] - AUCs_RF_Mech[, "AUC_Test"]
quantile(DiffMech, c(0.025, 0.975))


## Plot the distributions of the AUCs from both methods in the training data
MechanisticAUC_Train <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Train"])
AgnosticAUC_Train <- data.frame(AUC = AUCs_RF_Agnostic[, "AUC_Train"])

MechanisticAUC_Train$modelType <- "Mech"
AgnosticAUC_Train$modelType <- "Agnostic"

ModelCompareAUC_Train <- rbind(MechanisticAUC_Train, AgnosticAUC_Train)

## Plot the distributions of the AUCs from both methods in the testing data
MechanisticAUC_Test <- data.frame(AUC = AUCs_RF_Mech[, "AUC_Test"])
AgnosticAUC_Test <- data.frame(AUC = AUCs_RF_Agnostic[, "AUC_Test"])

MechanisticAUC_Test$modelType <- "Mech"
AgnosticAUC_Test$modelType <- "Agnostic"

ModelCompareAUC_Test <- rbind(MechanisticAUC_Test, AgnosticAUC_Test)

##########
## The N of important features (output)
ImportanceAgnostic <- data.frame(NofOutputFeatures = AUCs_RF_Agnostic[, "N_ImportanVariables"])
ImportanceMech <- data.frame(NofOutputFeatures = AUCs_RF_Mech[, "N_ImportanVariables"])

ImportanceAgnostic$modelType <- "Agnostic"
ImportanceMech$modelType <- "Mech"

ModelCompare_ImportantFeatures <- rbind(ImportanceAgnostic, ImportanceMech)

############
## Plots
My_Theme = theme(
  axis.title.x = element_text(size = 7),
  axis.text.x = element_text(size = 5),
  axis.title.y = element_text(size = 7),
  axis.text.y = element_text(size = 5),
  plot.title = element_text(size=9)
)

DiffHist_Agnostic <- ggplot(as.data.frame(DiffAgnostic), aes(DiffAgnostic, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0.10, 0.50)) +
  labs(title="Histogram of the difference between the training and testing data using the agnostic model") + My_Theme

DiffHist_Mech <- ggplot(as.data.frame(DiffMech), aes(DiffMech, fill = "red")) + 
  geom_histogram(bins = 25) +
  scale_x_continuous(limits = c(0.10, 0.50)) +
  labs(title="Histogram of the difference between the training and testing data using the mechanistic model") + My_Theme

AUC_Train_DistrHist <- ggplot(ModelCompareAUC_Train, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5, adjust=0.1) +
  scale_x_continuous(limits = c(0.6, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the training data") + My_Theme

AUC_Test_DistrHist <- ggplot(ModelCompareAUC_Test, aes(AUC, fill = modelType)) + 
  geom_density(alpha = 0.5) +
  scale_x_continuous(limits = c(0.6, 1)) +
  labs(title="AUC distribution of the agnostic and mechanistic RF models in the testing data") + My_Theme

BarPlotImportance <- ggplot(data=ModelCompare_ImportantFeatures, aes(x=NofOutputFeatures, color = modelType, fill = modelType)) +
  geom_bar(stat="bin") +
  scale_x_continuous(limits = c(90, 105)) +
  labs(title="Distribution of the number of important features (output)") + 
  facet_grid(~modelType, scale='free_x') +
  My_Theme

png("./Figs/RF/RF_BS_Stroma100_AgnIndvGns.png", width = 3000, height = 1500, res = 300)
(AUC_Train_DistrHist / AUC_Test_DistrHist) | (DiffHist_Agnostic / DiffHist_Mech) + plot_layout(widths = c(2, 1)) 
dev.off()


