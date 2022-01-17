
###############################################################################################
## Mohamed Omar
## 21/07/2019
## Goal: Discovery and validation of a small gene signature that can predict the metastasis potential in primary prostate cancer
################################################################################################

## Clean work space
rm(list = ls())

## Set the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")


## Load necessary packages
library(MetaIntegrator)
library(GEOquery)
library(pROC)
library(caret)
library(genefilter)
library(mltools)

########################################################
## Load Training data sets (4)
load("./Data/ProstateData.rda")

## Load Testing data set
## Load Testing data set
load("./Data/Dataset1.rda")
load("./Data/Dataset6.rda")


#######################################################

## Getting the phenotype data for each data set
#pheno4 <- pData(Dataset4)
#ProstateData$originalData$GSE46691$pheno <- pheno4
pheno2 <- Dataset2$pheno
pheno3 <- Dataset3$pheno
pheno4 <- Dataset4$pheno
Pheno5 <- Dataset5$pheno
Pheno6 <- Dataset6$pheno

################################

## Getting the expression data for each data set
## load expr4
#load("/Users/mohamedomar/Documents/Research/Projects/Prostate/Data/expr4.rda")
#ProstateData$originalData$GSE46691$expr <- expr4
expr1 <- Dataset1$expr
expr2 <- Dataset2$expr
expr3 <- Dataset3$expr
expr4 <- Dataset4$expr
expr5 <- Dataset5$expr
expr6 <- Dataset6$expr

## Checking if the expression data are normalized and log2 transformed
# boxplot(expr2[,1:15], outline= FALSE)
# boxplot(expr3[,1:15], outline = FALSE)
# boxplot(expr4[,1:15], outline= FALSE)

#################################################################
## Create a list containing training data sets
AllDataSets <- list(Dataset2, Dataset3, Dataset4, Dataset5)
names(AllDataSets) <- c("GSE55935", "GSE51066", "GSE46691", "GSE41408")

###################################################################

## Annotate expression

## Expr1
head(rownames(expr1))
rownames(expr1) <- Dataset1$keys
expr1 <- expr1[!is.na(rownames(expr1)), ]
#####################
## expr2
#head(rownames(expr2))
#rownames(expr2) <- Dataset2$keys
expr2 <- expr2[!is.na(rownames(expr2)), ]
dim(expr2)

# X2 <- expr2
# ffun <- filterfun(pOverA(p = 0.8, A = 100))
# filt2 <- genefilter(2^X2, ffun)
# expr2 <- expr2[filt2, ]
#####################
## expr3
#rownames(expr3) <- Dataset3$keys
#expr3 <- expr3[!is.na(rownames(expr3)), ]
dim(expr3)
# X3 <- expr3
# filt3 <- genefilter(2^X3, ffun)
# expr3 <- expr3[filt3, ]
# #######################
# dim(expr4)
# X4 <- expr4
# filt4 <- genefilter(2^X4, ffun)
# expr4 <- expr4[filt4, ]

######################
## expr5
head(rownames(expr5))

rownames(expr5) <- Dataset5$keys
expr5 <- expr5[!is.na(rownames(expr5)), ]
dim(expr5)

#####################
## expr6
# processed

# ############################################################
####################################################################
## Modify pheno2
# Remove cell lines
pheno2 <- pheno2[-c(47:54), ]
pheno2$Metastasis <- pheno2$`lymph node metastasis status:ch1`
pheno2$Metastasis[pheno2$Metastasis == 0] <- "No_Mets"
pheno2$Metastasis[pheno2$Metastasis == 1] <- "Mets"
pheno2 <- pheno2[!(pheno2$Metastasis == "NA"), ]

pheno2$Metastasis <- as.factor(pheno2$Metastasis)
table(pheno2$Metastasis)
# Modify expr2
expr2 <- expr2[,colnames(expr2) %in% rownames(pheno2)]
all(rownames(pheno2) == colnames(expr2))
## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset2$expr <- expr2
Dataset2$pheno <- pheno2

####### 
## Modify pheno4
pheno4$Metastasis <- pheno4$`metastatic event:ch1`
pheno4$Metastasis[pheno4$Metastasis == 0] <- "No_Mets"
pheno4$Metastasis[pheno4$Metastasis == 1] <- "Mets"
pheno4 <- pheno4[!(pheno4$Metastasis == "NA"), ]

pheno4$Metastasis <- as.factor(pheno4$Metastasis)
table(pheno4$Metastasis)
levels(pheno4$Metastasis)

## Modify sample names to match sample names of pheno4
head(colnames(expr4))
colnames(expr4) <- gsub(".CEL", "", colnames(expr4))
colnames(expr4) <- gsub(".+\\.", "", colnames(expr4))
expr4 <- expr4[,order(colnames(expr4))]

rownames(pheno4) <- pheno4$description
head(rownames(pheno4))
rownames(pheno4) <- gsub(".+\\.", "", rownames(pheno4))
pheno4 <- pheno4[order(rownames(pheno4)), ]

all(rownames(pheno4) == colnames(expr4))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset4$expr <- expr4
Dataset4$pheno <- pheno4

######
# Modify pheno3
pheno3$Metastasis <- pheno3$`metastatic event:ch1`
pheno3$Metastasis[pheno3$Metastasis == 0] <- "No_Mets"
pheno3$Metastasis[pheno3$Metastasis == 1] <- "Mets"

pheno3$Metastasis <- as.factor(pheno3$Metastasis)
table(pheno3$Metastasis)
levels(pheno3$Metastasis)

all(rownames(pheno3) == colnames(expr3))

## Finally, replace the expression and phenotype data in the dataset with the new modified versions
Dataset3$pheno <- pheno3
Dataset3$expr <- expr3

#####
## Modify Pheno5
table(Pheno5$`development of metastasis:ch1`)

Pheno5$Metastasis <- Pheno5$`development of metastasis:ch1`

Pheno5$Metastasis[Pheno5$Metastasis == "no"] <- "No_Mets"
Pheno5$Metastasis[Pheno5$Metastasis == "yes"] <- "Mets"

Pheno5$Metastasis <- as.factor(Pheno5$Metastasis)
table(Pheno5$Metastasis)
all(rownames(Pheno5) == colnames(expr5))

Dataset5$expr <- expr5
Dataset5$pheno <- Pheno5

######
# Pheno6
table(Pheno6$`extreme:ch1`)
Pheno6$Metastasis <- Pheno6$`extreme:ch1`

Pheno6$Metastasis[Pheno6$Metastasis == "Indolent"] <- "No_Mets"
Pheno6$Metastasis[Pheno6$Metastasis == "Lethal"] <- "Mets"

Pheno6$Metastasis <- as.factor(Pheno6$Metastasis)
table(Pheno6$Metastasis)
all(rownames(Pheno6) == colnames(expr6))

Dataset6$expr <- expr6
Dataset6$pheno <- Pheno6


#########################################################################
############################################################################

## Label samples (All samples need to be assigned labels in the $class vector, 1 for ‘disease’ or 0 for ‘control’)
Dataset2 <- classFunction(Dataset2, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset3 <- classFunction(Dataset3, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset4 <- classFunction(Dataset4, column = "Metastasis", diseaseTerms = c("Mets"))
Dataset5 <- classFunction(Dataset5, column = "Metastasis", diseaseTerms = c("Mets"))

############################################################################
############################################################################
#########################################################################
## The metaanalysis

## Creating the meta object
AllDataSets <- list(Dataset2, Dataset3, Dataset4, Dataset5)
names(AllDataSets) <- c("GSE55935", "GSE51066", "GSE46691", "GSE41408")

Prostate_meta <- list()
Prostate_meta$originalData <- AllDataSets

## Replace keys within each data set
Prostate_meta$originalData$GSE55935$keys <- rownames(expr2)
Prostate_meta$originalData$GSE51066$keys <- rownames(expr3)
Prostate_meta$originalData$GSE46691$keys <- rownames(expr4)
Prostate_meta$originalData$GSE41408$keys <- rownames(expr5)

## Check the meta object before the metaanalysis
checkDataObject(Prostate_meta, "Meta", "Pre-Analysis") ## If true, Proceed to the meta analysis


## Run the meta analysis
Prostate_metaanalysis <- runMetaAnalysis(Prostate_meta, runLeaveOneOutAnalysis = F, maxCores = 3)

## Filter out significant genes from the metaanalysis results (this will be the gene signature that separates Metas from No_Mets)
Prostate_metaanalysis <- filterGenes(Prostate_metaanalysis, isLeaveOneOut = F, effectSizeThresh = 0, FDRThresh = 0.15, numberStudiesThresh = 4, heterogeneityPvalThresh = 0.05)



## Stats of all genes
AllGenes <- Prostate_metaanalysis$metaAnalysis$pooledResults

## The OncoDx genes (15 genes)
OncoDx <- c("BGN", "COL1A1", "SFRP4", "FLNC", "GSN", "TPM2", "GSTM2", "TPX2", "LAMB3", "FAM13C", "KLK2", "AZGP1", "SRD5A2", "DUSP1", "FOS")

OncoDx_Genes <- c("BGN", "COL1A1", "SFRP4", "FLNC", "GSN", "TPM2", "GSTM2", "TPX2", "LAMB3", "FAM13C", "KLK2", "AZGP1", "SRD5A2", "DUSP1", "FOS")
save(OncoDx_Genes, file = "./Objs/OncoDxGenes.rda")

# Remocve reference genes
#OncoDx <- c("BGN", "COL1A1", "SFRP4", "FLNC", "GSN", "TPM2", "GSTM2", "TPX2", "FAM13C", "KLK2", "AZGP1", "SRD5A2")

summary(OncoDx %in% rownames(AllGenes))

# Keep only the stats of OncoDx genes (2 genes are absent)
OncoDx <- AllGenes[OncoDx, ]
#OncoDx <- na.omit(OncoDx)

# Replace the gene names with the OncoDx genes
Prostate_metaanalysis$filterResults$FDR0.15_es0_nStudies4_looaFALSE_hetero0.05$posGeneNames <- c("BGN", "COL1A1", "SFRP4", "TPX2")
Prostate_metaanalysis$filterResults$FDR0.15_es0_nStudies4_looaFALSE_hetero0.05$negGeneNames <- c("FLNC", "GSN", "TPM2", "GSTM2", "LAMB3", "FAM13C", "KLK2", "AZGP1", "SRD5A2", "DUSP1", "FOS")

OncoDx_Filter <- Prostate_metaanalysis$filterResults[[1]]
OncoDx_Filter_Summary <- summarizeFilterResults(metaObject = Prostate_metaanalysis, getMostRecentFilter(Prostate_metaanalysis))


## Create a summary ROC curve (Training data sets)
set.seed(333)
summaryROCPlot(metaObject = Prostate_metaanalysis, filterObject = OncoDx_Filter, bootstrapReps = 500)
###########################################################################################
#############################################################################

## The next step is validation on an indepndent data set (GSE116918)

# Get the validation data set
#valDataSet <- getGEOData("GSE116918")
valDataSet1 <- Dataset1
valDataSet2 <- Dataset6

# Acess the phenotype and expression data
val_pheno1 <- Dataset1$pheno
val_expr1 <- Dataset1$expr

val_pheno2 <- Pheno6
val_expr2 <- expr6

## Modify val_expr1
rownames(val_expr1) <- valDataSet1$keys
summary(is.na(rownames(val_expr1)))
val_expr1 <- val_expr1[!is.na(rownames(val_expr1)), ]
dim(val_expr1)


## Modify val_pheno1
val_pheno1$Metastasis <- val_pheno1$`met event (1=yes, 0=no):ch1`
val_pheno1$Metastasis[val_pheno1$Metastasis == 0] <- "No_Mets"
val_pheno1$Metastasis[val_pheno1$Metastasis == 1] <- "Mets"
val_pheno1$Metastasis <- as.factor(val_pheno1$Metastasis)
table(val_pheno1$Metastasis)
all(rownames(val_pheno1) == colnames(val_expr1))

valDataSet1$pheno <- val_pheno1
valDataSet1$expr <- val_expr1
valDataSet1$keys <- rownames(val_expr1)

## Label the samples
valDataSet1 <- classFunction(valDataSet1, column = "Metastasis", diseaseTerms = c("Mets"))
valDataSet2 <- classFunction(valDataSet2, column = "Metastasis", diseaseTerms = c("Mets"))

##############################################
## Examine the performance of OncoDx filter

## ROC plot
rocPlot(datasetObject = valDataSet1, filterObject = OncoDx_Filter)

# Validation dataset2
valDataSet2$keys <- rownames(val_expr2)
rocPlot(datasetObject = valDataSet2, filterObject = OncoDx_Filter)

# PRC plot
prcPlot(datasetObject = valDataSet1, filterObject = OncoDx_Filter)

# val dataset2
prcPlot(datasetObject = valDataSet2, filterObject =OncoDx_Filter)

# Violin plot
violinPlot(filterObject = OncoDx_Filter, datasetObject = valDataSet1, labelColumn = "Metastasis")

violinPlot(filterObject = OncoDx_Filter, datasetObject = valDataSet2, labelColumn = "Metastasis")

##############
## Calculate a signature score (Z score) and add it to the phenotype table
val_pheno1$score <- calculateScore(filterObject = OncoDx_Filter, datasetObject = valDataSet1)

# Val Dataset2
val_pheno2$score <- calculateScore(filterObject = OncoDx_Filter, datasetObject = valDataSet2)

## Find the best threshold (for further use for the classification in the Test data)
thr_test1 <- coords(roc(val_pheno1$Metastasis, val_pheno1$score, levels = c("No_Mets", "Mets")),"best", transpose = TRUE)["threshold"]
thr_test1

# Val Dataset2
thr_test2 <- coords(roc(val_pheno2$Metastasis, val_pheno2$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
thr_test2 


## Find the optimal trade off between the sensitivity and specificity (We want good sensitivity with preservation of a decent accuracy)
coords(roc(val_pheno1$Metastasis, val_pheno1$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")
# Val Dataset2
coords(roc(val_pheno2$Metastasis, val_pheno2$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")

## Predictions
# Val Dataset1
Test_Predictions1 <- ifelse(val_pheno1$score >= thr_test1, "Mets", "No_Mets")
confusionMatrix(as.factor(Test_Predictions1), val_pheno1$Metastasis, positive = "Mets")

mcc(preds = as.factor(Test_Predictions1), actuals = val_pheno1$Metastasis)

# Val Dataset2
Test_Predictions2 <- ifelse(val_pheno2$score >= thr_test2, "Mets", "No_Mets")
confusionMatrix(as.factor(Test_Predictions2), val_pheno2$Metastasis, positive = "Mets")

mcc(preds = as.factor(Test_Predictions2), actuals = val_pheno2$Metastasis)
############################################


# Survival (Meta-score)

## Load necessary packages
library(survival)
library(survminer)

## Turn the score column into numeric
val_pheno1[,"score"] <- as.numeric(val_pheno1[,"score"]) 
val_pheno2[,"score"] <- as.numeric(val_pheno2[, "score"])

## Divide the dataset into quantiles based on the Meta score risk.
#quantiles <- quantile(val_pheno1$score, probs=c(0.33333,0.66667))
# val_Pheno1
val_pheno1$Meta_Score <- val_pheno1$score
val_pheno1[which(val_pheno1$score <= thr_test1), "Meta_Score"] = "low"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
val_pheno1[which(val_pheno1$score > thr_test1),"Meta_Score"] = "high"

# val_pheno2
val_pheno2$Meta_Score <- val_pheno2$score
val_pheno2[which(val_pheno2$score <= thr_test2), "Meta_Score"] = "low"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
val_pheno2[which(val_pheno2$score > thr_test2),"Meta_Score"] = "high"


val_pheno1$Time <- as.numeric(val_pheno1$`follow-up time (met, months):ch1`)
val_pheno2$Time <- as.numeric(val_pheno2$`fup.month:ch1`)

val_pheno1$Event <- as.numeric(val_pheno1$`met event (1=yes, 0=no):ch1`)
val_pheno2$Event <- as.numeric(val_pheno2$Metastasis)
val_pheno2$Event[val_pheno2$Event == 2] <- 0 

## Keep only relevant information from the phenotype table
surv_data1 <- val_pheno1[,c("Time","Event","score","Meta_Score","gleason grade:ch1", "patient age (years):ch1", "psa (ng/ml):ch1", "t-stage:ch1")]
surv_data1$Meta_Score <- factor(surv_data1$Meta_Score, levels = c("low", "high"))

surv_data2 <- val_pheno2[, c("Time", "Event", "score", "Meta_Score", "gleason:ch1", "age:ch1", "cancer.percent:ch1")]
surv_data2$Meta_Score <- factor(surv_data2$Meta_Score, levels = c("low", "high"))

# Create a survival object
surv_data.surv1 <-  with(surv_data1, Surv(Time, Event == 1))
surv_data.surv2 <-  with(surv_data2, Surv(Time, Event == 1))

#Calculate p-value
survdifftest1 <- survdiff(surv_data.surv1 ~ Meta_Score, data = surv_data1)
survpvalue1 <- 1 - pchisq(survdifftest1$chisq, length(survdifftest1$n) - 1)

survdifftest2 <- survdiff(surv_data.surv2 ~ Meta_Score, data = surv_data2)
survpvalue2 <- 1 - pchisq(survdifftest2$chisq, length(survdifftest2$n) - 1)

## Create a linear test p-value
# surv_data_lin <- val_pheno1[,c("Time","Event","Meta_Score")]
# surv_data_lin[,"Meta_Score"] = as.vector(surv_data_lin[,"Meta_Score"])
# surv_data_lin[which(surv_data_lin[,"Meta_Score"]=="low"),"Meta_Score"] = 1
# surv_data_lin[which(surv_data_lin[,"Meta_Score"]=="int"),"Meta_Score"] = 2
# surv_data_lin[which(surv_data_lin[,"Meta_Score"]=="high"),"Meta_Score"] = 3
# surv_data_lin[ , "Meta_Score"] <- as.numeric(surv_data_lin[ , "Meta_Score"])
# survpvalue_linear <- summary(coxph(Surv(Time, as.numeric(Event))~Meta_Score, data=surv_data_lin))$sctest[3]
# survpvalue_linear <-format(as.numeric(survpvalue_linear), digits=3)


## Kaplan-Meier curve
png(filename = "./Figs/Meta/KM_Survival.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data1),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the gene signature Meta score in the Test data", xlab = "Time in months", ylab = "Survival")

dev.off()

###############

png(filename = "./Figs/Meta/KM_Survival2.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data2),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the gene signature Meta score in the Test data", xlab = "Time in months", ylab = "Survival")

dev.off()

#############################################################################
###### COX 

## First validation dataset
Cox_data1 <- surv_data1
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data1$Age <- ifelse(Cox_data1$`patient age (years):ch1` >= 70, ">=70 years", "<70 years")
Cox_data1$Age <- factor(Cox_data1$Age, levels = c("<70 years", ">=70 years"))

Cox_data1$PSA <- ifelse(Cox_data1$`psa (ng/ml):ch1` > 20, "High PSA (> 20)", "Low PSA (<= 20)")
Cox_data1$PSA <- factor(Cox_data1$PSA, levels = c("Low PSA (<= 20)", "High PSA (> 20)"))

Cox_data1$Gleason <- ifelse(Cox_data1$`gleason grade:ch1` == 6 , "Gleason = 6", "Gleason > 6")
Cox_data1$Gleason <- factor(Cox_data1$Gleason, levels = c("Gleason = 6", "Gleason > 6"))

Cox_data1$T_Stage <- ifelse(Cox_data1$`t-stage:ch1` == c("T1", "T2"), "T1/T2", "T3/T4/Unknown")
Cox_data1$T_Stage <- factor(Cox_data1$T_Stage, levels = c("T1/T2", "T3/T4/Unknown"))

## Univariate Cox model
# PFS
covariates1 <- c("Meta_Score", "Age", "PSA", "Gleason", "T_Stage")
univ_formulas1 <- sapply(covariates1,
                         function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models1 <- lapply( univ_formulas1, function(x){coxph(x, data = Cox_data1)})
# Extract data 
univ_results1 <- lapply(univ_models1,
                        function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                        })

UniCOX1 <- t(as.data.frame(univ_results1, check.names = F))
as.data.frame(UniCOX1)
########################################

## Multivariate

CoxModel1 <- coxph(Surv(Time, Event)~ Meta_Score+Age+PSA+Gleason+T_Stage, data=Cox_data1)

png(filename = "./Figs/Meta/CoxModel.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel1, data = Cox_data1, main = "Cox proportional hazards model")
dev.off()

############################################################################
############################################################################
## Second validation dataset
Cox_data2 <- surv_data2
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data2$Age <- ifelse(Cox_data2$`age:ch1` >= 70, ">=70 years", "<70 years")
Cox_data2$Age <- factor(Cox_data2$Age, levels = c("<70 years", ">=70 years"))

#Cox_data1$PSA <- ifelse(Cox_data1$`psa (ng/ml):ch1` > 20, "High PSA (> 20)", "Low PSA (<= 20)")
#Cox_data1$PSA <- factor(Cox_data1$PSA, levels = c("Low PSA (<= 20)", "High PSA (> 20)"))

Cox_data2$Gleason <- ifelse(Cox_data2$`gleason:ch1` == 6 , "Gleason = 6", "Gleason > 6")
Cox_data2$Gleason <- factor(Cox_data2$Gleason, levels = c("Gleason = 6", "Gleason > 6"))

Cox_data2$`cancer.percent:ch1` <- as.integer(Cox_data2$`cancer.percent:ch1`)
Cox_data2$CancerPercent <- ifelse(Cox_data2$`cancer.percent:ch1` >= 50, ">= 50", "< 50")
Cox_data2$CancerPercent <- factor(Cox_data2$CancerPercent, levels = c("< 50", ">= 50"))

## Univariate Cox model
# PFS
covariates2 <- c("Meta_Score", "Age", "Gleason", "CancerPercent")
univ_formulas2 <- sapply(covariates2,
                         function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models2 <- lapply( univ_formulas2, function(x){coxph(x, data = Cox_data2)})
# Extract data 
univ_results2 <- lapply(univ_models2,
                        function(x){ 
                          x <- summary(x)
                          p.value<-signif(x$wald["pvalue"], digits=2)
                          wald.test<-signif(x$wald["test"], digits=2)
                          beta<-signif(x$coef[1], digits=2);#coeficient beta
                          HR <-signif(x$coef[2], digits=2);#exp(beta)
                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                          HR <- paste0(HR, " (", 
                                       HR.confint.lower, "-", HR.confint.upper, ")")
                          res<-c(beta, HR, wald.test, p.value)
                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                        "p.value")
                          return(res)
                          #return(exp(cbind(coef(x),confint(x))))
                        })

UniCOX2 <- t(as.data.frame(univ_results2, check.names = F))
as.data.frame(UniCOX2)
########################################

## Multivariate

CoxModel2 <- coxph(Surv(Time, Event)~ Meta_Score+Age+Gleason+CancerPercent, data=Cox_data2)

png(filename = "./Figs/Meta/CoxModel2.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel2, data = Cox_data2, main = "Cox proportional hazards model")
dev.off()
