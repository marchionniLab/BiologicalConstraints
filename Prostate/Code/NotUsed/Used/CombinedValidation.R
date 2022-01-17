## The next step is validation on indepndent data set (GSE116918, and GSE16560)

# Get the validation data set
#valDataSet <- getGEOData("GSE116918")
valDataSet1 <- Dataset1
valDataSet2 <- Dataset6
valDataSet3 <- Dataset7

# Acess the phenotype and expression data
val_pheno1 <- Dataset1$pheno
val_expr1 <- Dataset1$expr

val_pheno2 <- Pheno6
val_expr2 <- expr6

val_pheno3 <- Pheno7
val_expr3 <- expr7



#####


#################################
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
## val_Pheno2 is already processed

## Label the samples
#valDataSet1 <- classFunction(valDataSet1, column = "Metastasis", diseaseTerms = c("Mets"))
#valDataSet2 <- classFunction(valDataSet2, column = "Metastasis", diseaseTerms = c("Mets"))
##############################################
allExpr <- list(val_expr1, val_expr2, val_expr3)


CommonGns <- Reduce("intersect", lapply(allExpr, rownames))

val_expr1 <- val_expr1[CommonGns, ]
val_expr2 <- val_expr2[CommonGns, ]
val_expr3 <- val_expr3[CommonGns, ]

valExpr <- cbind(val_expr1, val_expr2, val_expr3) 
valExpr <- normalizeBetweenArrays(valExpr, method = "quantile")

### Now lets examine the performance of our filter on the independent data set
#val_pheno1 <- val_pheno1[, c(36:45)]
#val_pheno2 <- val_pheno2[, c(47:56)]


val_pheno1 <- val_pheno1[, c("gleason grade:ch1", "patient age (years):ch1", "follow-up time (met, months):ch1", "Metastasis")]
colnames(val_pheno1) <- c("Gleason", "Age", "Time", "Metastasis")
val_pheno1$Age <- NULL

val_pheno2 <- val_pheno2[, c("gleason:ch1", "age:ch1", "fup.month:ch1", "Metastasis")]
colnames(val_pheno2) <- c("Gleason", "Age", "Time", "Metastasis")
val_pheno2$Age <- NULL

val_pheno3 <- val_pheno3[, c("tumour gleason:ch1", "total follow up (months):ch1", "Metastasis")]
colnames(val_pheno3) <- c("Gleason", "Time", "Metastasis")
val_pheno3$Gleason <- gsub("\\=.+", "", val_pheno3$Gleason)

#### Combine pheno
val_pheno <- rbind(val_pheno1, val_pheno2, val_pheno3)

## Check for consistent sample names
all(rownames(val_pheno) == colnames(valExpr))

###########
## Collect into one data set
valDataSet <- list()
valDataSet$pheno <- val_pheno
valDataSet$expr <- valExpr
valDataSet$keys <- rownames(valExpr)
valDataSet$formattedName <- c("Validation dataset")

## Class
valDataSet <- classFunction(valDataSet, column = "Metastasis", diseaseTerms = c("Mets"))


# Validation dataset
rocPlot(datasetObject = valDataSet, filterObject = New_filter)


################################

#### Using PRC plot

# val dataset
prcPlot(datasetObject = valDataSet, filterObject = New_filter)

##############
#### Calculate a signature score (Z score) and add it to the phenotype table

# Val Dataset1
val_pheno$score <- calculateScore(filterObject = New_filter, datasetObject = valDataSet)

#########################################

## ROC curve in the Testing data using pROC

roc(val_pheno$Metastasis, val_pheno$score, plot = TRUE, print.auc=TRUE, levels = c("No_Mets", "Mets"), col="blue", lwd=2, grid=TRUE, auc = TRUE, ci = TRUE)

#### Find the best threshold (for further use for the classification in the Test data)

# Val Dataset1
thr_test <- coords(roc(val_pheno$Metastasis, val_pheno$score, levels = c("No_Mets", "Mets"),),"best", transpose = TRUE)["threshold"]
thr_test 



#### Find the optimal trade off between the sensitivity and specificity (We want good sensitivity with preservation of a decent accuracy)

# Val Dataset
coords(roc(val_pheno$Metastasis, val_pheno$score, levels = c("No_Mets", "Mets")), transpose = TRUE, "local maximas")

#############################
##### Predictions
Test_Predictions <- ifelse(val_pheno$score <= thr_test, "Mets", "No_Mets")
confusionMatrix(as.factor(Test_Predictions), val_pheno$Metastasis, positive = "Mets", mode = "everything")

mcc(preds = as.factor(Test_Predictions), actuals = val_pheno$Metastasis)


#############

Meta <- list() 
Meta$originalData$GSE116918 <- valDataSet1
Meta$originalData$GSE16560 <- valDataSet2
Meta$originalData$GSE70769 <- valDataSet3

Meta$originalData$GSE16560$keys <- rownames(Meta$originalData$GSE16560$expr)

Meta$originalData$GSE116918 <- classFunction(Meta$originalData$GSE116918, column = "Metastasis", diseaseTerms = c("Mets"))
Meta$originalData$GSE16560 <- classFunction(Meta$originalData$GSE16560, column = "Metastasis", diseaseTerms = c("Mets"))
Meta$originalData$GSE70769 <- classFunction(Meta$originalData$GSE70769, column = "Metastasis", diseaseTerms = c("Mets"))

set.seed(333)
summaryROCPlot(metaObject = Meta, filterObject = New_filter)

set.seed(333)
pooledROCPlot(metaObject = Meta, filterObject = New_filter)

##################################################

## Turn the score column into numeric
val_pheno[,"score"] <- as.numeric(val_pheno[,"score"]) 

## Divide the dataset into quantiles based on the Meta score risk.
#quantiles <- quantile(val_pheno1$score, probs=c(0.33333,0.66667))
# val_Pheno1
val_pheno$Meta_Score <- val_pheno$score
val_pheno[which(val_pheno$score <= thr_test), "Meta_Score"] = "high"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
val_pheno[which(val_pheno$score > thr_test),"Meta_Score"] = "low"


val_pheno$Time <- as.numeric(val_pheno$Time)

val_pheno$Event <- as.numeric(val_pheno$Metastasis)
val_pheno$Event[val_pheno$Event == 2] <- 0 

## Keep only relevant information from the phenotype table
surv_data <- val_pheno[,c("Time","Event","score","Meta_Score","Gleason")]
surv_data$Meta_Score <- factor(surv_data$Meta_Score, levels = c("low", "high"))


# Create a survival object
surv_data.surv <-  with(surv_data, Surv(Time, Event == 1))

#Calculate p-value
survdifftest <- survdiff(surv_data.surv ~ Meta_Score, data = surv_data)
survpvalue <- 1 - pchisq(survdifftest$chisq, length(survdifftest$n) - 1)



## Kaplan-Meier curve
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the gene signature Meta score in the Test data", xlab = "Time in months", ylab = "Survival")



#############################################################################
###### COX 

## First validation dataset
Cox_data <- surv_data
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data$Age <- ifelse(Cox_data$Age >= 70, ">=70 years", "<70 years")
Cox_data$Age <- factor(Cox_data$Age, levels = c("<70 years", ">=70 years"))

Cox_data$Gleason <- ifelse(Cox_data$Gleason == 6 , "Gleason = 6", "Gleason > 6")
Cox_data$Gleason <- factor(Cox_data$Gleason, levels = c("Gleason = 6", "Gleason > 6"))

## Univariate Cox model
# PFS
covariates <- c("Meta_Score", "Gleason")
univ_formulas <- sapply(covariates,
                         function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = Cox_data)})
# Extract data 
univ_results <- lapply(univ_models,
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

UniCOX <- t(as.data.frame(univ_results, check.names = F))
as.data.frame(UniCOX)
########################################

## Multivariate

CoxModel <- coxph(Surv(Time, Event)~ Meta_Score+Gleason, data=Cox_data)

png(filename = "./Figs/Meta/CoxModel.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel, data = Cox_data, main = "Cox proportional hazards model")
dev.off()

############################################################################
#############################################################################
##### Survival by gene
PositiveGenes <- PositiveGenes[PositiveGenes %in% rownames(valExpr)]
NegativeGenes <- NegativeGenes[NegativeGenes %in% rownames(valExpr)]

# Keep only the Good AND Bad Genes (Associated with no metastasis // Metastasis)
exprUP <- valExpr[PositiveGenes, ]
exprDown <- valExpr[NegativeGenes, ]

# Z-transform exprGood
exprUP <- t(scale(t(exprUP)))
exprDown <- t(scale(t(exprDown)))

# create a merged pdata and Z-scores object
CoxDataUp <- data.frame(val_pheno, t(exprUP))
CoxDataDown <- data.frame(val_pheno, t(exprDown))


#################################################################################
## Define optimal cutpoints for each gene (converting the absolute expression into categorical low/high expression)

## For the positive/Up genes
CutPoint_UPGenes <- surv_cutpoint(data = CoxDataUp, time = "Time", event = "Event", variables = PositiveGenes)
CutPoint_UPGenes

SurvData_UpGenes <- surv_categorize(CutPoint_UPGenes)

## For the negative/Down genes
CutPoint_DownGenes <- surv_cutpoint(data = CoxDataDown, time = "Time", event = "Event", variables = NegativeGenes)
CutPoint_DownGenes

SurvData_DownGenes <- surv_categorize(CutPoint_DownGenes)


########################################################################  
## Fit Up genes

Fit_TMSB10 <- survfit(Surv(Time, Event) ~ TMSB10 , data = SurvData_UpGenes)
Fit_ASPN <- survfit(Surv(Time, Event) ~ ASPN, data = SurvData_UpGenes)
Fit_P4HA2 <- survfit(Surv(Time, Event) ~ P4HA2, data = SurvData_UpGenes)
Fit_PTDSS1 <- survfit(Surv(Time, Event) ~ PTDSS1, data = SurvData_UpGenes)
Fit_NCOA2 <- survfit(Surv(Time, Event) ~ NCOA2, data = SurvData_UpGenes)



##################################################################
## Fit down genes

Fit_AZGP1 <- survfit(Surv(Time, Event) ~ AZGP1, data = SurvData_DownGenes)
Fit_CPA3 <- survfit(Surv(Time, Event) ~ CPA3, data = SurvData_DownGenes)
Fit_BRE <- survfit(Surv(Time, Event) ~ BRE, data = SurvData_DownGenes)
Fit_UFM1 <- survfit(Surv(Time, Event) ~ UFM1, data = SurvData_DownGenes)
Fit_DPT <- survfit(Surv(Time, Event) ~ DPT, data = SurvData_DownGenes)
Fit_CHRNA2 <- survfit(Surv(Time, Event) ~ CHRNA2, data = SurvData_DownGenes)
Fit_NT5E <- survfit(Surv(Time, Event) ~ NT5E, data = SurvData_DownGenes)
Fit_LTBR <- survfit(Surv(Time, Event) ~ LTBR, data = SurvData_DownGenes)
Fit_EDN3 <- survfit(Surv(Time, Event) ~ EDN3, data = SurvData_DownGenes)
Fit_COL4A5 <- survfit(Surv(Time, Event) ~ COL4A5, data = SurvData_DownGenes)

##################################################################

## Plot Up genes

Up_Plot1 <- ggsurvplot(Fit_TMSB10,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "TMSB10")

Up_Plot2 <- ggsurvplot(Fit_ASPN,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "ASPN")

Up_Plot3  <-  ggsurvplot(Fit_P4HA2,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "P4HA2")

Up_Plot4 <- ggsurvplot(Fit_PTDSS1,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "PTDSS1")


Up_Plot5 <-  ggsurvplot(Fit_NCOA2,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "NCOA2")


Up_PlotList <- list(Up_Plot1, Up_Plot2, Up_Plot3, Up_Plot4, Up_Plot5)
names(Up_PlotList) <- PositiveGenes

SplotUp <- arrange_ggsurvplots(Up_PlotList, title = "Survival plots of up-regulated genes",ncol = 2, nrow = 3)
ggsave("UpGenes3.pdf", SplotUp, width = 40, height = 30, units = "cm")

###############################################################
## Plot Bad genes

Down_Plot1 <- ggsurvplot(Fit_AZGP1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "AZGP1")

Down_Plot2 <- ggsurvplot(Fit_CPA3,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "CPA3")

Down_Plot3 <- ggsurvplot(Fit_BRE,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "BRE")


Down_Plot4 <- ggsurvplot(Fit_UFM1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "UFM1")

Down_Plot5 <- ggsurvplot(Fit_DPT,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "DPT")



Down_Plot6 <- ggsurvplot(Fit_CHRNA2,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "CHRNA2")


Down_Plot7 <- ggsurvplot(Fit_NT5E,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "NT5E")

Down_Plot8 <- ggsurvplot(Fit_LTBR,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "LTBR")

Down_Plot9 <- ggsurvplot(Fit_EDN3,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "EDN3")


Down_Plot10 <- ggsurvplot(Fit_COL4A5,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "COL4A5")




Down_PlotList <- list(Down_Plot1, Down_Plot2, Down_Plot3, Down_Plot4, Down_Plot5, Down_Plot6, Down_Plot7, Down_Plot8, Down_Plot9, Down_Plot10)
names(Down_PlotList) <- NegativeGenes

SplotDown <- arrange_ggsurvplots(Down_PlotList, title = "Survival plots of down-regulated genes", ncol = 2, nrow = 5)
ggsave("DownGenes3.pdf", SplotDown, width = 40, height = 30, units = "cm")
