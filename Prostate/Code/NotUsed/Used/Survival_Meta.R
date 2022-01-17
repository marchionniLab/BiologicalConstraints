## Mohamed Omar
## 21/07/2019
## Goal: Survival analysis in the Testing data sets using the Meta score from the gene signature predicting metastasis
################################################################################################

## Clean working space
rm(list = ls())

## Set the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

## Load necessary packages
library(survival)
library(survminer)

############################################################################

## Load val_pheno1 
load("./Objs/val_pheno1_Zscore.rda")

val_pheno3[,c(1:5)] <- NULL

#val_pheno6$`tumour gleason:ch1` <- gsub("\\=.+", "", val_pheno6$`tumour gleason:ch1`)

## Turn the score column into numeric
val_pheno1[,"score"] <- as.numeric(val_pheno1[,"score"]) 
val_pheno2[,"score"] <- as.numeric(val_pheno2[, "score"])
val_pheno3[,"score"] <- as.numeric(val_pheno3[, "score"])
#val_pheno6[,"score"] <- as.numeric(val_pheno6[, "score"])


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

# val_pheno3
val_pheno3$Meta_Score <- val_pheno3$score
val_pheno3[which(val_pheno3$score <= thr_test3), "Meta_Score"] = "low"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
val_pheno3[which(val_pheno3$score > thr_test3),"Meta_Score"] = "high"

# val_pheno6
#val_pheno6$Meta_Score <- val_pheno6$score
#val_pheno6[which(val_pheno6$score <= thr_test6), "Meta_Score"] = "low"
#val_pheno1[which(val_pheno1$score > quantiles[1] &  val_pheno1$score < quantiles[2]),"Meta_Score"] = "int"
#val_pheno6[which(val_pheno6$score > thr_test6),"Meta_Score"] = "high"


val_pheno1$Time <- as.numeric(val_pheno1$`follow-up time (met, months):ch1`)
val_pheno2$Time <- as.numeric(val_pheno2$`fup.month:ch1`)
val_pheno3$Time <- as.numeric(val_pheno3$met_time)
#val_pheno6$Time <- as.numeric(val_pheno6$`total follow up (months):ch1`)


val_pheno1$Event <- as.numeric(val_pheno1$`met event (1=yes, 0=no):ch1`)

val_pheno2$Event <- as.numeric(val_pheno2$Metastasis)
val_pheno2$Event[val_pheno2$Event == 2] <- 0 

val_pheno3$Event <- as.numeric(val_pheno3$Metastasis)
val_pheno3$Event[val_pheno3$Event == 2] <- 0 
table(val_pheno3$Event)

#val_pheno6$Event <- as.numeric(val_pheno6$Metastasis)-1
#table(val_pheno6$Event)

## Keep only relevant information from the phenotype table
surv_data1 <- val_pheno1[,c("Time","Event","score","Meta_Score","gleason grade:ch1", "patient age (years):ch1", "psa (ng/ml):ch1", "t-stage:ch1")]
surv_data1$Meta_Score <- factor(surv_data1$Meta_Score, levels = c("low", "high"))

surv_data2 <- val_pheno2[, c("Time", "Event", "score", "Meta_Score", "gleason:ch1", "age:ch1", "cancer.percent:ch1")]
surv_data2$Meta_Score <- factor(surv_data2$Meta_Score, levels = c("low", "high"))

surv_data3 <- val_pheno3[, c("Time", "Event", "score", "Meta_Score", "Re.Gleasonsum", "age.y", "ERG", "re.PositiveMargin", "re.Stage", "PTEN_LOSS")]
surv_data3$Meta_Score <- factor(surv_data3$Meta_Score, levels = c("low", "high"))

#surv_data6 <- val_pheno6[, c("Time", "Event", "score", "Meta_Score", "tumour gleason:ch1", "psa at diag:ch1", "tumour %:ch1", "positive surgical margins (psm):ch1", "clinical stage:ch1", "biochemical relapse (bcr):ch1")]
#surv_data6$Meta_Score <- factor(surv_data6$Meta_Score, levels = c("low", "high"))

# Create a survival object
surv_data.surv1 <-  with(surv_data1, Surv(Time, Event == 1))
surv_data.surv2 <-  with(surv_data2, Surv(Time, Event == 1))
surv_data.surv3 <-  with(surv_data3, Surv(Time, Event == 1))
#surv_data.surv6 <-  with(surv_data6, Surv(Time, Event == 1))


#Calculate p-value
survdifftest1 <- survdiff(surv_data.surv1 ~ Meta_Score, data = surv_data1)
survpvalue1 <- 1 - pchisq(survdifftest1$chisq, length(survdifftest1$n) - 1)

survdifftest2 <- survdiff(surv_data.surv2 ~ Meta_Score, data = surv_data2)
survpvalue2 <- 1 - pchisq(survdifftest2$chisq, length(survdifftest2$n) - 1)

survdifftest3 <- survdiff(surv_data.surv3 ~ Meta_Score, data = surv_data3)
survpvalue3 <- 1 - pchisq(survdifftest3$chisq, length(survdifftest3$n) - 1)

# survdifftest6 <- survdiff(surv_data.surv6 ~ Meta_Score, data = surv_data6)
# survpvalue6 <- 1 - pchisq(survdifftest6$chisq, length(survdifftest6$n) - 1)

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
png(filename = "./Figs/Meta/KM_Survival_GSE116918.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data1),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the gene signature Meta score in GSE116918", xlab = "Time in months", ylab = "Survival")

dev.off()

###############

png(filename = "./Figs/Meta/KM_Survival_GSE16560.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data2),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the gene signature Meta score in the Test data", xlab = "Time in months", ylab = "Survival")

dev.off()

###################
png(filename = "./Figs/Meta/KM_Survival_JHU.png", width = 2000, height = 2000, res = 300)
ggsurvplot(survfit(Surv(Time, Event) ~ Meta_Score, data = surv_data3),
           risk.table = "nrisk_cumevents",
           pval = TRUE,
           ggtheme = theme_bw(),
           surv.median.line = "hv",
           pval.method = TRUE, 
           risk.table.y.text.col = TRUE,
           risk.table.y.text = TRUE, palette = c("green", "red"), title = "Survival by the gene signature Meta score in JHU cohort", xlab = "Time in months", ylab = "Survival")

dev.off()

#################

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
UniCOX1 <- as.data.frame(UniCOX1)
UniCOX1
write.csv(UniCOX1, file = "./Objs/Meta/UniCox_GSE116918.csv")
########################################

## Multivariate

CoxModel1 <- coxph(Surv(Time, Event)~ Meta_Score+Age+PSA+Gleason+T_Stage, data=Cox_data1)

png(filename = "./Figs/Meta/CoxModel_GSE116918.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel1, data = Cox_data1, main = "Cox proportional hazards model in GSE116918")
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

#######################################################
#######################################################
## Third validation dataset
Cox_data3 <- surv_data3
#Cox_data$score <- as.numeric(Cox_data$score)
#Cox_CutPoint <- surv_cutpoint(data = Cox_data, time = "Time", event = "Event", variables = "score")
#Cox_CutPoint
#Cox_data_Categ <- surv_categorize(Cox_CutPoint)
#CoxData_Final <- cbind(Cox_data_Categ, Cox_data[4:7])
#CoxData_Final$score <- factor(CoxData_Final$score, levels = c("low", "high"))

Cox_data3$Age <- ifelse(Cox_data3$age.y >= 65, ">=65 years", "<65 years")
Cox_data3$Age <- factor(Cox_data3$Age, levels = c("<65 years", ">=65 years"))

#Cox_data1$PSA <- ifelse(Cox_data1$`psa (ng/ml):ch1` > 20, "High PSA (> 20)", "Low PSA (<= 20)")
#Cox_data1$PSA <- factor(Cox_data1$PSA, levels = c("Low PSA (<= 20)", "High PSA (> 20)"))

Cox_data3$Gleason <- ifelse(Cox_data3$Re.Gleasonsum <= 7 , "Gleason <= 7", "Gleason > 7")
Cox_data3$Gleason <- factor(Cox_data3$Gleason, levels = c("Gleason <= 7", "Gleason > 7"))

Cox_data3$ERG <- as.factor(Cox_data3$ERG)
table(Cox_data3$ERG)
levels(Cox_data3$ERG)

Cox_data3$Margin <- as.factor(Cox_data3$re.PositiveMargin)
table(Cox_data3$Margin)
levels(Cox_data3$Margin) <- c("Negative", "Negative", "Positive", "Positive")


Cox_data3$Stage <- as.factor(Cox_data3$re.Stage)
table(Cox_data3$Stage)
levels(Cox_data3$Stage) <- c("< T3", "< T3", "< T3", "< T3", "< T3", "< T3", "T3", "T3")

Cox_data3$PTEN <- as.factor(Cox_data3$PTEN_LOSS)
table(Cox_data3$PTEN)
levels(Cox_data3$PTEN)

## Univariate Cox model
# PFS
covariates3 <- c("Meta_Score", "Age", "Gleason", "ERG", "Margin", "Stage", "PTEN")
univ_formulas3 <- sapply(covariates3,
                         function(x) as.formula(paste('Surv(Time, Event)~', x)))

univ_models3 <- lapply( univ_formulas3, function(x){coxph(x, data = Cox_data3)})
# Extract data 
univ_results3 <- lapply(univ_models3,
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

UniCOX3 <- t(as.data.frame(univ_results3, check.names = F))
UniCOX3 <- as.data.frame(UniCOX3)
UniCOX3
write.csv(UniCOX3, file = "./Objs/Meta/UniCox_JHU.csv")

########################################

## Multivariate

CoxModel3 <- coxph(Surv(Time, Event)~ Meta_Score+Age+Gleason+ERG+Margin+Stage+PTEN, data=Cox_data3)

png(filename = "./Figs/Meta/CoxModel3_JHU.png", width = 3000, height = 2000, res =400)
ggforest(CoxModel3, data = Cox_data3, main = "Cox proportional hazards model in JHU cohort")
dev.off()

####################################################################
