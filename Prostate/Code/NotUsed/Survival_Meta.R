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
