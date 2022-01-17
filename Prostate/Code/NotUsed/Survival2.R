## Clean working space
rm(list = ls())

## Set the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

## Load Dataset1
load("./Data/Dataset6.rda")

## Load val_pheno1 
load("./Objs/val_pheno1_Zscore.rda")

## Get expression matrix
expr6 <- Dataset6$expr

## Annotate expr1
#rownames(expr6) <- Dataset6$keys
#summary(is.na(rownames(expr6)))
#expr1 <- expr1[!(rownames(expr1) == ""), ]
#expr1 <- expr1[!is.na(rownames(expr1)), ]
#dim(expr1)

## Get the phenotype
pheno6 <- Dataset6$pheno

table(pheno6$`extreme:ch1`)
pheno6$Metastasis <- pheno6$`extreme:ch1`

pheno6$Metastasis[pheno6$Metastasis == "Indolent"] <- "No_Mets"
pheno6$Metastasis[pheno6$Metastasis == "Lethal"] <- "Mets"

pheno6$Metastasis <- as.factor(pheno6$Metastasis)
table(pheno6$Metastasis)
all(rownames(pheno6) == colnames(expr6))

Dataset6$expr <- expr6
Dataset6$pheno <- pheno6

## Keep only the relevant information (Metastasis Event and Time)

pheno6$Event <- as.numeric(pheno6$Metastasis)
pheno6$Event[pheno6$Event == 2] <- 0 

pheno6$Time <- as.numeric(pheno6$`fup.month:ch1`)

Phenotype <- pheno6[,-c(1:44)]


# Keep only the Good AND Bad Genes (Associated with no metastasis // Metastasis)
exprUP <- expr6[rownames(expr6)%in%PositiveGenes, ]
exprDown <- expr6[rownames(expr6)%in%NegativeGenes, ]

# Z-transform exprGood
exprUP <- t(scale(t(exprUP)))
exprDown <- t(scale(t(exprDown)))

# create a merged pdata and Z-scores object
CoxDataUp <- data.frame(Phenotype, t(exprUP))
CoxDataDown <- data.frame(Phenotype, t(exprDown))


#################################################################################
## Define optimal cutpoints for each gene (converting the absolute expression into categorical low/high expression)
PositiveGenes <- rownames(exprUP)
NegativeGenes <- rownames(exprDown)

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
ggsave("UpGenes2.pdf", SplotUp, width = 40, height = 30, units = "cm")

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
ggsave("DownGenes2.pdf", SplotDown, width = 40, height = 30, units = "cm")

