## Clean working space
rm(list = ls())

## Set the working directory
setwd("/Users/mohamedomar/Documents/Research/Projects/Prostate")

## Load Dataset1
load("./Data/Dataset1.rda")

## Load val_pheno1 
load("./Objs/val_pheno1_Zscore.rda")

## Get expression matrix
expr1 <- Dataset1$expr

## Annotate expr1
head(rownames(expr1))
rownames(expr1) <- Dataset1$keys
expr1 <- expr1[!is.na(rownames(expr1)), ]

rownames(expr1) <- gsub("\\,.+", "", rownames(expr1))
rownames(expr1) <- gsub("\\-.+", "", rownames(expr1))

expr1 <- aggregate(expr1[,], list(Gene = rownames(expr1)), FUN = mean)
rownames(expr1) <- expr1$Gene
expr1$Gene <- NULL

expr1 <- as.matrix(expr1)
dim(expr1)
#

## Get the phenotype
pheno1 <- Dataset1$pheno
pheno1$Metastasis <- pheno1$`met event (1=yes, 0=no):ch1`
pheno1$Metastasis[pheno1$Metastasis == 0] <- "No_Mets"
pheno1$Metastasis[pheno1$Metastasis == 1] <- "Mets"
pheno1$Metastasis <- as.factor(pheno1$Metastasis)
table(pheno1$Metastasis)
all(rownames(pheno1) == colnames(expr1))

Dataset1$pheno <- pheno1
Dataset1$expr <- expr1
Dataset1$keys <- rownames(expr1)


## Keep only the relevant information (Metastasis Event and Time)

pheno1$Event <- as.numeric(pheno1$`met event (1=yes, 0=no):ch1`)
pheno1$Time <- as.numeric(pheno1$`follow-up time (met, months):ch1`)

Phenotype <- pheno1[,-c(1:44)]


# Keep only the Good AND Bad Genes (Associated with no metastasis // Metastasis)
exprUP <- expr1[rownames(expr1) %in% PositiveGenes, ]
exprDown <- expr1[rownames(expr1) %in% NegativeGenes, ]

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

Fit_ASNS<- survfit(Surv(Time, Event) ~ ASNS , data = SurvData_UpGenes)
Fit_CADPS <- survfit(Surv(Time, Event) ~ CADPS, data = SurvData_UpGenes)
Fit_CAMK2N1 <- survfit(Surv(Time, Event) ~ CAMK2N1, data = SurvData_UpGenes)
Fit_ENSA <- survfit(Surv(Time, Event) ~ ENSA, data = SurvData_UpGenes)
Fit_FOXH1 <- survfit(Surv(Time, Event) ~ FOXH1, data = SurvData_UpGenes)

Fit_GABRD <- survfit(Surv(Time, Event) ~ GABRD, data = SurvData_UpGenes)
Fit_GNPTAB <- survfit(Surv(Time, Event) ~ GNPTAB, data = SurvData_UpGenes)
Fit_HLTF <- survfit(Surv(Time, Event) ~ HLTF, data = SurvData_UpGenes)
Fit_KBTBD2 <- survfit(Surv(Time, Event) ~ KBTBD2, data = SurvData_UpGenes)
Fit_MRPL11 <- survfit(Surv(Time, Event) ~ MRPL11, data = SurvData_UpGenes)
Fit_NCOA2 <- survfit(Surv(Time, Event) ~ NCOA2, data = SurvData_UpGenes)
Fit_OLR1 <- survfit(Surv(Time, Event) ~ OLR1, data = SurvData_UpGenes)
Fit_PDAP1 <- survfit(Surv(Time, Event) ~ PDAP1, data = SurvData_UpGenes)
Fit_PTPN9 <- survfit(Surv(Time, Event) ~ PTPN9, data = SurvData_UpGenes)
Fit_RFTN1 <- survfit(Surv(Time, Event) ~ RFTN1, data = SurvData_UpGenes)
Fit_RNF19A <- survfit(Surv(Time, Event) ~ RNF19A, data = SurvData_UpGenes)
Fit_SOX4 <- survfit(Surv(Time, Event) ~ SOX4, data = SurvData_UpGenes)
Fit_TOP2A <- survfit(Surv(Time, Event) ~ TOP2A, data = SurvData_UpGenes)
Fit_TRIAP1 <- survfit(Surv(Time, Event) ~ TRIAP1, data = SurvData_UpGenes)
Fit_ZNF512B <- survfit(Surv(Time, Event) ~ ZNF512B, data = SurvData_UpGenes)



##################################################################
## Fit down genes

Fit_AZGP1 <- survfit(Surv(Time, Event) ~ AZGP1, data = SurvData_DownGenes)
Fit_CBLL1 <- survfit(Surv(Time, Event) ~ CBLL1, data = SurvData_DownGenes)
Fit_CHRNA2 <- survfit(Surv(Time, Event) ~ CHRNA2, data = SurvData_DownGenes)
Fit_CPA3 <- survfit(Surv(Time, Event) ~ CPA3, data = SurvData_DownGenes)
Fit_DCXR <- survfit(Surv(Time, Event) ~ DCXR, data = SurvData_DownGenes)
Fit_DPT <- survfit(Surv(Time, Event) ~ DPT, data = SurvData_DownGenes)
Fit_EDN3 <- survfit(Surv(Time, Event) ~ EDN3, data = SurvData_DownGenes)
Fit_KCTD14 <- survfit(Surv(Time, Event) ~ KCTD14, data = SurvData_DownGenes)
Fit_KLF4 <- survfit(Surv(Time, Event) ~ KLF4, data = SurvData_DownGenes)
Fit_KRT9 <- survfit(Surv(Time, Event) ~ KRT9, data = SurvData_DownGenes)

Fit_LUZP2 <- survfit(Surv(Time, Event) ~ LUZP2, data = SurvData_DownGenes)
Fit_NT5E <- survfit(Surv(Time, Event) ~ NT5E, data = SurvData_DownGenes)
Fit_PART1 <- survfit(Surv(Time, Event) ~ PART1, data = SurvData_DownGenes)
Fit_PTN <- survfit(Surv(Time, Event) ~ PTN, data = SurvData_DownGenes)
Fit_PTPRN2 <- survfit(Surv(Time, Event) ~ PTPRN2, data = SurvData_DownGenes)
Fit_SIDT1 <- survfit(Surv(Time, Event) ~ SIDT1, data = SurvData_DownGenes)
Fit_UFM1 <- survfit(Surv(Time, Event) ~ UFM1, data = SurvData_DownGenes)
Fit_WNT8B <- survfit(Surv(Time, Event) ~ WNT8B, data = SurvData_DownGenes)

##################################################################

## Plot Up genes

Up_Plot1 <- ggsurvplot(Fit_ASNS,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "ASNS")

Up_Plot2 <- ggsurvplot(Fit_CADPS,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "CADPS")

Up_Plot3  <-  ggsurvplot(Fit_CAMK2N1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "CAMK2N1")

Up_Plot4 <- ggsurvplot(Fit_ENSA,
                       risk.table = FALSE,
                       pval = TRUE,
                       ggtheme = theme_minimal(),
                       risk.table.y.text.col = FALSE,
                       risk.table.y.text = FALSE, title = "ENSA")


Up_Plot5 <-  ggsurvplot(Fit_FOXH1,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "FOXH1")

Up_Plot6 <-  ggsurvplot(Fit_GABRD,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "GABRD")


Up_Plot7 <-  ggsurvplot(Fit_GNPTAB,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "GNPTAB")


Up_Plot8 <-  ggsurvplot(Fit_HLTF,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "HLTF")


Up_Plot9 <-  ggsurvplot(Fit_KBTBD2,
                        risk.table = FALSE,
                        pval = TRUE,
                        ggtheme = theme_minimal(),
                        risk.table.y.text.col = FALSE,
                        risk.table.y.text = FALSE, title = "KBTBD2")


Up_Plot10 <-  ggsurvplot(Fit_MRPL11,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "MRPL11")


Up_Plot11 <-  ggsurvplot(Fit_NCOA2,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "NCOA2")


Up_Plot12 <-  ggsurvplot(Fit_OLR1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "OLR1")


Up_Plot13 <-  ggsurvplot(Fit_PDAP1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "PDAP1")


Up_Plot14 <-  ggsurvplot(Fit_PTPN9,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "PTPN9")


Up_Plot15 <-  ggsurvplot(Fit_RFTN1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "RFTN1")

Up_Plot16 <-  ggsurvplot(Fit_RNF19A,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "RNF19A")

Up_Plot17 <-  ggsurvplot(Fit_SOX4,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "SOX4")

Up_Plot18 <-  ggsurvplot(Fit_TOP2A,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "TOP2A")

Up_Plot19 <-  ggsurvplot(Fit_TRIAP1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "TRIAP1")

Up_Plot20 <-  ggsurvplot(Fit_ZNF512B,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "ZNF512B")



Up_PlotList <- list(Up_Plot1, Up_Plot2, Up_Plot3, Up_Plot4, Up_Plot5, Up_Plot6, Up_Plot7, Up_Plot8, Up_Plot9, Up_Plot10, Up_Plot11, Up_Plot12, Up_Plot13, Up_Plot14, Up_Plot15, Up_Plot16, Up_Plot17, Up_Plot18, Up_Plot19, Up_Plot20)
names(Up_PlotList) <- PositiveGenes

SplotUp <- arrange_ggsurvplots(Up_PlotList, title = "Survival plots of up-regulated genes (GSE116918)",ncol = 2, nrow = 5)
ggsave("./Figs/Meta/UpGenes_GSE116918.pdf", SplotUp, width = 40, height = 30, units = "cm")

###############################################################
## Plot Bad genes

Down_Plot1 <- ggsurvplot(Fit_AZGP1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "AZGP1")

Down_Plot2 <- ggsurvplot(Fit_CBLL1,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "CBLL1")

Down_Plot3 <- ggsurvplot(Fit_CHRNA2,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "CHRNA2")


Down_Plot4 <- ggsurvplot(Fit_CPA3,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "CPA3")

Down_Plot5 <- ggsurvplot(Fit_DCXR,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "DCXR")



Down_Plot6 <- ggsurvplot(Fit_DPT,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "DPT")


Down_Plot7 <- ggsurvplot(Fit_EDN3,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "EDN3")

Down_Plot8 <- ggsurvplot(Fit_KCTD14,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "KCTD14")

Down_Plot9 <- ggsurvplot(Fit_KLF4,
                         risk.table = FALSE,
                         pval = TRUE,
                         ggtheme = theme_minimal(),
                         risk.table.y.text.col = FALSE,
                         risk.table.y.text = FALSE, title = "KLF4")


Down_Plot10 <- ggsurvplot(Fit_KRT9,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "KRT9")

Down_Plot11 <- ggsurvplot(Fit_LUZP2,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "LUZP2")

Down_Plot12 <- ggsurvplot(Fit_NT5E,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "NT5E")

Down_Plot13 <- ggsurvplot(Fit_PART1,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "PART1")

Down_Plot14 <- ggsurvplot(Fit_PTN,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "PTN")

Down_Plot15 <- ggsurvplot(Fit_PTPRN2,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "PTPRN2")

Down_Plot16 <- ggsurvplot(Fit_SIDT1,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "SIDT1")

Down_Plot17 <- ggsurvplot(Fit_UFM1,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "UFM1")

Down_Plot18 <- ggsurvplot(Fit_WNT8B,
                          risk.table = FALSE,
                          pval = TRUE,
                          ggtheme = theme_minimal(),
                          risk.table.y.text.col = FALSE,
                          risk.table.y.text = FALSE, title = "WNT8B")




Down_PlotList <- list(Down_Plot1, Down_Plot2, Down_Plot3, Down_Plot4, Down_Plot5, Down_Plot6, Down_Plot7, Down_Plot8, Down_Plot9, Down_Plot10, Down_Plot11, Down_Plot12, Down_Plot13, Down_Plot14, Down_Plot15, Down_Plot16, Down_Plot17, Down_Plot18)
names(Down_PlotList) <- NegativeGenes

SplotDown <- arrange_ggsurvplots(Down_PlotList, title = "Survival plots of down-regulated genes (GSE116918)", ncol = 3, nrow = 6)
ggsave("./Figs/Meta/DownGenes_GSE116918.pdf", SplotDown, width = 40, height = 30, units = "cm")

