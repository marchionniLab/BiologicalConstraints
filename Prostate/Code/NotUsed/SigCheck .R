
library(SigCheck)

## Load Testing data set
load("./Data/Dataset1.rda")

valDataSet <- Dataset1

# Acess the phenotype and expression data
val_pheno1 <- Dataset1$pheno
val_expr1 <- Dataset1$expr

## Modify val_pheno1
val_pheno1$Metastasis <- val_pheno1$`met event (1=yes, 0=no):ch1`
val_pheno1$Metastasis[val_pheno1$Metastasis == 0] <- "No_Mets"
val_pheno1$Metastasis[val_pheno1$Metastasis == 1] <- "Mets"
val_pheno1$Metastasis <- as.factor(val_pheno1$Metastasis)
table(val_pheno1$Metastasis)
all(rownames(val_pheno1) == colnames(val_expr1))


FData <- data.frame(Genes = rownames(val_expr1), GeneSymbol = valDataSet$keys)

ValData <- ExpressionSet(assayData = val_expr1, phenoData = AnnotatedDataFrame(val_pheno1), 
                         featureData = AnnotatedDataFrame(FData), annotation = "GPL25318")

ValData$`follow-up time (met, months):ch1` <- as.numeric(ValData$`follow-up time (met, months):ch1`)


Check <- sigCheck(ValData, classes="Metastasis",
                   signature=knownSignatures$cancer$ADORNO,
                   annotation="GeneSymbol", scoreMethod = "classifier")



Pos <- c("TMSB10", "IQGAP3", "CST1", "RPRML", "CST2", "RHOT1", "STC2", "FOXH1", "PTDSS1", "HES6", "TRIM59")
Neg <- c("AZGP1", "NT5DC1", "SLC13A2", "CNGB1", "CPA3", "WNT6", "KCTD14", "PTPRN2", "BRE", "TTLL9", "ADRBK2", "KLHDC7A", "MIB2", "UFM1", "CCK", "KIAA1210", "POTEG")
NewClassifier <- c(Pos, Neg)

Check <- sigCheck(ValData, classes="Metastasis", survival = "follow-up time (met, months):ch1",
                  signature=NewClassifier,
                  annotation="GeneSymbol", scoreMethod = "High", threshold = 0.90)



Check

Known <- c(OncoDx_Filter$posGeneNames, OncoDx_Filter$negGeneNames)
Known <- list(Known = Known)

Check_Compare <- sigCheckKnown(Check, known = Known$Known)
Check_Compare

sigCheckPlot(Check_Compare)






