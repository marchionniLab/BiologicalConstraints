

EmtGenes <- c("CDH1", "VIM", "CDH2", "FN1", "DSP", "PKP1", "MAP7", "ESRP1", 
              "KRT6A", "C1orf116", "SERPINB5", "LASS3", "LAD1", "PERP",
              "JUP", "DSG3", "DSC3", "OVOL1", "GRHL2", "ERBB3", "MARVELD3",
              "EVPL", "TRIM29", "FAM83B", "ESRP2", "CDS1", "EPN3", "EPHA1", 
              "KLF5", "ABLIM1", "PPL", "GRHL3", "HAS3", "ZNF750", "A2ML1", 
              "PRRG4", "MAL2", "RAB25", "TTC22", "IRF6", "GRHL1", "VIM", 
              "ANXA6", "LIX1L", "HLX", "SYT11", "GNAI2", "FYN", "ATP10A", 
              "CNRIP1", "EMILIN1", "DAB2", "OLFML2B", "COL6A1", "ZEB2", 
              "COL6A2", "ST3GAL2", "OLFML3", "ACVRL1",  "CMTM3", "PMP22", 
              "PCOLCE", "GPR124", "TIMP2", "LAMA4", "CDH2", "DACT1", "FN1", 
              "FSTL1", "FAP", "CALD1", "COL5A1", "VCAN", "COL5A2", "GLT8D2", 
              "POSTN", "COL1A1", "ADAM12", "COL6A3", "SPARC", "ZNF469", "COL1A2", 
              "NID2", "COL3A1", "PDGFRB")

EmtGenes <- EmtGenes[!duplicated(EmtGenes)]

MetSuppressors <- c("NME1", "CD82", "KISS1", "BRMS1", "MAP2K4", "CD44", "SERPINB5", "SDPR", "MED23", "KLF17")

myTSPs2 <- expand.grid(EmtGenes, MetSuppressors)
myTSPs2 <- as.matrix(myTSPs2)


OncoGenesEMT <- c("SKI", "MIR135B",	"RAB22A",	"FASN",	"FHL2",	"ELK1",	"PPM1D", "AXL",	"HMGA2",	"CTSZ",	"PIK3R1",
                  "GLI2",	"MUC1",	"HOXA9", "MIR130B", "MYCN", "TBX2",	"PTP4A2",	"ID1",	"BCL6",	"MIR106B", "GATA1",
                  "AKT1",	"PAK5",	"JUNB",	"TBX3",	"CUL4A",	"TXN", "CAMK1D",	"PRKCA",	"USP22",	"BMI1",	"KLF8",
                  "EIF3I",	"MYC",	"RAF1",	"MIR221",	"NEAT1",	"BRD4",	"TRIM28",	"BRAF",	"EPS8",	"FGFR2", "FRAT1",
                  "MTDH",	"PTP4A3",	"EPCAM",	"PDGFRB",	"MIR191",	"PRKCI",	"MIR520G",	"YWHAZ",	"GLO1",	"TLE1",	"BMP7",
                  "ARHGEF2",	"KDM5B",	"HRAS",	"MIR517C", "MIR224",	"SQSTM1",	"ROS1",	"NCOA3", "AR",	"AURKA", "CCND1",
                  "GSK3A",	"MIR92B",	"WWTR1", "ETS1", "SNAI1",	"EIF5A2",	"CIP2A", "PDGFB",	"BIRC2", "ETV1",	"GMNN",
                  "PAK1",	"ERBB2", "CDK14",	"S100A4",	"MIR93", "ZNF217", "UHRF1",	"EEF1D", "NANOGP8",	"RAC1",	"HSPA4",
                  "MIR301A", "HOXD9",	"SOX4",	"GATA6", "MKL2",	"ACTN4", "ROCK1",	"FOSL1",	"SATB1",	"BCL2",	"EGFR",
                  "ALK",	"MIR373",	"SETDB1",	"MIR10B",	"MET",	"MIR663A", "TFCP2",	"HOTAIR",	"LETMD1",	"MUC4",	"SOX2",
                  "HSPB1", "RHOC",	"MYB", "FOXQ1",	"ETV4",	"ERG",	"CTNNB1",	"MITF",	"EIF4E",	"HSPA5", "MIR1281",
                  "CXCR4", "MIR96",	"MIR19A",	"PIK3CA",	"BRF2",	"TYMS",	"SREBF1",	"BCL2L1",	"KIT", "GAB2", "EML4",
                  "YWHAG", "LIN28A", "PIM2", "FZD2", "MIR21", "CD24", "LEF1", "SIX1",	"WNT1", "ELK3", "MTOR",
                  "MCL1", "CRK",	"JUN",	"TAZ", "FGFR1",	"IRS2", "GNA13", "ID2",	"NOTCH4", "CBLB",	"MLLT3",
                  "IGF1R",	"MIR506",	"KRAS",	"FOXM1", "ZEB1AS1",	"AKTIP",	"MDM2",	"ELAVL1",	"MACC1", "BIRC5",	"CRYAB",
                  "MDM4", "S100A8",	"CEACAM6", "AKT2", "MALAT1",	"ITGA3", "YY1",	"UCA1",	"CRKL", "MIR454", "SRC",
                  "MAP3K7",	"CREB1", "MYD88",	"LCN2",	"MIR485",	"NANOG",	"PAX2",	"YBX1", "UBE3C",	"ALDH1A1",	"PLAGL2",
                  "TWIST1", "GOLPH3",	"JAK2", "PLAC8")

TumorSuppEMT <- c("CREBBP",	"LIMA1", "SMAD2", "FBP1", "GSK3B", "PKD1", "WNT11", "TWIST2",	"ESRP1", "MIR888", "IFT88",
                  "GPC3",	"MIR34C", "ANXA1",	"BMP4", "EPAS1",	"NDRG2", "MIR136", "ING4",	"FLNA", "VEGFA", "TUSC7",
                  "NUMB",	"RUNX2",	"DAB2",	"MIR2182", "MIR1011",	"ESR1",	"TDGF1", "PARP1", "FOXA2", "UHRF2", "MIR186",
                  "MIR1242", "MIR30C1", "ARHGEF12",	"DKK1", "MIR29A",	"MIR146A", "IGF1", "SPOP", "LRIG1", "ANGPTL4", "CDKN2A",
                  "TP53BP2", "HPGD", "PHLDA2", "AXIN2", "GKN2",	"ATM", "MIR497", "DNMT3B",	"KRT19",	"MIR449A", "MEG3",
                  "MIR25", "ABCG2", "MIR143",	"ESRRB", "MIR193A", "FBXW7", "PTPN6", "MIR200B", "SOCS3",	"PAWR",	"EPHB2",
                  "AFAP1L2", "MIR196B", "ZYX", "TGFBR3", "HRG",	"GJB2",	"LATS1", "BRMS1",	"AXIN1", "MIRLET7B", "CDH13",
                  "LOX", "HIPK2", "MIR137", "MIR241", "ROR2",	"FOXC1", "FHL1", "MIR2181",	"PRKAA2", "BMP2",	"MIR185",
                  "HDAC3", "RNH1", "ESR2",	"EGLN3", "CD44", "FOXO4", "CSF2", "MIR1012", "ZFP36", "RNF111",	"MIR26B",
                  "MIRLET7A1", "NEDD4L","NOTCH2",	"TGFBR2",	"MIRLET7G", "KDM3A", "FAS",	"RASAL2", "MARVELD1", "DAB2IP", "ERF",
                  "VDR", "H2AFX", "SCUBE2", "MIR375", "MIR206", "WWOX", "DICER1",	"NF1", "TNFSF12", "NFKB1", "AJAP1",
                  "NR1I2", "MIR29B1", "ISG15","TSC1", "MIR132", "MIR211", "BRD7", "CXCL12", "TSC2", "CEBPA", "MIR15A",
                  "CEACAM1", "SIRT3", "BBC3", "HNF4A", "PDCD4",	"CD82", "DDR2", "CCNDBP1", "THBD",	"MIR33A", "EED",
                  "CTNNBIP1",	"DKK3", "MIR203A", "STAT5A", "IL17A", "ITGB3", "FBLN1",	"SCRIB", "TP53", "PIN1", "NOTCH3",
                  "IGFBP3", "MIR152", "SMAD4", "SLC9A3R1", "ITGA5", "WISP3", "VSNL1",	"PROX1", "MIR1241", "MAP3K4",	"PTPA",
                  "MIR122", "MIR141", "PRKAA1", "S100A2", "MTUS1", "DAPK1", "MIR106A", 	"HIF1A", "CDH5", "MIR26A2", "NDRG1",
                  "IGFBP7", "SIAH1",	"TXNIP", "MIR91", "CXCR2",	"KL",	"UIMC1", "ERRFI1", "RNF8", "TIMP3", "IRF8",
                  "VHL", "GSN", "STK11",	"CLU", "NR2C2",	"MIR1291", "MIR10A", "MIR29C", "SOD2", "SFRP2",	"CXCL14",
                  "YPEL3", "KDM5A", "AHR", "BATF2", "NUAK1", "RHOB", "EPB41L3",	"MIR424",	"MIR181A2", "GLS2",	"MIR148A",
                  "MIR494", "INPP4B",	"CLDN1", "MIR26A1",	"RACK1", "MIR181A1",	"MIR1243", 	"KLK6", "BRCA1", "TP53BP1",	"MIR187",
                  "CTNND1", "RB1",	"TCF4", "MIR205", "MIR214", "DNMT1", "MIR204", "ITGB1",	"EGR1", "TP53INP1", "CFTR",
                  "EPHB3", 	"LEFTY1",	"MIR1941",	"PEBP1", "STAT1", "PCDH9", "MIR217", "FOXO3", "PTEN",	"TRIM62", "HIC1",
                  "MIR145",	"PCGF2", "MIR181B1", "KDM6A",	"TNFAIP8L2", "SUFU", "MIR200A",	"MIR134",	"IL17RD",	"MIR326",	"MIR100",
                  "MIR487B", "MIRLET7D", "CDH11",	"NR4A1",	"AGTR1", 	"MIR30A",	"HOXB13")



EmtGenes <- expand.grid(OncoGenesEMT, TumorSuppEMT)
EmtGenes <- as.matrix(EmtGenes)
save(EmtGenes, file = "/Volumes/Macintosh/Dropbox (MechPred)/MechPred/USER/Mohamed/MechanisticModels/Genes/EMT_New.rda")

