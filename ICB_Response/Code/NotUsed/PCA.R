
### Prinicpal component analysis on the NP dataset

library(ggbiplot)

# PCA for the samples
x <- t(Expr_NP)
pc <- prcomp(x, center  = T)
names(pc)

plot(pc$x[, 1], pc$x[, 2], col = Pheno_NP$SR, main = "PCA", xlab = "PC1", ylab = "PC2")


ggbiplot(pc, groups = Pheno_NP$SR, ellipse = TRUE, 
         circle = F)

