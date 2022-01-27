
rm(list = ls())
gc()

## Run RF

############################################################################
### Load library
require(switchBox)
require(Biobase)
require(limma)
require(pROC)
require(caret)
require(RColorBrewer)
require(ggplot2)
require(reshape)
require(plotROC)
library(enrichR)
library(mltools)
library(xtable)
library(boot)
library(patchwork)
library(randomForest)

## ---------------------

l1 = load("../Objs/list.data.rda")
l2 = load("../Objs/list.mech.rda")


run_RF = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, filter=FALSE, classes=NULL){
  
  if(is.null(classes)){
    classes = levels(as.factor(trainGroup))
  }
  
  genes = intersect(rownames(trainMat), rownames(testMat))
  
  trainMat = trainMat[genes, ]
  testMat = testMat[genes, ]
  
  mechTSPs = mechTSPs[which(mechTSPs[, 1] %in% genes & mechTSPs[, 2] %in% genes), ]
  
  if(! nrow(mechTSPs) > 0){
    stop("ERROR: no mechanistic pairs found in data")
  }
  
  mechGenes = intersect(genes, unique(as.vector(mechTSPs)))
  
  trainMat = trainMat[genes, ]
  testMat = testMat[genes, ]
  
  tmp <- as.vector(table(trainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  
  fit.rf = tuneRF(x=t(trainMat), 
                  y=trainGroup, 
                  mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, 
                  sampsize = sampsizes)
  
  irf = randomForest::importance(fit.rf, scale=FALSE, type = 2)
  irf = irf[order(irf, decreasing=TRUE), ]
  irf = irf[irf > 0]
  
  nvar = length(irf)
  
  pred_train = predict(fit.rf, newdata=t(trainMat), type="vote")
  pred_test = predict(fit.rf, newdata=t(testMat), type="vote")

  roc.train = roc(trainGroup, pred_train[, 2], 
                  plot = F, print.auc = TRUE, 
                  levels = classes, 
                  direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  roc.test = roc(testGroup, pred_test[, 2], 
                 plot = F, print.auc = TRUE, 
                 levels = classes, 
                 direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  
  list(mechTSPs=mechTSPs, mechGenes=mechGenes, fit=fit.rf, pred.train=pred_train, pred.test=pred_test,
       roc.train=roc.train, roc.test=roc.test, nvar=nvar)
  
}

list.run = lapply(list.data, function(x){
  
  results = Reduce(rbind, lapply(list.mech, function(mechTSPs){
    
    R = run_RF(trainMat=x$trainMat, 
                 testMat=x$testMat, 
                 trainGroup=x$trainGroup, 
                 testGroup=x$testGroup, 
                 mechTSPs=mechTSPs)
    
    c(
      candidateGenes=length(R$mechGenes),
      auc_train=round(as.numeric(R$roc.train$auc), 3),
      auc_test=round(as.numeric(R$roc.test$auc), 3)
    )
    
  }))
  
  df = data.frame(names(list.mech), as.data.frame(results))
  colnames(df)[1] = c("mechanism")
  rownames(df) = 1:nrow(df)
  
  df
  
})


save(list.run, file="../Objs/list.run.rf.rda")









