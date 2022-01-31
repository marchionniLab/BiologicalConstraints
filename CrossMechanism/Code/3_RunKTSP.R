
rm(list = ls())
gc()

## Run KTSP

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
library(rutils)

## ---------------------

l1 = load("../Objs/list.data.rda")
l2 = load("../Objs/list.mech.rda")


run_ktsp = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, krange, featNo, filter=FALSE, classes=NULL){
  
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
  
  if(filter){
    filter_fun = SWAP.Filter.Wilcoxon
  }else{
    filter_fun = NULL
  }
  
  ktsp.fit = SWAP.Train.KTSP(inputMat=trainMat, 
                              phenoGroup=trainGroup, 
                              krange=krange, 
                              FilterFunc=filter_fun, 
                              featureNo=featNo, 
                              RestrictedPairs=mechTSPs) 
  
  stats.train = SWAP.KTSP.Statistics(inputMat=trainMat, classifier=ktsp.fit, CombineFunc=sum)
  stats.test = SWAP.KTSP.Statistics(inputMat=testMat, classifier=ktsp.fit, CombineFunc=sum)
  
  roc.train = roc(trainGroup, stats.train$statistics, 
                  plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", 
                  levels = classes, direction = "<", col="blue", lwd=2, grid=TRUE)
  roc.test = roc(testGroup, stats.test$statistics, 
                  plot = F, ci = T, print.thres = "all", print.auc = TRUE, print.auc.col="black", 
                  levels = classes, direction = "<", col="blue", lwd=2, grid=TRUE)
  
  list(mechTSPs=mechTSPs, ktsp.fit=ktsp.fit, stats.train=stats.train, stats.test=stats.test,
       roc.train=roc.train, roc.test=roc.test)
  
}

list.R = utils.lapply_i(list.data, function(x, i, data_title){
  
  print(sprintf("======= %s ==========", data_title))
  
  utils.lapply_i(list.mech, function(mechTSPs, j, mech_title){
    
    print(sprintf("[%s]", mech_title))
    
    run_ktsp(trainMat=x$trainMat, 
                 testMat=x$testMat, 
                 trainGroup=x$trainGroup, 
                 testGroup=x$testGroup, 
                 mechTSPs=mechTSPs, 
                 krange=1:50, featNo=100)
    
  })
  
})

list.run = lapply(list.R, function(Rlist){
  
  results = Reduce(rbind, lapply(Rlist, function(R){
    
    c(
      candidatePairs=nrow(R$mechTSPs),
      selPairs=nrow(R$ktsp.fit$TSPs),
      genes=length(unique(as.vector(R$ktsp.fit$TSPs))),
      auc_train=round(as.numeric(R$roc.train$auc), 3),
      auc_test=round(as.numeric(R$roc.test$auc), 3)
    )
    
  }))
  
  df = data.frame(names(list.mech), as.data.frame(results))
  colnames(df)[1] = c("mechanism")
  rownames(df) = 1:nrow(df)
  
  df
  
})

save(list.R, file="../Objs/list.R.ktsp.rda")
save(list.run, file="../Objs/list.run.ktsp.rda")









