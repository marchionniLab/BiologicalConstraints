
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
library(rutils)

## ---------------------

l1 = load("../Objs/list.data.rda")
l2 = load("../Objs/list.mech.rda")


## ---------------------------------------------
## run RF at the gene level

run_RF_genes = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, filter=FALSE, classes=NULL){
  
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


## ---------------------------------------------
## run RF at the pair level

# x = list.data$Bladder
# data_title = "Bladder"
# trainMat=x$trainMat; testMat=x$testMat; trainGroup=x$trainGroup; testGroup=x$testGroup
# mechTSPs = list.mech$TF_MIR

run_RF_pairs = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, filter=FALSE, classes=NULL, data_title){
  
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
  
  print(length(unique(as.vector(mechTSPs))))
  
  if(data_title != "ZBladder"){
    
    fit.ktsp = SWAP.Train.KTSP(inputMat=trainMat, 
                               phenoGroup=trainGroup, 
                               krange=2:min(nrow(mechTSPs), 25), 
                               FilterFunc=NULL, 
                               featureNo=length(unique(as.vector(mechTSPs))), 
                               RestrictedPairs=mechTSPs)
    
  }else{
    
    # bladder
    
    fit.ktsp = SWAP.Train.KTSP(inputMat=trainMat[unique(as.vector(mechTSPs)), ], 
                               phenoGroup=trainGroup, 
                               krange=37, 
                               FilterFunc=SWAP.Filter.Wilcoxon, 
                               featureNo=length(unique(as.vector(mechTSPs))))
    
  }
    
  # if(length(unique(as.vector(mechTSPs))) > 500){
  #   
  #   fit.ktsp = SWAP.Train.KTSP(inputMat=trainMat, 
  #                              phenoGroup=trainGroup, 
  #                              krange=1:min(nrow(mechTSPs), 50), 
  #                              FilterFunc=SWAP.Filter.Wilcoxon, 
  #                              featureNo=500, 
  #                              RestrictedPairs=mechTSPs)
  # 
  # }
  
  trainV = 1 * SWAP.KTSP.Statistics(trainMat, fit.ktsp)$comparisons ## samples x pairs
  testV = 1 * SWAP.KTSP.Statistics(testMat, fit.ktsp)$comparisons
  
  tmp <- as.vector(table(trainGroup))
  num_classes <- length(tmp)
  min_size <- tmp[order(tmp,decreasing=FALSE)[1]]
  sampsizes <- rep(min_size,num_classes)
  
  fit.rf = tuneRF(x=trainV, 
                  y=trainGroup, 
                  mtryStart = 1, ntreeTry=500, stepFactor = 1, improve=0.05, trace=F, plot=F, doBest=T, 
                  sampsize = sampsizes)
  
  irf = randomForest::importance(fit.rf, scale=FALSE, type = 2)
  irf = irf[order(irf, decreasing=TRUE), ]
  irf = irf[irf > 0]
  
  nvar = length(irf)
  
  pred_train = predict(fit.rf, newdata=trainV, type="vote")
  pred_test = predict(fit.rf, newdata=testV, type="vote")
  
  roc.train = roc(trainGroup, pred_train[, 2], 
                  plot = F, print.auc = TRUE, 
                  levels = classes, 
                  direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  roc.test = roc(testGroup, pred_test[, 2], 
                 plot = F, print.auc = TRUE, 
                 levels = classes, 
                 direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  
  list(mechTSPs=mechTSPs, fit.ktsp=fit.ktsp, trainV=trainV, testV=testV, fit=fit.rf, pred.train=pred_train, pred.test=pred_test,
       roc.train=roc.train, roc.test=roc.test, nvar=nvar)
  
}

## ---------------------------------------------

## gene level run

list.R.genes =  utils.lapply_i(list.data, function(x, i, data_title){
  
  print(sprintf("======= %s ==========", data_title))
  
  utils.lapply_i(list.mech, function(mechTSPs, j, mech_title){
    
    print(sprintf("[%s]", mech_title))
    
    run_RF_genes(trainMat=x$trainMat, 
                     testMat=x$testMat, 
                     trainGroup=x$trainGroup, 
                     testGroup=x$testGroup, 
                     mechTSPs=mechTSPs)
    
  })
  
})

list.run.genes = utils.lapply_i(list.R.genes, function(Rlist, i, data_title){
  
  results = Reduce(rbind, utils.lapply_i(Rlist, function(R, j, mech_title){
    
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

## pair level run

print("====== PAIRS =========")

list.R.pairs = utils.lapply_i(list.data, function(x, i, data_title){
  
  print(sprintf("======= %s ==========", data_title))
  
  utils.lapply_i(list.mech, function(mechTSPs, j, mech_title){
    
    print(sprintf("[%s]", mech_title))
    
    run_RF_pairs(trainMat=x$trainMat, 
                     testMat=x$testMat, 
                     trainGroup=x$trainGroup, 
                     testGroup=x$testGroup, 
                     mechTSPs=mechTSPs,
                     data_title=data_title)
    
  })
  
})


list.run.pairs = utils.lapply_i(list.R.pairs, function(Rlist, i, data_title){
  
  results = Reduce(rbind, utils.lapply_i(Rlist, function(R, j, mech_title){
    
    c(
      candidatePairs=nrow(R$mechTSPs),
      selPairs=nrow(R$fit.ktsp$TSPs),
      auc_train=round(as.numeric(R$roc.train$auc), 3),
      auc_test=round(as.numeric(R$roc.test$auc), 3)
    )
    
  }))
  
  df = data.frame(names(list.mech), as.data.frame(results))
  colnames(df)[1] = c("mechanism")
  rownames(df) = 1:nrow(df)
  
  df
  
})

## save

save(list.R.genes, list.R.pairs, file="../Objs/list.R.rf.2.rda") 
save(list.run.genes, list.run.pairs, file="../Objs/list.run.rf.2.rda")









