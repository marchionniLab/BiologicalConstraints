
rm(list = ls())
gc()

## Run XGB

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
library(genefilter)
library(xgboost)

## ---------------------

l1 = load("../Objs/list.data.rda")
l2 = load("../Objs/list.mech.rda")


## ---------------------------------------------
## run XGB at the gene level

xgb_param <- list(
  # General Parameters
  booster            = "gbtree",          # default = "gbtree"
  silent             = 1,                 # default = 0
  # Booster Parameters
  eta                = 0.1,           #0.1   # default = 0.3, range: [0,1]
  gamma              = 0,             #0  # default = 0,   range: [0,∞]
  max_depth          = 1,             # 1
  min_child_weight   = 1,             #1   # default = 1,   range: [0,∞]
  subsample          = 0.4,           #0.4     # default = 1,   range: (0,1]
  colsample_bytree   = 1,           #1      # default = 1,   range: (0,1]
  colsample_bylevel  = 1,    #1  # default = 1,   range: (0,1]
  lambda             = 0,             #0  # default = 1
  alpha              = 0,           # 0    # default = 0
  # Task Parameters
  objective          = "binary:logistic",   # default = "reg:linear"
  eval_metric        = "auc"
)

run_XGB_genes = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, filter=FALSE, classes=NULL){
  
  if(is.null(classes)){
    classes = levels(as.factor(trainGroup))
  }
  
  trainGroup.bin = factor(trainGroup, levels=classes, labels=0:1)
  testGroup.bin = factor(testGroup, levels=classes, labels=0:1)
  
  trainGroup.bin = as.numeric(as.character(trainGroup.bin))
  testGroup.bin = as.numeric(as.character(testGroup.bin))
  
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
  
  idsplit = createDataPartition(trainGroup.bin, p=0.7, list=F)[, 1]
  
  dat.train.xgb = xgb.DMatrix(t(trainMat[, idsplit]), label=trainGroup.bin[idsplit])
  dat.validate.xgb = xgb.DMatrix(t(trainMat[, -idsplit]), label=trainGroup.bin[-idsplit])
  dat.test.xgb = xgb.DMatrix(t(testMat), label=testGroup.bin)
  
  watchlist = list(train=dat.train.xgb, test=dat.validate.xgb)
  
  fit.xgb = xgb.train(xgb_param, data=dat.train.xgb, watchlist=watchlist, 
                      nrounds = 500, early_stopping_rounds = 50, 
                      scale_pos_weight = sum(trainGroup.bin == 0)/sum(trainGroup.bin == 1))

  nvar = -1
  tryCatch({
    
    imp = xgb.importance(model=fit.xgb)
    imp = imp[order(imp$Gain, decreasing=TRUE), ]
    imp = imp[imp$Gain > 0, ]
    nvar = length(imp$Gain)
    
  }, error = function(e){ })

  pred_train = predict(fit.xgb, dat.train.xgb)
  pred_test = predict(fit.xgb, dat.test.xgb)
  
  roc.train = roc(trainGroup[idsplit], pred_train, 
                  plot = F, print.auc = TRUE, 
                  levels = classes, 
                  direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  roc.test = roc(testGroup, pred_test, 
                 plot = F, print.auc = TRUE, 
                 levels = classes, 
                 direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  
  list(mechTSPs=mechTSPs, mechGenes=mechGenes, idsplit=idsplit, fit=fit.xgb, 
       pred.train=pred_train, pred.test=pred_test,
       roc.train=roc.train, roc.test=roc.test, nvar=nvar)
  
}


## ---------------------------------------------
## run XGB at the pair level

run_XGB_pairs = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, filter=FALSE, classes=NULL){
  
  if(is.null(classes)){
    classes = levels(as.factor(trainGroup))
  }
  
  trainGroup.bin = factor(trainGroup, levels=classes, labels=0:1)
  testGroup.bin = factor(testGroup, levels=classes, labels=0:1)
  
  trainGroup.bin = as.numeric(as.character(trainGroup.bin))
  testGroup.bin = as.numeric(as.character(testGroup.bin))
  
  genes = intersect(rownames(trainMat), rownames(testMat))
  
  trainMat = trainMat[genes, ]
  testMat = testMat[genes, ]
  
  mechTSPs = mechTSPs[which(mechTSPs[, 1] %in% genes & mechTSPs[, 2] %in% genes), ]
  
  if(! nrow(mechTSPs) > 0){
    stop("ERROR: no mechanistic pairs found in data")
  }
  
  print(length(unique(as.vector(mechTSPs))))
  
  # if(length(unique(as.vector(mechTSPs))) < 500){
  
    fit.ktsp = SWAP.Train.KTSP(inputMat=trainMat, 
                               phenoGroup=trainGroup, 
                               krange=2:min(nrow(mechTSPs), 25), 
                               FilterFunc=NULL, 
                               featureNo=length(unique(as.vector(mechTSPs))), 
                               RestrictedPairs=mechTSPs)
  # }else{
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

  idsplit = createDataPartition(trainGroup.bin, p=0.7, list=F)[, 1]
  
  dat.train.xgb = xgb.DMatrix(trainV[idsplit, ], label=trainGroup.bin[idsplit])
  dat.validate.xgb = xgb.DMatrix(trainV[-idsplit, ], label=trainGroup.bin[-idsplit])
  dat.test.xgb = xgb.DMatrix(testV, label=testGroup.bin)
  
  watchlist = list(train=dat.train.xgb, test=dat.validate.xgb)
  
  fit.xgb = xgb.train(xgb_param, data=dat.train.xgb, watchlist=watchlist, 
                      nrounds = 500, early_stopping_rounds = 50, 
                      scale_pos_weight = sum(trainGroup.bin == 0)/sum(trainGroup.bin == 1))
  
  
  nvar = -1
  tryCatch({
    
    imp = xgb.importance(model=fit.xgb)
    imp = imp[order(imp$Gain, decreasing=TRUE), ]
    imp = imp[imp$Gain > 0, ]
    nvar = length(imp$Gain)
    
  }, error = function(e){ })
  
  pred_train = predict(fit.xgb, dat.train.xgb)
  pred_test = predict(fit.xgb, dat.test.xgb)
  
  roc.train = roc(trainGroup[idsplit], pred_train, 
                  plot = F, print.auc = TRUE, 
                  levels = classes, 
                  direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  roc.test = roc(testGroup, pred_test, 
                 plot = F, print.auc = TRUE, 
                 levels = classes, 
                 direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  
  list(mechTSPs=mechTSPs, fit.ktsp=fit.ktsp, idsplit=idsplit, trainV=trainV, testV=testV, 
       fit=fit.xgb, pred.train=pred_train, pred.test=pred_test,
       roc.train=roc.train, roc.test=roc.test, nvar=nvar)
  
}

## ---------------------------------------------

## gene level run

list.R.genes =  utils.lapply_i(list.data, function(x, i, data_title){
  
  print(sprintf("======= %s ==========", data_title))
  
  utils.lapply_i(list.mech, function(mechTSPs, j, mech_title){
    
    print(sprintf("[%s]", mech_title))
    
    run_XGB_genes(trainMat=x$trainMat, 
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
      auc_test=round(as.numeric(R$roc.test$auc), 3),
      variables=R$nvar
    )
    
  }))
  
  df = data.frame(names(list.mech), as.data.frame(results))
  colnames(df)[1] = c("mechanism")
  rownames(df) = 1:nrow(df)
  
  df
  
})

## pair level run

list.R.pairs = utils.lapply_i(list.data, function(x, i, data_title){
  
  print(sprintf("======= %s ==========", data_title))
  
  utils.lapply_i(list.mech, function(mechTSPs, j, mech_title){
    
    print(sprintf("[%s]", mech_title))
    
    run_XGB_pairs(trainMat=x$trainMat, 
                     testMat=x$testMat, 
                     trainGroup=x$trainGroup, 
                     testGroup=x$testGroup, 
                     mechTSPs=mechTSPs)
    
  })
  
})


list.run.pairs = utils.lapply_i(list.R.pairs, function(Rlist, i, data_title){
  
  results = Reduce(rbind, utils.lapply_i(Rlist, function(R, j, mech_title){
    
    c(
      candidatePairs=nrow(R$mechTSPs),
      selPairs=nrow(R$fit.ktsp$TSPs),
      auc_train=round(as.numeric(R$roc.train$auc), 3),
      auc_test=round(as.numeric(R$roc.test$auc), 3),
      variables=R$nvar
    )
    
  }))
  
  df = data.frame(names(list.mech), as.data.frame(results))
  colnames(df)[1] = c("mechanism")
  rownames(df) = 1:nrow(df)
  
  df
  
})

## save

save(list.R.genes, list.R.pairs, file="../Objs/list.R.xgb.2.rda")
save(list.run.genes, list.run.pairs, file="../Objs/list.run.xgb.2.rda")









