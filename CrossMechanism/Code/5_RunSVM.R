
rm(list = ls())
gc()

## Run SVM

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

## ---------------------

l1 = load("../Objs/list.data.rda")
l2 = load("../Objs/list.mech.rda")


## ---------------------------------------------
## run SVM at the gene level

run_SVM_genes = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, filter=FALSE, classes=NULL){
  
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
  
  # grid.svm = expand.grid(degree = 3, scale = 0.01, C = 0.25)
  
  control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)
  
  set.seed(333)
  fit.svmPoly <- train(group ~ ., data=data.frame(group=trainGroup, t(trainMat)), method="svmPoly", trControl=control, metric = "ROC")  
  
  grid.svm = expand.grid(degree = fit.svmPoly$finalModel@kernelf@kpar$degree, 
                         scale = fit.svmPoly$finalModel@kernelf@kpar$scale, 
                         C = fit.svmPoly$finalModel@param$C)
  
  
  fit.svm = train(group ~ ., data=data.frame(group=trainGroup, t(trainMat)), method="svmPoly", 
                  trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), 
                  tuneGrid = grid.svm, metric = "ROC")
  
  nvar = -1
  tryCatch({
    isvm = varImp(fit.svm, scale = TRUE)
    isvm = isvm$importance
    isvm = isvm[order(isvm$Progression, decreasing = TRUE),]
    isvm = isvm[isvm$Progression > 0, ]
    nvar <- nrow(isvm)
  }, error = function(e){ })

  pred_train = predict(fit.svm, newdata=data.frame(t(trainMat)), type="prob")
  pred_test = predict(fit.svm, newdata=data.frame(t(testMat)), type="prob")
  
  roc.train = roc(trainGroup, pred_train[, 2], 
                  plot = F, print.auc = TRUE, 
                  levels = classes, 
                  direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  roc.test = roc(testGroup, pred_test[, 2], 
                 plot = F, print.auc = TRUE, 
                 levels = classes, 
                 direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  
  list(mechTSPs=mechTSPs, mechGenes=mechGenes, fit=fit.svm, pred.train=pred_train, pred.test=pred_test,
       roc.train=roc.train, roc.test=roc.test, nvar=nvar)
  
}


## ---------------------------------------------
## run SVM at the pair level

# x = list.data$Bladder
# data_title = "Bladder"
# trainMat=x$trainMat; testMat=x$testMat; trainGroup=x$trainGroup; testGroup=x$testGroup
# mechTSPs = list.mech$TF_MIR

run_SVM_pairs = function(trainMat, testMat, trainGroup, testGroup, mechTSPs, filter=FALSE, classes=NULL){
  
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
  
  fit.ktsp = SWAP.Train.KTSP(inputMat=trainMat, 
                             phenoGroup=trainGroup, 
                             krange=2:min(nrow(mechTSPs), 25), 
                             FilterFunc=NULL, 
                             featureNo=length(unique(as.vector(mechTSPs))), 
                             RestrictedPairs=mechTSPs)

  trainV = 1 * SWAP.KTSP.Statistics(trainMat, fit.ktsp)$comparisons ## samples x pairs
  testV = 1 * SWAP.KTSP.Statistics(testMat, fit.ktsp)$comparisons

  control <- trainControl(method="repeatedcv", number=10, repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary, allowParallel = T)
  
  set.seed(333)
  fit.svmPoly <- train(group ~ ., data=data.frame(group=trainGroup, trainV), method="svmPoly", trControl=control, metric = "ROC")  
  
  grid.svm = expand.grid(degree = fit.svmPoly$finalModel@kernelf@kpar$degree, 
                         scale = fit.svmPoly$finalModel@kernelf@kpar$scale, 
                         C = fit.svmPoly$finalModel@param$C)
  
  fit.svm = train(group ~ ., data=data.frame(group=trainGroup, trainV), method="svmPoly", 
                  trControl=trainControl(method = "none", classProbs = TRUE, summaryFunction = twoClassSummary), 
                  tuneGrid = grid.svm, metric = "ROC")
  
  nvar = -1
  tryCatch({
    isvm = varImp(fit.svm, scale = TRUE)
    isvm = isvm$importance
    isvm = isvm[order(isvm$Progression, decreasing = TRUE),]
    isvm = isvm[isvm$Progression > 0, ]
    nvar <- nrow(isvm)
  }, error = function(e){  })
  
  pred_train = predict(fit.svm, newdata=data.frame(trainV), type="prob")
  pred_test = predict(fit.svm, newdata=data.frame(testV), type="prob")
  
  roc.train = roc(trainGroup, pred_train[, 2], 
                  plot = F, print.auc = TRUE, 
                  levels = classes, 
                  direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  roc.test = roc(testGroup, pred_test[, 2], 
                 plot = F, print.auc = TRUE, 
                 levels = classes, 
                 direction = "<", col = "blue", lwd = 2, grid = TRUE, auc = TRUE, ci = TRUE)
  
  list(mechTSPs=mechTSPs, fit.ktsp=fit.ktsp, trainV=trainV, testV=testV, fit=fit.svm, pred.train=pred_train, pred.test=pred_test,
       roc.train=roc.train, roc.test=roc.test, nvar=nvar)
  
}

## ---------------------------------------------

## gene level run

list.R.genes =  utils.lapply_i(list.data, function(x, i, data_title){
  
  print(sprintf("======= %s ==========", data_title))
  
  utils.lapply_i(list.mech, function(mechTSPs, j, mech_title){
    
    print(sprintf("[%s]", mech_title))
    
    run_SVM_genes(trainMat=x$trainMat, 
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
      # variables=R$nvar
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
    
    run_SVM_pairs(trainMat=x$trainMat, 
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
      auc_test=round(as.numeric(R$roc.test$auc), 3)
    )
    
  }))
  
  df = data.frame(names(list.mech), as.data.frame(results))
  colnames(df)[1] = c("mechanism")
  rownames(df) = 1:nrow(df)
  
  df
  
})

## save

save(list.R.genes, list.R.pairs, file="../Objs/list.R.svm.2.rda")
save(list.run.genes, list.run.pairs, file="../Objs/list.run.svm.2.rda")









