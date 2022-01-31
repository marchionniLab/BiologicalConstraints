
library(knitr)
library(rutils)

list.boot = list()

# bladder

l1 = load("../../Bladder/Objs/KTSP/bootobjectKTSP_RAND.rda")

list.boot[["KTSP"]] = 
  list(
    list(boot=bootobject_74, n=74),
    list(boot=bootobject_100, n=100),
    list(boot=bootobject_200, n=200),
    list(boot=bootobject_500, n=500)
  )

rm(list=l1)

l2 = load("../../Bladder/Objs/RF/RFBootObjects_RAND.rda")

list.boot[["RF"]] = 
  list(
    list(boot=bootobjectMech, n=0),
    list(boot=bootobjectAgnostic_74, n=74),
    list(boot=bootobjectAgnostic_100, n=100),
    list(boot=bootobjectAgnostic_200, n=200),
    list(boot=bootobjectAgnostic_500, n=500)
  )

rm(list=l2)

l3 = load("../../Bladder/Objs/SVM/SVMBootObjects_RAND.rda")

list.boot[["SVM"]] = 
  list(
    list(boot=bootobjectMech, n=0),
    list(boot=bootobjectAgnostic_74, n=74),
    list(boot=bootobjectAgnostic_100, n=100),
    list(boot=bootobjectAgnostic_200, n=200),
    list(boot=bootobjectAgnostic_500, n=500)
  )

rm(list=l3)

l4 = load("../../Bladder/Objs/XGB/XGBBootObjects_RAND.rda")

list.boot[["XGB"]] = 
  list(
    list(boot=bootobjectMech, n=0),
    list(boot=bootobjectAgnostic_74, n=74),
    list(boot=bootobjectAgnostic_100, n=100),
    list(boot=bootobjectAgnostic_200, n=200),
    list(boot=bootobjectAgnostic_500, n=500)
  )

rm(list=l4)

# ------------------

df.list = utils.lapply_i(list.boot, function(blist, i, y){
  
  print(y)
  
  if(y == "KTSP"){
    
    df = Reduce(rbind, lapply(blist, function(bobj){
      
      nB = bobj$n
      B = bobj$boot
      
      Bdf = data.frame(B$t)
      colnames(Bdf) = c("rand_n", "rand_auc_train", "rand_auc_test", "mech_n", "mech_auc_train", "mech_auc_test", "rand_diff", "mech_diff")
      
      df = rbind(
        data.frame(auc=Bdf$rand_auc_train, type="random", n=nB, data="train"),
        data.frame(auc=Bdf$mech_auc_train, type="mech", n=nB, data="train"),
        data.frame(auc=Bdf$rand_auc_test, type="random", n=nB, data="test"),
        data.frame(auc=Bdf$mech_auc_test, type="mech", n=nB, data="test")
      )
      
      df
      
    }))
    
    df
    
  }else{
    
    df = Reduce(rbind, lapply(blist, function(bobj){
      
      nB = bobj$n
      B = bobj$boot
      
      Bdf = data.frame(B$t)
      colnames(Bdf) = c("auc_train", "auc_test", "n")

      # df = rbind(
      #   data.frame(auc=Bdf$auc_train, type=ifelse(nB==0, "mech", "random"), n=Bdf$n, data="train"),
      #   data.frame(auc=Bdf$auc_test, type=ifelse(nB==0, "mech", "random"), n=Bdf$n, data="test")
      # )
      
            
      df = rbind(
        data.frame(auc=Bdf$auc_train, type=ifelse(nB==0, "mech", "random"), n=nB, data="train"),
        data.frame(auc=Bdf$auc_test, type=ifelse(nB==0, "mech", "random"), n=nB, data="test")
      )
      
      df
      
    }))
    
    df
    
  }
  
  print(table(df$n, df$type))
  
  df$n = factor(df$n)
  df$data = factor(df$data, levels=c("train", "test"))
  df$type = factor(df$type, levels=c("random", "mech"))
  
  df$set = paste(df$type, ", n=", df$n, sep="")
  df$set[which(df$type == "mech")] = "mech"
  df$set = factor(df$set)
  
  df
  
})

save(list.boot, df.list, file="rand_boot_bladder.rda")


