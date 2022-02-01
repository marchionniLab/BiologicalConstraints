
sink("0.log.txt", split=TRUE)

# source("1_AssembleMech.R")

rm(list = ls())
gc()

## =======================================

print(sprintf("---------- %s -----------", "KTSP"))

tryCatch({

source("3_RunKTSP.R")

}, error = function(e){print(e)})
  
rm(list = ls())
gc()


## =======================================

# print(sprintf("---------- %s -----------", "RF"))
# 
# tryCatch({
# 
# source("4_RunRF.R")
# 
# }, error = function(e){print(e)})
# 
# rm(list = ls())
# gc()


## =======================================

print(sprintf("---------- %s -----------", "SVM"))

tryCatch({

source("5_RunSVM.R")

}, error = function(e){print(e)})
  
rm(list = ls())
gc()


## =======================================


## =======================================

# print(sprintf("---------- %s -----------", "XGB"))
# 
# tryCatch({
#   
#   source("6_RunXGB.R")
#   
# }, error = function(e){print(e)})
# 
# rm(list = ls())
# gc()


## =======================================

sink()



