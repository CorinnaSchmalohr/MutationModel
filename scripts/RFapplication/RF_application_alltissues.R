# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[args]
print(tissue)

load(paste0("data/rdata/", tissue, "/completeData.RData"))

data$precedingDimer = NULL
data$followingDimer = NULL
data$trimer = NULL
data$pentamer = NULL
data$septamer = NULL
data$context = NULL
data$mutated = as.factor(data$mutated)
data$inexon = NULL

library(ranger)
# random Forest #####
chroms = unique(removed$chr)
dir.create(paste0("fig/", tissue, "/RF_variable_importance/"), 
           showWarnings=F)
dir.create(paste0("data/rdata/", tissue, "/RFmodel/"), 
           showWarnings=F)
perf = lapply(chroms, function(cr){
   print(cr)
   rf = ranger(mutated ~ ., data = data, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  8,
               holdout = T, case.weights = as.integer(removed$chr != cr),
               respect.unordered.factors = 'partition', 
               scale.permutation.importance = T, probability = T)
   save(rf, file = paste0("data/rdata/", tissue, "/RFmodel/", cr, ".RData"))
   
   return(c(rf$prediction.error))
})
save(perf, file = paste0("data/rdata/", tissue, "/RFmodel/perf.RData"))
#####



imp = sapply(chroms, function(cr){
   print(cr)
   load(paste0("data/rdata/", tissue, "/RFmodel/", cr, ".RData"))
   return(rf$variable.importance)
})
save(imp, file = paste0("data/rdata/", tissue, "/RFmodel/imp.RData"))

