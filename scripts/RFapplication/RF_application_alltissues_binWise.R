# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[args]
print(tissue)
nThreads = 24
load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))

data$context = NULL
data$pentamer = NULL
data$trimer = NULL
data$septamer = NULL
data$inexon = NULL
data$replDirection = as.integer(factor(as.character(data$replDirection),
                            levels=c("left", "unknown", "right")))

library(ranger)
# random Forest #####
chroms = unique(removed$chr)
dir.create(paste0("fig/", tissue, "/RF_variable_importance/"),
           showWarnings=F)
dir.create(paste0("data/rdata/", tissue, "/RFmodel/"),
           showWarnings=F)
print("growing forests")
temp = lapply(chroms, function(cr){
   cat(cr, ' ')
   subData = data[removed$chr != cr,]
   rf = ranger(mutated ~ ., data = subData, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  nThreads,
               respect.unordered.factors = 'partition',
               scale.permutation.importance = T, probability = T, verbose=F)
   save(rf, file = paste0("data/rdata/", tissue, "/RFmodel/",
                          cr, "_withBinWise.RData"))
   return(NULL)
})
cat('\n')
#####


# get variable importances ####
print("var Importances")
imp = sapply(chroms, function(cr){
   cat(cr, ' ')
   load(paste0("data/rdata/", tissue, "/RFmodel/",
               cr, "_withBinWise.RData"))
   return(rf$variable.importance)
})
save(imp, file = paste0("data/rdata/", tissue,
                        "/RFmodel/imp_withBinWise.RData"))
cat('\n')
#####


# predict on held-out chromosomes #####
print("predictions")
predictions = lapply(chroms, function(cr){
   cat(cr, ' ')
   load(paste0("data/rdata/", tissue, "/RFmodel/", cr, "_withBinWise.RData"))
   p = predict(rf, data = data[removed$chr == cr,], num.threads=nThreads, verbose=F)
   temp = data.frame(p$predictions,data[removed$chr == cr,"mutated"])
   colnames(temp) = c("pred_0", "pred_1", "labels")
   return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/rdata/", tissue,
                                "/RFmodel/predictions_withBinWise.RData"))
cat('\n')
#####


# calculate variable importance p-values ####
print("pvals")
nPermutations = 10000
maxData = 50000
if(nrow(data) > maxData){
   samp = sample(1:nrow(data), maxData, replace=F)
   data = data[samp,]
   save(samp, data,file = paste0("data/rdata/", tissue, 
                            "/RFmodel/allChrs_withBinWise_sample.RData"))
}
print("true vals")
rf = ranger(mutated ~ ., data = data, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads = nThreads,
            respect.unordered.factors = 'partition',
            scale.permutation.importance = T, probability = T, verbose=F)
save(rf, file = paste0("data/rdata/", tissue,
                       "/RFmodel/allChrs_withBinWise.RData"))
print("permutations")
pvals = importance_pvalues(rf, method="altmann", num.permutations=nPermutations,
                           mutated ~ ., data = data,
                           seed = 1234, num.threads = nThreads,
                           respect.unordered.factors = 'partition',
                           scale.permutation.importance = T,
                           probability = T, verbose=T)
save(pvals, file = paste0("data/rdata/", tissue,
                          "/RFmodel/imp_pvals_withBinWise.RData"))
#####

print("done")
