# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
tissue = tissues[args]
print(tissue)
dir.create(paste0("data/rdata/", tissue, "/permutations"), showWarnings = F)

nPermutations = 10000
maxData = 50000
nThreads = 24
load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))

data$context = NULL
data$pentamer = NULL
data$trimer = NULL
data$septamer = NULL
data$inexon = NULL

library(ranger)

# random Forest on entire data  #####
print("pvals")
rf = ranger(mutated ~ ., data = data, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads = nThreads,
            respect.unordered.factors = 'partition',
            scale.permutation.importance = T, probability = T)
save(rf, file = paste0("data/rdata/", tissue, 
                       "/RFmodel/allChrs_withBinWise.RData"))
realval = rf$variable.importance


if(nrow(data) > maxData){
   samp = sample(1:nrow(data), maxData, replace=F)
   dataSubset = data[samp,]
   rfSubset = ranger(mutated ~ ., data = dataSubset, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads = nThreads,
               respect.unordered.factors = 'partition',
               scale.permutation.importance = T, probability = T)
   save(rfSubset, dataSubset, samp, file = paste0("data/rdata/", tissue, 
                                            "/RFmodel/allChrs_withBinWise_Subset.RData"))
   realval = rfSubset$variable.importance
}


perms = sapply(1:nPermutations, function(i){
   cat(i,' ')
   if(nrow(data) > maxData){
      permData  = dataSubset
   } else{
      permData = data
   }
   permData[,"mutated"] = sample(permData[,"mutated"])
   imp = ranger(mutated ~ ., data = permData, importance = 'permutation',
                write.forest = F,  num.threads = nThreads,
                respect.unordered.factors = 'partition',
                verbose = F,scale.permutation.importance = T,
                probability = T)$variable.importance
   save(imp, file=paste0("data/rdata/", tissue, 
                         "/permutations/perm_",i,".RData"))
   return(imp)
})

imp_pvals <- sapply(rownames(perms), function(i) {
   (sum(perms[i, ] >= realval[i]))/(ncol(perms))
})
imp_pvals[imp_pvals==0] = 1/nPermutations
save(imp_pvals, file = paste0("data/rdata/", tissue, 
                          "/RFmodel/imp_pvals_withBinWise_manual.RData"))
#####
