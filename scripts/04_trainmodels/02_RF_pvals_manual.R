args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]
print(tissue)
library(ranger)
nThreads = 24
nPermutations = 10000
maxData = 30000
dir.create("data/rdata/RFmodel/permutations",showWarnings=F)
dir.create(paste0("data/rdata/RFmodel/permutations/", tissue),showWarnings=F)

# prepare data #####
load(paste0(paste0("data/procData/traindata/traindata_processed_",
                   tissue, ".RData")))
if(nrow(dat) > maxData){
   samp = sample(1:nrow(dat), maxData, replace=F)
   dat = dat[samp,]
   save(samp, dat,file = paste0("data/rdata/RFmodel/",
                                tissue, "_subDataForPvals.RData"))
}
# chroms = unique(datchroms)
#####


# calculate true importance values ####
print("true vals")
rf = ranger(mutated ~ ., data = dat, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads = nThreads,
            respect.unordered.factors = 'partition',
            scale.permutation.importance = T, probability = T, verbose=F)
save(rf, file = paste0("data/rdata/RFmodel/", tissue,
                       "_RFonSubDataForPvals.RData"))
realval = rf$variable.importance
######


# perform permutations #####
perms = sapply(1:nPermutations, function(i){
   cat(i,' ')
   permData = dat
   permData[,"mutated"] = sample(permData[,"mutated"])
   imp = ranger(mutated ~ ., data = permData, importance = 'permutation',
                write.forest = F,  num.threads = nThreads,
                respect.unordered.factors = 'partition',
                verbose = F,scale.permutation.importance = T,
                probability = T)$variable.importance
   save(imp, file=paste0("data/rdata/RFmodel/permutations/", tissue, 
                         "/perm_",i,".RData"))
   return(imp)
})
######


# calculate p-values #####
imp_pvals <- sapply(rownames(perms), function(i) {
   (sum(perms[i, ] >= realval[i]))/(ncol(perms))
})
imp_pvals[imp_pvals==0] = 1/nPermutations
#####

