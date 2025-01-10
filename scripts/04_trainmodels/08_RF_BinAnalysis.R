library(ranger)
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, "1kb" = 500, "100bp" = 50) # Range to the left and right around the pos

for(t in tissues){
  print(t)
  
  lapply(seq(length(ranges), 1), function(r){
  
    load(file=paste0("data/procData/traindata/traindata_processed_",
                     t, "_BinAnalysis_",names(ranges[r]),".RData"))
    chroms = unique(datchroms_bin)
    
    ## Generate random forrest with permutation importance ##
    print("growing forests")
    temp = lapply(chroms, function(cr){
      trainData = dat_bin[datchroms_bin != cr,]
      rf = ranger(mutated ~ ., data = trainData, importance = 'permutation',
                  write.forest = T, seed = 1234,
                  respect.unordered.factors = 'partition',
                  scale.permutation.importance = T, 
                  probability = T, verbose=F)
      save(rf, file = paste0("data/rdata/RFmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
      return(NULL)
    })
    
    
    # get permutation importances ####
    print("var Importances")
    imp = sapply(chroms, function(cr){
      load(paste0("data/rdata/RFmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
      return(rf$variable.importance)
    })
    save(imp, file = paste0("data/rdata/RFmodel/BinAnalysis/", t, "_",names(ranges[r]),
                            "_importances.RData"))
  
  })
  
}

