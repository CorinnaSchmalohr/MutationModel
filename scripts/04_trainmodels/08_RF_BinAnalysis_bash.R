args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("liver", "colon", "skin", "brain","breast", "esophagus","kidney", "lung","ovary","prostate")
t = tissues[args]
#t = "liver"
print(t)
#dir.create("/cellnet/MutationModel/data/rdata/RFmodel/BinAnalysis/", showWarnings = F)


library(ranger)
# Range to the left and right around the pos
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, "1kb" = 500, "100bp" = 50) 

lapply(seq(length(ranges)), function(r){
  
  load(file=paste0("/cellnet/MutationModel/data/procData/traindata/traindata_processed_",
                   t, "_BinAnalysis_",names(ranges[r]),".RData"))
  chroms = unique(datchroms_bin)
  
  
  ### Generate random forrest with permutation importance ##
  #print("growing forests")
  #temp = lapply(chroms, function(cr){
  #  cat(cr, ' ')
  #  trainData = dat_bin[datchroms_bin != cr,]
  #  rf = ranger(mutated ~ ., data = trainData, importance = 'permutation',
  #              write.forest = T, seed = 1234,
  #              respect.unordered.factors = 'partition',
  #              scale.permutation.importance = T, 
  #              probability = T, verbose=F)
  #  save(rf, file = paste0("/cellnet/MutationModel/data/rdata/RFmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
  #  return(NULL)
  #})
  #cat('\n')
  #
  #
  #
  ## get permutation importances ####
  #print("var Importances")
  #imp = sapply(chroms, function(cr){
  #  cat(cr, ' ')
  #  load(paste0("/cellnet/MutationModel/data/rdata/RFmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
  #  return(rf$variable.importance)
  #})
  #save(imp, file = paste0("/cellnet/MutationModel/data/rdata/RFmodel/BinAnalysis/", t, "_",names(ranges[r]),"_importances.RData"))

  
  # grow forest with impurity_corrected #
  print("growing forests with impurity_corrected")
  temp = lapply(chroms, function(cr){
    cat(cr, ' ')
    trainData = dat_bin[datchroms_bin != cr,]
    rf = ranger(mutated ~ ., data = trainData, importance = 'impurity_corrected',
                write.forest = T, seed = 1234, 
                respect.unordered.factors = 'partition',
                probability = T, verbose=T)
    save(rf, file = paste0("/cellnet/MutationModel/data/rdata/RFmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
    return(NULL)
  })
  
  # and extract these importances
  print("extracting importances")
  imp = sapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("/cellnet/MutationModel/data/rdata/RFmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
    return(rf$variable.importance)
  })
  save(imp, file = paste0("/cellnet/MutationModel/data/rdata/RFmodel/BinAnalysis/", t, "_",names(ranges[r]),"_importances.RData"))
  
})

