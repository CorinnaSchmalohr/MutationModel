args = as.numeric(commandArgs(trailingOnly=T))
library(ranger)
library(plyr)
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
nThreads = 12

data_size = c("subsampled", "full")
size = data_size[args]
print(size)

# Merge all downsampled or full tissue data
all_dat = sapply(tissues, function(tissue){
  if(size == "subsampled"){
    load(paste0("data/procData/traindata/traindata_processed_",tissue, "_subsampled.RData"))
  } else {
    load(paste0("data/procData/traindata/traindata_processed_",tissue, ".RData"))
  }
  dat = cbind(dat, datchroms)
  return(as.data.frame(dat))
})
all_dat = do.call(rbind.fill, all_dat)
all_dat = all_dat[ ,colSums(is.na(all_dat)) == 0] # removes all cols with NAs

all_datchroms = all_dat$datchroms
all_dat = all_dat[,1:ncol(all_dat)-1]
if(size == "subsampled"){
  save(all_dat,all_datchroms, file=paste0("data/procData/traindata/traindata_processed_subsampled_allTissues.RData"))
} else {
  save(all_dat,all_datchroms, file=paste0("data/procData/traindata/traindata_processed_allTissues.RData"))
}

#load(file=paste0("data/procData/traindata/traindata_processed_allTissues.RData"))
chroms = unique(all_datchroms)

# Train on all data with CWCV #####
# grow forest with impurity_corrected
print("growing forests with impurity_corrected")
temp = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = all_dat[all_datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData, importance = "impurity_corrected",
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=T)
  if(size == "subsampled"){
    save(rf, file = paste0("data/rdata/RFmodel/allTissues_",cr, "_subsampled.RData"))
  } else {
    save(rf, file = paste0("data/rdata/RFmodel/allTissues_",cr, ".RData"))
  }
  return(NULL)
})

# and extract these importances
print("extracting importances")
imp = sapply(chroms, function(cr){
  cat(cr, ' ')
  load(paste0("data/rdata/RFmodel/allTissues_",cr, ".RData"))
  return(rf$variable.importance)
})
if(size == "subsampled"){
  save(imp, file = paste0("data/rdata/RFmodel/allTissues_importances_subsampled_gini.RData"))
} else {
  save(imp, file = paste0("data/rdata/RFmodel/allTissues_importances_gini.RData"))
}



# calculate p-values based on whole data with impurity_corrected #####
#print("compute rf on whole data")
#rf = ranger(mutated ~ ., data = all_dat, importance='impurity_corrected',
#            write.forest = T, seed = 1234, 
#            respect.unordered.factors = 'partition',
#            probability = T, verbose=T)
#save(rf, file = paste0("data/rdata/RFmodel/allTissues_wholeDataRF.RData"))
#print("permutations")
#nPermutations = 10000
#pvals = importance_pvalues(rf, method="altmann", num.permutations=nPermutations,
#                           mutated ~ ., data = all_dat,
#                           seed = 1234,
#                           respect.unordered.factors = 'partition',
#                           scale.permutation.importance = T,
#                           probability = T, verbose=T)
#save(pvals, file = paste0("data/rdata/RFmodel/allTissues_pVals.RData"))
#print("done")
#####

######


# predict on held-out chromosomes #####
print("predictions")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  if(size == "subsampled"){
    load(paste0("data/rdata/RFmodel/allTissues_",cr, "_subsampled.RData"))
  } else {
    load(paste0("data/rdata/RFmodel/allTissues_",cr, ".RData"))
  }
  testData = all_dat[all_datchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
if(size == "subsampled"){
  save(predictions, file = paste0("data/rdata/RFmodel/allTissues_predictions_subsampled.RData"))
} else {
  save(predictions, file = paste0("data/rdata/RFmodel/allTissues_predictions.RData"))
}
cat('\n')

# predict on single-tissue data set chromosome-wise ####
print("predictions on tissues")
for(t in tissues){
  print(t)
  load(paste0("data/procData/traindata/traindata_processed_",t, ".RData"))
  dat = cbind(dat[colnames(dat) %in% colnames(all_dat)], mutated = dat$mutated)
  
  predictions = lapply(chroms, function(cr){
    cat(cr, ' ')
    if(size == "subsampled"){
      load(paste0("data/rdata/RFmodel/allTissues_",cr, "_subsampled.RData"))
    } else {
      load(paste0("data/rdata/RFmodel/allTissues_",cr, ".RData"))
    }
    testData = dat[datchroms == cr,]
    p = predict(rf, data = testData, verbose=F)
    temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
  })
  names(predictions) = chroms
  if(size == "subsampled"){
    save(predictions, file = paste0("data/rdata/RFmodel/allTissues_predictions_on_", t, "_subsampled.RData"))
  } else {
    save(predictions, file = paste0("data/rdata/RFmodel/allTissues_predictions_on_", t, ".RData"))
  }
}

cat('\n')
#####


print("done")
