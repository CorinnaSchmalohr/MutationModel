args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)
nThreads = 16
# nPermutations = 10000
# maxData = 50000
dir.create("data/Modeling/exomeTrainData/",showWarnings=F)
dir.create("data/Modeling/exomeTrainData/RF",showWarnings=F)

print(tissue)

# load data for this tissue######
load(paste0("data/MutTables/exomeTrainData/", 
            tissue, "_Muts_mapped_processed.RData"))
chroms = unique(datchroms)
#####

# grow forest with impurity_corrected #####
print("growing forests with impurity_corrected")
temp = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData, importance = "impurity_corrected",
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=T)
  save(rf, file = paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", 
                         cr, ".RData"))
  return(NULL)
})

# and extract the importances
print("extracting importances")
imp = sapply(chroms, function(cr){
  cat(cr, ' ')
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", cr, ".RData"))
  return(rf$variable.importance)
})
save(imp, file = paste0("data/Modeling/exomeTrainData/RF/", tissue,
                        "_importances_gini.RData"))
######


# predict on held-out chromosomes #####
print("predictions")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData,
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=T)
  save(rf, file = paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", 
                         cr, "_forPrediction.RData"))
  testData = dat[datchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/Modeling/exomeTrainData/RF/", tissue,
                                "_predictions.RData"))
cat('\n')
#####

print("done")
