.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)
nThreads = 56
# nPermutations = 10000
# maxData = 50000
dir.create("data/Modeling/WholeGenomeData/",showWarnings=F)
dir.create("data/Modeling/WholeGenomeData/RF",showWarnings=F)


# load data for this tissue######
load("data/MutTables/exomeTrainData/allTissues_Muts_mapped_processed.RData")
chroms = unique(datchroms)
#####


# create final model #####
rf = ranger(mutated ~ ., data = dat, importance = "impurity_corrected",
            write.forest = F, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
importance = rf$variable.importance
rf = ranger(mutated ~ ., data = dat,
            write.forest = T, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
save(rf, importance, file = "data/Modeling/exomeTrainData/RF/allTissues_finalModel.RData")

#####


# grow forest with impurity_corrected #####
print("growing forests with impurity_corrected")
imp = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData, importance = "impurity_corrected",
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=T)
  return(rf$variable.importance)
})
save(imp, file = paste0("data/Modeling/exomeTrainData/RF/allTissues_importances_gini.RData"))
######


# predict on held-out chromosomes #####
print("predictions")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData,
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=F)
  save(rf, file = paste0("data/Modeling/exomeTrainData/RF/alltissues_", 
                         cr, "_forPrediction.RData"))
  testData = dat[datchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
save(predictions, file = "data/Modeling/exomeTrainData/RF/allTissues_predictions.RData")
cat('\n')
#####




print("done")
