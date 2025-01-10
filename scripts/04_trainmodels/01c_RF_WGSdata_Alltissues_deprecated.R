.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)
nThreads = 28

# load data for this tissue######
load("data/MutTables/WholeGenomeData/WGSMuts_allTissues_mapped_processed.RData")
chroms = unique(datchroms)
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
save(imp, file = paste0("data/Modeling/WholeGenomeData/RF/allTissues_importances_gini.RData"))
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
  save(rf, file = paste0("data/Modeling/WholeGenomeData/RF/alltissues_", 
                         cr, "_forPrediction.RData"))
  testData = dat[datchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
save(predictions, file = "data/Modeling/WholeGenomeData/RF/allTissues_predictions.RData")
cat('\n')
#####


# create final model #####
print("get final model importance")
rf = ranger(mutated ~ ., data = dat, importance = "impurity_corrected",
            write.forest = F, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
importance = rf$variable.importance
print("get final model for prediction")
rf = ranger(mutated ~ ., data = dat,
            write.forest = T, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
save(rf, importance, file = "data/Modeling/WholeGenomeData/RF/allTissues_finalModel.RData")
#####

print("done")
