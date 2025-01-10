.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)


dir.create("data/Modeling/WholeExomeData/RF", showWarnings = F)
load(paste0("data/Modeling/exomeTrainData/RF/allTissues_finalModel.RData"))

nThreads = 28
dumpVar = sapply(1:50, function(i){
   print(i)
   # load data
   load(paste0("data/MutTables/WholeExomeData/exomeMuts_part", i, "_allTissues.RData"))
   testDat = data$pred
   yhat = predict(rf, testDat, type = "response", num.threads = nThreads)
   
   predictions = data.frame(datpos,  
                            prediction = do.call(c,predictions))
   save(predictions, 
        file = paste0("data/Modeling/WholeExomeData/RF/exomeMuts_RFpredictions_part",
                      i, ".RData"))
   return(NA)
})
print("done")
