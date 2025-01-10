
library(ranger)
nThreads = 16
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "luad","ovary",
            "prostate", "skin")

dumpVar = sapply(tissues, function(tissue){
   # load RF for this tissue
   load(paste0("data/rdata/RFmodel/", tissue, "_", "chr1.RData"))
   
   # iterate through 4 parts of chr1
   for(i in 1:4){
      cat(i,' ')
      # load data
      load(paste0("data/procData/Chr1arm/Chr1arm_allPos_part",i,"_",
                  tissue,"_processed.RData")) #dat
      load(paste0("data/procData/Chr1arm/Chr1arm_allPos_part",i,"_", 
                  tissue, ".RData")) #data
      toExclude = data$pred$ConsensusExcludable == 1 | 
         data$pred$repeatMasker ==1 | 
         data$pred$tandemRepeatFinder == 1
      pos = data$muts[!toExclude,]
      # fix column name mismatch
      colnames(dat)[colnames(dat) == "aPhased_repeats"] = "aPhased_repeates"
      # predict
      p = predict(rf, data = dat, num.threads=nThreads, verbose=F)
      # save
      predictions = data.frame(pos,  prediction = p$predictions[,2])
      save(predictions, file = paste0("data/procData/Chr1arm/Chr1arm_predictions_part",i,
                                      "_", tissue, ".RData"))
   }
   return(NA)
})


print("done")
