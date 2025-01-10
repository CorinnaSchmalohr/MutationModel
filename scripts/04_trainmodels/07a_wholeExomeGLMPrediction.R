args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
#tissue = tissues[args]
#library(ranger)
#nThreads = 16
for(tissue in tissues){
  print(tissue)
  dumpVar = sapply(1:31, function(i){
     print(i)
     # load data
     load(paste0("data/procData/exomeMuts/exomeMuts_part",i,"_",
                 tissue,"_processed.RData"))
     chrs = unique(datpos$chr)
     
     # Final model prediction 
     load(paste0("data/rdata/GLMmodel/", tissue, "_sig.RData"))
     yhat = predict(logR, newdata = dat, type = "response")
     predictions = data.frame(datpos, prediction = yhat)
     
     
     # CWCV prediction with RF #####
     #predictions = sapply(chrs, function(cr){
     #   load(paste0("data/rdata/RFmodel/", tissue, "_", cr,".RData"))
     #   subdata = dat[datpos$chr == cr,]
     #   p = predict(rf, data = subdata, num.threads=nThreads, verbose=F)
     #   return(p$predictions[,2])
     #}, simplify = F)  
     #predictions = data.frame(datpos,  
     #                         prediction = do.call(c,predictions))
     #save(predictions, 
     #     file = paste0("data/procData/exomeMuts/exomeMuts_RFpredictions_part",
     #                   i, "_", tissue, ".RData"))
     #####
     
     # CWCV prediction with glm #####
     #predictions = sapply(chrs, function(cr){
     #  load(paste0("data/rdata/GLMmodel/", tissue, "_", cr,"_sig.RData"))
     #  subdata = dat[datpos == cr,]
     #  yhat = predict(logR, newdata = subdata, type = "response")
     #  return(yhat)
     #}, simplify = F)  
     #predictions = data.frame(datpos,  
     #                         prediction = do.call(c,predictions))
     #####
    
     save(predictions, 
          file = paste0("data/procData/exomeMuts/exomeMuts_GLMpredictions_part",
                        i, "_", tissue, ".RData"))
     return(NA)
     return(NA)
  })
  print("done")
}

