
# preparation #####
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
# load("data/rdata/context_pPerChr.RData")
load("data/rdata/pentamers.RData")
#####


# predict mutation based on sequence context for whole exome #####
dumpVar = sapply(tissues, function(tissue){
   print(tissue)
   
   # get TP positions per sample CV
   load(paste0("data/rdata/", tissue, "/MutsWithTumorIDs.RData"))

   # take the 50 samples with the most muts and use them for CV
   samps = table(MutsWIDs$TumorID)
   samps = names(sort(samps, decreasing=T)[1:50])
   TPpos = sapply(samps, function(samp){
      sub = MutsWIDs[MutsWIDs$TumorID %in% samp,]
      pos = paste(sub$chr, sub$pos, sep = "_")
      return(pos)
   }, simplify = F)
   #rm(data, toExclude, MutsWIDs, samps)
   rm(MutsWIDs, samps)
   
   # load sample-wise trained context data
   load(paste0("data/rdata/contextPatientCV/", tissue,
               "_pChrs.RData"))
   
   dumpVar2 = sapply(1:31, function(i){
      cat(i, ' ')
      # load data
      #load(paste0("data/procData/exomeMuts/exomeMuts_part",i,"_",
      #            tissue,"_processed.RData")) # dat and datpos
     load(paste0("data/procData/exomeMuts/exomeMuts_part",i,
                 "_", tissue, ".RData"))
     datpos = as.data.frame(cbind(chr = data$muts$chr, context = data$muts$context))

      rm(data)
      # predictions
      indx = cbind(pent2context[datpos$context], datpos$chr)
      allPos = paste(datpos$chr, datpos$pos, sep = "_")
      predictions = sapply(names(pChrs), function(samp){
         pred = data.frame(pred = pChrs[[samp]][indx], 
                           label = allPos %in% TPpos[[samp]])
         save(pred, file = paste0("data/rdata/contextPatientCV/exomeMuts_part",
                                         i,"_", tissue,"_",samp,"_results.RData"))
         return(NA)
      }, simplify = F)
      
      return(NA)
   }, simplify = F)
   cat('\n')
   return(NA)
})
#####
