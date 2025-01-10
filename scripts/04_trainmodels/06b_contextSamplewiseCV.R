# preparation #####
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
dir.create("data/rdata/contextPatientCV", showWarnings=F)
load("data/rdata/pentamers.RData")
# load baseline pentamer frequencies per chromosome 
load("data/rdata/chrFreqs.RData")
######


# iterate through tissues and to patient-wise CV
dumpVar = sapply(tissues, function(tissue){
   print(tissue)
   # load data for this tissue
   load(paste0("data/rdata/", tissue, "/MutsWithTumorIDs.RData"))

   # take the 50 samples with the most muts and use them for CV
   samps = table(MutsWIDs$TumorID)
   samps = names(sort(samps, decreasing=T)[1:50])

   # get context-based predictions
   pChrs = sapply(samps, function(samp){
      cat(samp, ' ')
      TPsel = which(MutsWIDs$TumorID == samp)
      selContext = table(MutsWIDs$context[TPsel])
      TNsel = do.call(c,sapply(1:length(selContext), function(i){
         cont = names(selContext)[i]
         temp = which(MutsWIDs$context == cont & MutsWIDs$mutated == 0)
         sample(temp,
                size=min(selContext[i], length(temp)))
      }, simplify=F))
      sel = c(TPsel, TNsel)
      trainDat = MutsWIDs[-sel,]
      # testDat = MutsWIDs[sel,]
      
      # train Context
      pChr = sapply(paste0("chr", 1:22), function(cr){
         sub = trainDat$context[trainDat$chr != cr & trainDat$mutated == 1]
         fr1 = table(factor(sub, levels=pents[,1]))
         fr2 = table(factor(sub, levels=pents[,2]))
         fr = as.numeric(fr1+fr2)
         res = fr/(rowSums(chrFreqs) - chrFreqs[,cr])
         return(res)
      })
      return(pChr)
      # # apply to testData
      # indx = cbind(pent2context[testDat$context], testDat$chr)
      # predictions = pChr[indx]
      # 
      # # return predictions
      # temp = data.frame(pred = predictions,  label = testDat$mutated)
      # return(temp)
   }, simplify = F)
   cat('\n')
   save(pChrs, file = paste0("data/rdata/contextPatientCV/", tissue,
                                  "_pChrs.RData"))
   return(NA)
})
print("done")
#####

