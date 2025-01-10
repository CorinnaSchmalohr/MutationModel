# get mutation rate and nPositions, percTP, correlation etc. (dataInfos) #####
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
dataInfos = sapply(tissues, function(tissue){
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue, ".RData"))
  chrCount = sapply(table(datchroms), as.integer)
  chrPerc = sapply(split(dat$mutated, datchroms),function(x){mean(x==1)})
  cors = cor(dat[sapply(dat, is.numeric)], use = "pair")
  return(list(nMuts = nrow(dat), 
              percTP = mean(dat$mutated == 1), 
              nMutsPerChr = chrCount,
              percTPperChr = chrPerc, 
              cors = cors))
}, simplify = F)
save(dataInfos, file = "data/rdata/dataInfos.RData")
#####