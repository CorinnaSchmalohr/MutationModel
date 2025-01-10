args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "luad", "ovary", "prostate", "skin")
t = tissues[args]
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "none" = 0) # Range to the left and right around the pos
#dir.create(path = "./data/rdata/GLMmodel/BinAnalysis/",showWarnings=F)
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")


print(t)

lapply(seq(length(ranges), 1), function(r){
  print(r)
  load(file=paste0("data/procData/traindata/traindata_processed_",
                   t, "_BinAnalysis_",names(ranges[r]),".RData"))
  chroms = unique(datchroms_bin)
  featureN = colnames(dat_bin)[colnames(dat_bin) != "mutated"]
  # glm #####
  print("training models")
  temp = lapply(chroms, function(cr){
    trainData = dat_bin[datchroms_bin != cr,]
    logR = glm(formula = mutated ~ ., data = trainData, 
               family = binomial(link = "logit"))
    save(logR, file = paste0("data/rdata/GLMmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
    return(NULL)
  })
  #####
  
  # calculate variable  p-values ####
  print("pvals")
  pvals <- sapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("data/rdata/GLMmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData")) 
    p = summary(logR)$coefficients[,"Pr(>|z|)"]
    return(p[featureN])
  })
  pvals <- as.data.frame(pvals)
  rownames(pvals) = featureN
  save(pvals, file = paste0("data/rdata/GLMmodel/BinAnalysis/", t, "_",names(ranges[r]), "_pvals.RData"))
  #####
  
  # get variable importances ####
  print("var Importances")
  imp = sapply(chroms, function(cr){
    load(paste0("data/rdata/GLMmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
    return(logR$coefficients)
  })
  save(imp, file = paste0("data/rdata/GLMmodel/", t, "_",names(ranges[r]),"_importances.RData"))
  #####
})

