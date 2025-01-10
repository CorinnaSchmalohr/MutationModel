args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("liver", "colon", "skin", "brain","breast", "esophagus","kidney", "luad","ovary","prostate")
t = tissues[args]
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, "1kb" = 500, "100bp" = 50) # Range to the left and right around the pos
#dir.create(path = "./data/rdata/GLMmodel/BinAnalysis/",showWarnings=F)


print(t)

lapply(seq(length(ranges), 1), function(r){
  
  load(file=paste0("data/procData/traindata/traindata_processed_",
                   t, "_BinAnalysis_",names(ranges[r]),".RData"))
  chroms = unique(datchroms_bin)
  
  # glm #####
  print("training models")
  temp = lapply(chroms, function(cr){
    cat(cr, ' ')
    trainData = dat_bin[datchroms_bin != cr,]
    logR = glm(formula = mutated ~ ., data = trainData, 
               family = binomial(link = "logit"))
    save(logR, file = paste0("data/rdata/GLMmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
    return(NULL)
  })
  cat('\n')
  #####
  
  # calculate variable  p-values ####
  print("pvals")
  pvals <- sapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("data/rdata/GLMmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData")) 
    temp = drop1(object = logR, test = "LRT")
    drop1_features <- setNames(temp$`Pr(>Chi)`, rownames(temp))
    return(drop1_features)
  })
  pvals <- pvals[-1,] # remove first empty row
  pvals <- as.data.frame(pvals)
  save(pvals, file = paste0("data/rdata/GLMmodel/", t, "_",names(ranges[r]), "_pvals.RData"))
  #####
  
  # get variable importances ####
  print("var Importances")
  imp = sapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("data/rdata/GLMmodel/BinAnalysis/", t, "_", names(ranges[r]), "_", cr, ".RData"))
    return(logR$coefficients)
  })
  save(imp, file = paste0("data/rdata/GLMmodel/", t, "_",names(ranges[r]),"_importances.RData"))
  cat('\n')
  #####
})
