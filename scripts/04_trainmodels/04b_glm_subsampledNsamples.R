args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]

# prepare data for this tissue######
load(paste0("data/procData/traindata/traindata_processed_",
            tissue, "_subsampledNsamples.RData"))
chroms = unique(datchroms)
#####

# CWCV glm #####
print("CWCV")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  
  print("training model")
  trainData = dat[datchroms != cr,]
  logR = glm(formula = mutated ~ ., data = trainData, 
             family = binomial(link = "logit"))
  save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, ".RData"))
  
  print("p-values")
  pvals = coef(summary(logR))[,4][-1]
  sigFreatures = names(pvals[pvals < 0.05])
  
  print("training significant models")
  trainData = cbind(trainData[sigFreatures], mutated = trainData$mutated)
  logR = glm(formula = mutated ~ ., data = trainData, 
             family = binomial(link = "logit"))
  save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, "_sig.RData"))
  
  print("predictions on left out chromosome")
  testData = cbind(dat[datchroms == cr,sigFreatures], mutated = dat[datchroms == cr,]$mutated)
  yhat = predict(logR, newdata = testData, type = "response")
  temp = data.frame(pred  = yhat,label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/rdata/GLMmodel/", tissue,
                                "_predictions_subsampledNsamples_sig.RData"))
cat('\n')

# List all CWCV significant predictors
feature_selection = sapply(chroms, function(cr){
  cat(cr, ' ')
  load(paste0("data/rdata/GLMmodel/", tissue, "_", cr, "_sig.RData"))
  return(sigFreatures)
})
names(feature_selection) = chroms
save(feature_selection, file = paste0("data/rdata/GLMmodel/", tissue,"_CWCV_subsampledNsamples_feature_selection.RData"))
cat('\n')
#####


# Final glm #####
print("training final glm")
logR = glm(formula = mutated ~ ., data = dat, family = binomial(link = "logit"))
save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples.RData"))

print("p-values")
pvals = coef(summary(logR))[,4][-1]
sigFreatures = names(pvals[pvals < 0.05])
# check if sig. features are also present in at least 60% of the selected CWCV features
feature_selection = table(do.call(c,feature_selection))
feature_selection = feature_selection[feature_selection > length(chroms)*0.6]
sigFreatures_selection = sigFreatures[sigFreatures %in% names(feature_selection)]

print("training final model on significant features")
trainData = cbind(dat[sigFreatures_selection], mutated = dat$mutated)
logR = glm(formula = mutated ~ ., data = trainData, family = binomial(link = "logit"))
save(logR, sigFreatures, file = paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_sig.RData"))
#####

print("done with this tissue")



# Old model calculation #####
#for(tissue in tissues){
#  print(tissue)
#  # load data for this tissue
#  load(paste0("data/procData/traindata/traindata_processed_",
#              tissue, "_subsampledNsamples.RData"))
#  chroms = unique(datchroms)
#  
#  
#  # glm 
#  print("training models")
#  temp = lapply(chroms, function(cr){
#    cat(cr, ' ')
#    trainData = dat[datchroms != cr,]
#    logR = glm(formula = mutated ~ ., data = trainData, 
#               family = binomial(link = "logit"))
#    save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, ".RData"))
#    return(NULL)
#  })
#  cat('\n')
#  
#  
#  # get variable importances 
#  print("var Importances")
#  imp = sapply(chroms, function(cr){
#    cat(cr, ' ')
#    load(paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, ".RData"))
#    return(logR$coefficients)
#  })
#  save(imp, file = paste0("data/rdata/GLMmodel/", tissue,"_importances_subsampledNsamples.RData"))
#  cat('\n')
#  
#  # predict on held-out chromosomes 
#  print("predictions")
#  predictions = lapply(chroms, function(cr){
#    cat(cr, ' ')
#    load(paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, ".RData"))
#    testData = dat[datchroms == cr,]
#    yhat = predict(logR, newdata = testData, type = "response")
#    temp = data.frame(pred  = yhat,label = testData$mutated)
#    return(temp)
#  })
#  names(predictions) = chroms
#  save(predictions, file = paste0("data/rdata/GLMmodel/", tissue,
#                                  "_predictions_subsampledNsamples.RData"))
#  
#  
#  # calculate variable  p-values 
#  print("pvals")
#  pvals <- sapply(chroms, function(cr){
#    cat(cr, ' ')
#    load(paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, ".RData"))
#    temp = drop1(object = logR, test = "LRT")
#    drop1_features <- setNames(temp$`Pr(>Chi)`, rownames(temp))
#    return(drop1_features)
#  })
#  pvals <- pvals[-1,] # remove first empty row
#  pvals <- as.data.frame(pvals)
#  save(pvals, file = paste0("data/rdata/GLMmodel/", tissue,"_pvals_subsampledNsamples.RData"))
#  
#  # Remove non-significant features and retrain glm 
#  meanPvals = rowMeans(pvals)
#  sigFreatures = names(meanPvals[meanPvals < 0.05])
#  dat = cbind(dat[sigFreatures], mutated = dat$mutated)
#  
#  
#  # glm 
#  print("training significant models")
#  temp = lapply(chroms, function(cr){
#    cat(cr, ' ')
#    trainData = dat[datchroms != cr,]
#    logR = glm(formula = mutated ~ ., data = trainData, 
#               family = binomial(link = "logit"))
#    save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, "_sig.RData"))
#    return(NULL)
#  })
#  cat('\n')
#  
#  # get variable importances 
#  print("var Importances")
#  imp = sapply(chroms, function(cr){
#    cat(cr, ' ')
#    load(paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, "_sig.RData"))
#    return(logR$coefficients)
#  })
#  save(imp, file = paste0("data/rdata/GLMmodel/", tissue,"_importances_subsampledNsamples_sig.RData"))
#  cat('\n')
#  
#  # predict on held-out chromosomes 
#  print("predictions")
#  predictions = lapply(chroms, function(cr){
#    cat(cr, ' ')
#    load(paste0("data/rdata/GLMmodel/", tissue, "_subsampledNsamples_", cr, "_sig.RData"))
#    testData = dat[datchroms == cr,]
#    yhat = predict(logR, newdata = testData, type = "response")
#    temp = data.frame(pred  = yhat,label = testData$mutated)
#    return(temp)
#  })
#  names(predictions) = chroms
#  save(predictions, file = paste0("data/rdata/GLMmodel/", tissue,
#                                  "_predictions_subsampledNsamples_sig.RData"))
#  cat('\n')
#  
#  print("done with this tissue")
#}
#####
