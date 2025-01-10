args = as.numeric(commandArgs(trailingOnly=T))
library(plyr)
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")

data_size = c("subsampled", "full")
size = data_size[args]
print(size)

# Merge all down-sampled or full tissue data #####
all_dat = sapply(tissues, function(tissue){
  if(size == "subsampled"){
    load(paste0("data/procData/traindata/traindata_processed_",tissue, "_subsampled.RData"))
  } else {
    load(paste0("data/procData/traindata/traindata_processed_",tissue, ".RData"))
  }
  dat = cbind(dat, datchroms)
  return(as.data.frame(dat))
})
all_dat = do.call(rbind.fill, all_dat)
all_dat = all_dat[ ,colSums(is.na(all_dat)) == 0] # removes all cols with NAs

all_datchroms = all_dat$datchroms
all_dat = all_dat[,1:ncol(all_dat)-1]
if(size == "subsampled"){
  save(all_dat,all_datchroms, file=paste0("data/procData/traindata/traindata_processed_subsampled_allTissues.RData"))
} else {
  save(all_dat,all_datchroms, file=paste0("data/procData/traindata/traindata_processed_allTissues.RData"))
}
chroms = unique(all_datchroms)
#####

# CWCV glm #####
print("CWCV")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  
  print("training model")
  trainData = all_dat[all_datchroms != cr,]
  logR = glm(formula = mutated ~ ., data = trainData, 
             family = binomial(link = "logit"))
  if(size == "subsampled"){
    save(logR, file = paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled.RData"))
  } else {
    save(logR, file = paste0("data/rdata/GLMmodel/allTissues_",cr, ".RData"))
  }
  
  print("p-values")
  pvals = coef(summary(logR))[,4][-1]
  sigFreatures = names(pvals[pvals < 0.05])
  
  print("training significant models")
  trainData = cbind(trainData[sigFreatures], mutated = trainData$mutated)
  logR = glm(formula = mutated ~ ., data = trainData, 
             family = binomial(link = "logit"))
  if(size == "subsampled"){
    save(logR, sigFreatures, file = paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled_sig.RData"))
  } else {
    save(logR, sigFreatures, file = paste0("data/rdata/GLMmodel/allTissues_",cr, "_sig.RData"))
  }
  
  print("predictions on left out chromosome")
  testData = cbind(all_dat[all_datchroms == cr, sigFreatures], mutated = all_dat[all_datchroms == cr,]$mutated)
  yhat = predict(logR, newdata = testData, type = "response")
  temp = data.frame(pred  = yhat,label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
if(size == "subsampled"){
  save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_subsampled_sig.RData"))
} else {
  save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_sig.RData"))
}
cat('\n')

# predict on single-tissue data set chromosome-wise #####
print("predictions on tissues")
dumpVar = lapply(tissues, function(t){  
  print(t)
  load(paste0("data/procData/traindata/traindata_processed_",t, ".RData"))
  dat = cbind(dat[colnames(dat) %in% colnames(all_dat)], mutated = dat$mutated)
  
  predictions = lapply(chroms, function(cr){
    cat(cr, ' ')
    if(size == "subsampled"){
      load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled_sig.RData"))
    } else {
      load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_sig.RData"))
    }
    
    testData = dat[datchroms == cr,]
    yhat = predict(logR, newdata = testData, type = "response")
    temp = data.frame(pred  = yhat,label = testData$mutated)
    return(temp)
  })
  names(predictions) = chroms
  if(size == "subsampled"){
    save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_on_", t, "_subsampled.RData"))
  } else {
    save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_on_", t, ".RData"))
  }
  return(NULL)
})
cat('\n')
#####

# List all CWCV significant predictors
feature_selection = sapply(chroms, function(cr){
  cat(cr, ' ')
  if(size == "subsampled"){
    load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled_sig.RData"))
  } else {
    load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_sig.RData"))
  }
  return(sigFreatures)
})
names(feature_selection) = chroms
if(size == "subsampled"){
  save(feature_selection, "data/rdata/GLMmodel/allTissues_CWCV_feature_selection_subsampled.RData")
} else {
  save(feature_selection, "data/rdata/GLMmodel/allTissues_CWCV_feature_selection.RData")
}
cat('\n')
#####


# Final glm #####
print("training final glm")
logR = glm(formula = mutated ~ ., data = all_dat, family = binomial(link = "logit"))
if(size == "subsampled"){
  save(logR, file = "data/rdata/GLMmodel/allTissues_subsampled.RData")
} else {
  save(logR, file = "data/rdata/GLMmodel/allTissues.RData")
}
print("p-values")
pvals = coef(summary(logR))[,4][-1]
sigFreatures = names(pvals[pvals < 0.05])
# check if sig. features are also present in at least 60% of the selected CWCV features
feature_selection = table(do.call(c,feature_selection))
feature_selection = feature_selection[feature_selection > length(chroms)*0.6]
sigFreatures_selection = sigFreatures[sigFreatures %in% names(feature_selection)]

print("training final model on significant features")
trainData = cbind(all_dat[sigFreatures_selection], mutated = all_dat$mutated)
logR = glm(formula = mutated ~ ., data = trainData, family = binomial(link = "logit"))
if(size == "subsampled"){
  save(logR, file = "data/rdata/GLMmodel/allTissues_subsampled_sig.RData")
} else {
  save(logR, file = "data/rdata/GLMmodel/allTissues_sig.RData")
}
#####
print("done")











## Train on all data with CWCV #####
## glm #####
#print("training models")
#temp = lapply(chroms, function(cr){
#  cat(cr, ' ')
#  trainData = all_dat[all_datchroms != cr,]
#  logR = glm(formula = mutated ~ ., data = trainData, 
#             family = binomial(link = "logit"))
#  if(size == "subsampled"){
#    save(logR, file = paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled.RData"))
#  } else {
#    save(logR, file = paste0("data/rdata/GLMmodel/allTissues_",cr, ".RData"))
#  }
#  return(NULL)
#})
#cat('\n')
######
#
## calculate variable  p-values ####
#print("pvals")
#pvals <- sapply(chroms, function(cr){
#  cat(cr, ' ')
#  if(size == "subsampled"){
#    load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled.RData"))
#  } else {
#    load(paste0("data/rdata/GLMmodel/allTissues_",cr, ".RData"))
#  }
#  temp = drop1(object = logR, test = "LRT")
#  drop1_features <- setNames(temp$`Pr(>Chi)`, rownames(temp))
#  return(drop1_features)
#})
#pvals <- pvals[-1,] # remove first empty row
#pvals <- as.data.frame(pvals)
#if(size == "subsampled"){
#  save(pvals, file = paste0("data/rdata/GLMmodel/allTissues_pvals_subsampled.RData"))
#} else {
#  save(pvals, file = paste0("data/rdata/GLMmodel/allTissues_pvals.RData"))
#}
######
#
## Remove non-significant features and retrain glm #####
#meanPvals = rowMeans(pvals)
#sigFreatures = names(meanPvals[meanPvals < 0.05])
#all_dat = cbind(all_dat[sigFreatures], mutated = all_dat$mutated)
#
#
## Train on all data with CWCV #####
## glm #####
#print("training model on sig. features")
#temp = lapply(chroms, function(cr){
#  cat(cr, ' ')
#  trainData = all_dat[all_datchroms != cr,]
#  logR = glm(formula = mutated ~ ., data = trainData, 
#             family = binomial(link = "logit"))
#  if(size == "subsampled"){
#    save(logR, file = paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled_sig.RData"))
#  } else {
#    save(logR, file = paste0("data/rdata/GLMmodel/allTissues_",cr, "_sig.RData"))
#  }
#  return(NULL)
#})
#cat('\n')
######
#
## get variable importances ####
#print("Significant var Importances")
#imp = sapply(chroms, function(cr){
#  cat(cr, ' ')
#  if(size == "subsampled"){
#    load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled_sig.RData"))
#  } else {
#    load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_sig.RData"))
#  }
#  return(logR$coefficients)
#})
#if(size == "subsampled"){
#  save(imp, file = paste0("data/rdata/GLMmodel/allTissues_importances_subsampled_sig.RData"))
#} else {
#  save(imp, file = paste0("data/rdata/GLMmodel/allTissues_importances_sig.RData"))
#}
#cat('\n')
######
#
## predict on held-out chromosomes #####
#print("sig. predictions")
#predictions = lapply(chroms, function(cr){
#  cat(cr, ' ')
#  if(size == "subsampled"){
#    load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled_sig.RData"))
#  } else {
#    load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_sig.RData"))
#  }
#  testData = all_dat[all_datchroms == cr,]
#  yhat = predict(logR, newdata = testData, type = "response")
#  temp = data.frame(pred  = yhat,label = testData$mutated)
#  return(temp)
#})
#names(predictions) = chroms
#if(size == "subsampled"){
#  save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_subsampled_sig.RData"))
#} else {
#  save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_sig.RData"))
#}
#cat('\n')
######
#
#
#
## predict on single-tissue data set chromosome-wise ####
#print("predictions on tissues")
#for(t in tissues){
#  print(t)
#  load(paste0("data/procData/traindata/traindata_processed_",t, ".RData"))
#  dat = cbind(dat[colnames(dat) %in% colnames(all_dat)], mutated = dat$mutated)
#  
#  predictions = lapply(chroms, function(cr){
#    cat(cr, ' ')
#    if(size == "subsampled"){
#      load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_subsampled_sig.RData"))
#    } else {
#      load(paste0("data/rdata/GLMmodel/allTissues_",cr, "_sig.RData"))
#    }
#    
#    testData = dat[datchroms == cr,]
#    yhat = predict(logR, newdata = testData, type = "response")
#    temp = data.frame(pred  = yhat,label = testData$mutated)
#    return(temp)
#  })
#  names(predictions) = chroms
#  if(size == "subsampled"){
#    save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_on_", t, "_subsampled.RData"))
#  } else {
#    save(predictions, file = paste0("data/rdata/GLMmodel/allTissues_predictions_on_", t, ".RData"))
#  }
#}
#
#cat('\n')
######
#
#
#print("done")
#
