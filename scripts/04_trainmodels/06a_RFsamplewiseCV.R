args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "luad","ovary",
            "prostate", "skin")
tissue = tissues[args]
print(tissue)
#library(ranger)
#nThreads = 16
dir.create("data/rdata/GLMmodelPatientCV", showWarnings=F)


# load data for this tissue######
load(paste0("data/procData/traindata/traindata_processed_",
            tissue, ".RData"))
load(paste0("data/rdata/", tissue, "/MutsWithTumorIDs.RData"))
load(paste0("data/procData/traindata/traindata_", tissue, ".RData"))
# ids1 = paste(MutsWIDs$chr, MutsWIDs$pos, sep = "_")
# ids2 = paste(data$muts$chr, data$muts$pos, sep = "_")
# all.equal(ids1, ids2)
toExclude = data$pred$ConsensusExcludable == 1 | data$pred$repeatMasker == 1 | data$pred$tandemRepeatFinder == 1 | 
  data$pred$mappability_100mer < 1 |data$pred$mappabiliy_40mer < 1 | data$pred$mappability_24mer < 1
MutsWIDs = MutsWIDs[!toExclude,]
rm(data, toExclude)
#####


# take the 50 with the most muts and use them for CV #####
samps = table(MutsWIDs$TumorID)
samps = names(sort(samps, decreasing=T)[1:50])
#####



# glm #####
print("training models")
GLMresults = sapply(samps, function(samp){
  cat(samp, ' ')
  # train glm
  TPsel = which(MutsWIDs$TumorID == samp)
  selContext = table(MutsWIDs$context[TPsel])
  TNsel = do.call(c,sapply(1:length(selContext), function(i){
       cont = names(selContext)[i]
       temp = which(MutsWIDs$context == cont & MutsWIDs$mutated == 0)
       sample(temp, size=min(selContext[i], length(temp)))
  }, simplify=F))
  sel = c(TPsel, TNsel)
  trainData = dat[-sel,]
  logR = glm(formula = mutated ~ ., data = trainData, 
              family = binomial(link = "logit"))

  # calculate variable  p-values ####
  print("pvals")
  temp = drop1(object = logR, test = "LRT")
  drop1_features <- setNames(temp$`Pr(>Chi)`, rownames(temp))
  pvals <- drop1_features[-1] # remove first empty row
  pvals <- as.data.frame(pvals)
  rm(temp)
  cat('\n')
  
  # Remove non-significant features and retrain glm #####
  meanPvals = rowMeans(pvals)
  sigFreatures = names(meanPvals[meanPvals < 0.05])
  trainData = cbind(trainData[sigFreatures], mutated = trainData$mutated)
  
  # glm ##
  print("training significant models")
  logR = glm(formula = mutated ~ ., data = trainData, 
             family = binomial(link = "logit"))
  save(logR, file = paste0("data/rdata/GLMmodelPatientCV/", tissue, "_",samp, "_sig.RData"))
  cat('\n')
  
  # predict on held-out chromosomes #
  print("significant predictions")
  testData = dat[sel,]
  yhat = predict(logR, newdata = testData, type = "response")
  temp = data.frame(pred = yhat,label = testData$mutated)

  return(list(pvals = pvals, predictions = temp))
})
cat('\n')
save(GLMresults, file = paste0("data/rdata/GLMmodelPatientCV/", tissue, "_GLMresults.RData"))

#####



# random Forest with permutation importance #####
#print("growing forests")
#RFresults = sapply(samps, function(samp){
#   cat(samp, ' ')
#   # train RF
#   TPsel = which(MutsWIDs$TumorID == samp)
#   selContext = table(MutsWIDs$context[TPsel])
#   TNsel = do.call(c,sapply(1:length(selContext), function(i){
#      cont = names(selContext)[i]
#      temp = which(MutsWIDs$context == cont & MutsWIDs$mutated == 0)
#      sample(temp,
#             size=min(selContext[i], length(temp)))
#   }, simplify=F))
#   sel = c(TPsel, TNsel)
#   trainData = dat[-sel,]
#   rf = ranger(mutated ~ ., data = trainData, importance = 'permutation',
#               write.forest = T, seed = 1234, num.threads =  nThreads,
#               respect.unordered.factors = 'partition',
#               scale.permutation.importance = T, 
#               probability = T, verbose=F)
#   save(rf, file = paste0("data/rdata/RFmodelPatientCV/", tissue, "_", 
#                          samp, ".RData"))
#   
#   # predict on held-out patient
#   testData = dat[sel,]
#   p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
#   temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
#   # return importance and predictions
#   return(list(importance = rf$variable.importance, predictions = temp))
#}, simplify = F)
#cat('\n')
#save(RFresults, file = paste0("data/rdata/RFmodelPatientCV/", tissue,
#                        "_RFresults.RData"))
######

print("done")
