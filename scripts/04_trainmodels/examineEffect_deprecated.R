tissue = "luad"

# 0 = NONCODING or NA
# 1 = SYNONYMOUS 
# 2 = NONSYNONYMOUS
# 3 = START-LOST, STOP-GAIN or STOP-LOSS

# examine mutation effect
load(paste0("data/procData/traindata_", tissue, ".RData"))
dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
datchroms = data$muts$chr
toExclude = dat$ConsensusExcludable == 1 | dat$repeatMasker ==1 | dat$tandemRepeatFinder == 1
dat = dat[!toExclude,]
dat$ConsensusExcludable = NULL
dat$repeatMasker = NULL
dat$tandemRepeatFinder = NULL
datchroms = datchroms[!toExclude]
chroms = unique(datchroms)

temp = split(round(dat$effect, digits = 2), dat$mutated)
lvls = sort(unique(round(dat$effect, digits = 2)))
temp = sapply(temp, function(x){table(factor(x, levels=lvls))})
barplot(t(temp), beside=T, las = 1, xlab = "effect score", ylab = "count")
legend("topright", legend=c("TN", "TP"), fill=c("grey20", "grey80") )
#####


# test RF with and without effect score #####
load(paste0("data/procData/traindata_processed_",
       tissue, ".RData"))
library(ROCR)
library(ranger)
chroms = unique(datchroms)
nThreads = 10
res1 = lapply(chroms, function(cr){
   cat(cr, ' ')
   trainData = dat[datchroms != cr,]
   rf = ranger(mutated ~ ., data = trainData, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  nThreads,
               respect.unordered.factors = 'partition',
               scale.permutation.importance = T, probability = T, verbose=F)
   testData = dat[datchroms == cr,]
   p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
   ROCRpredict = prediction(predictions=p$predictions[,2], labels=testData$mutated)
   auc = performance(ROCRpredict,"auc")@y.values[[1]]
   roc = performance(ROCRpredict,"tpr", "fpr")
   return(perf)
})
res2 = lapply(chroms, function(cr){
   cat(cr, ' ')
   trainData = dat[datchroms != cr,]
   trainData$effect = NULL
   rf = ranger(mutated ~ ., data = trainData, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  nThreads,
               respect.unordered.factors = 'partition',
               scale.permutation.importance = T, probability = T, verbose=F)
   testData = dat[datchroms == cr,]
   p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
   ROCRpredict = prediction(predictions=p$predictions[,2], labels=testData$mutated)
   perf = performance(ROCRpredict,"auc")
   return(perf)
})
rocs1 = sapply(res1, function(x){x@y.values[[1]]})
rocs2 = sapply(res2, function(x){x@y.values[[1]]})
boxplot(cbind(rocs1, rocs2), las = 1, ylab = "AUROC", 
        names = c("with effect score", "without effect score"))
######