tissues = c("brain","breast", "colon","esophagus","kidney","liver", "luad","ovary"
            "prostate", "skin")
library(ranger)
nThreads = 10
ntree = 500
dir.create("data/rdata/RFmodel/",showWarnings=F)

tissue = "luad"

# prepare data for this tissue######
load(paste0("data/procData/traindata_", tissue, ".RData"))
dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
# dat$DNAse_clustered_100bp = NULL
# dat$methylation_100bp = NULL
datchroms = data$muts$chr
chroms = unique(datchroms)
# use only chr1-iteration for this test.
cr = "chr1"
trainData = dat[datchroms != cr,]
testData = dat[datchroms == cr,]
missingTrain = apply(trainData,1, function(x){sum(is.na(x))})
missingTest = apply(testData,1, function(x){sum(is.na(x))})
testData = testData[missingTest == 0,]
#####

# method 1: remove rows with missing data #####
tempData = trainData[missingTrain == 0,]
rf1 = ranger(mutated ~ ., data = tempData, importance='none',num.trees=ntree,
            write.forest = T, seed = 234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',probability = T, verbose=F)
pred1 = predict(object=rf1, data=testData)
rm(tempData)
#####

# method 2: impute with random values sampled from same column #####
nImpute = 5
rfList = lapply(1:nImpute, FUN = function(i){
   cat(i,' ')
   tempData = lapply(1:ncol(trainData),FUN=function(j){ 
      x = trainData[,j]
      missingX = is.na(x)
      x[missingX] = sample(x[!missingX],size = sum(missingX), replace = T)
      return(x) 
   })
   names(tempData) = colnames(trainData)
   tempData = as.data.frame(tempData)
   tempRF = ranger(mutated ~ ., data = tempData, importance='none',
                   num.trees=ceiling(ntree/nImpute),
                   write.forest = T, seed = 234, num.threads =  nThreads,
                   respect.unordered.factors = 'partition',probability = T, verbose=F)
})
combineRanger = function(rfList) {
   res = rfList[[1]]
   res$num.trees = sum(sapply(rfList,function(x){x$num.trees}))
   res$forest$num.trees = sum(sapply(rfList,function(x){x$forest$num.trees}))
   res$call$num.trees = res$num.trees
   res$forest$child.nodeIDs = sapply(rfList,function(x){x$forest$child.nodeIDs})
   res$forest$split.varIDs = sapply(rfList,function(x){x$forest$split.varIDs})
   res$forest$split.values = sapply(rfList,function(x){x$forest$split.values})
   if (!is.null(res$forest$terminal.class.counts)) {
      res$forest$terminal.class.counts = sapply(rfList,function(x){x$forest$terminal.class.counts})
   }
   res
}
rf2 = combineRanger(rfList)
pred2 = predict(object=rf2, data=testData)
rm(rfList)
#####

# method 3: impute with 0 #####
tempData = trainData
tempData[is.na(tempData)] = 0
rf3 = ranger(mutated ~ ., data = tempData, importance='none',num.trees=ntree,
             write.forest = T, seed = 234, num.threads =  nThreads,
             respect.unordered.factors = 'partition',probability = T, verbose=F)
pred3 = predict(object=rf3, data=testData)
#####

# method 4: impute with mean #####
tempData = lapply(1:ncol(trainData),FUN=function(j){ 
   x = trainData[,j]
   if(!is.numeric(x)){return(x)}
   x[is.na(x)] = mean(x,na.rm = T)
   return(x) 
})
names(tempData) = colnames(trainData)
tempData = as.data.frame(tempData)
rf4 = ranger(mutated ~ ., data = tempData, importance='none',num.trees=ntree,
             write.forest = T, seed = 234, num.threads =  nThreads,
             respect.unordered.factors = 'partition',probability = T, verbose=F)
pred4 = predict(object=rf4, data=testData)
#####


# compare Performances #####
library(ROCR)
perfs = lapply(1:4, function(i){
   p = get(paste0("pred", i))$prediction[,2]
   perf = prediction(predictions=p, labels = testData$mutated)
   roc = performance(perf, "tpr", "fpr")
   roc = cbind("TPR" = roc@x.values[[1]], "FPR" = roc@y.values[[1]])
   pr = performance(perf, "prec", "rec")
   pr = cbind("Recall" = pr@x.values[[1]], "Precision" = pr@y.values[[1]])
   
   auc = performance(perf, "auc")@y.values[[1]]
   return(list(roc = roc, pr = pr, auc = auc))
})
# plot ROC
sapply(1:4, function(i){
   if(i== 1){
      plot(perfs[[i]]$roc, type = "l", col = i, lty = i, lwd = 2)
   }else{
      lines(perfs[[i]]$roc, type = "l", col = i, lty = i, lwd = 2)
   }
})
legend("bottomright", col = 1:4, legend=c("remove", "impute", "zero", "mean"),
       lty = 1:4, lwd = 2)
sapply(1:4, function(i){
   if(i== 1){
      plot(perfs[[i]]$pr, type = "l", col = i, lty = i, xlim = c(0,0.2), lwd = 2)
   }else{
      lines(perfs[[i]]$pr, type = "l", col = i, lty = i, lwd = 2)
   }
})
legend("bottomright", col = 1:4, legend=c("remove", "impute", "zero", "mean"),
       lty = 1:4, lwd = 2)

aucs = sapply(perfs, function(x){x$auc})
barplot(aucs, names.arg=c("remove", "impute", "zero", "mean"), col = 1:4)
#####

