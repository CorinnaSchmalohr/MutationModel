args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]
library(ranger)
library(ROCR)
library(Biostrings)
print(tissue)
nThreads = 12

chroms = paste0("chr", 1:22)

# prepare trimers and count per tissue
bases = c("A", "C", "G", "T")
trimers = sort(apply(expand.grid(bases, bases, bases),1,paste0, collapse = ""))
compl = cbind(trimers,
              as.character(reverseComplement(DNAStringSet(trimers))))
compl = t(apply(compl,1,sort))
compl = unique(compl)
rownames(compl) = paste(compl[,1], compl[,2], sep = "/")
#####


# load test Data
tissue2predict = tissue
cat("testData: ", tissue2predict, "\n")
load(paste0("data/procData/traindata/traindata_processed_",
            tissue2predict, ".RData"))
testDat = dat
testChroms = datchroms
testTrimers = paste0(testDat$context_precedingBase, 
                     testDat$context_ref,
                     testDat$context_followingBase)
rm(dat, datchroms)

cat("trainData: ")
ROCs = sapply(tissues, function(predictingTissue){
  cat(predictingTissue,' ')
  
  # load train Data
  load(paste0("data/procData/traindata/traindata_processed_",
              predictingTissue, ".RData"))
  trainDat = dat
  trainChroms = datchroms
  rm(dat, datchroms)
  
  # make sure we have same set of predictors
  predictors = colnames(trainDat)
  testDat = testDat[,colnames(testDat) %in% predictors]
  toAdd = predictors[!predictors %in% colnames(testDat)]
  for(x in toAdd){
    testDat[x] = rep(mean(trainDat[,x]), nrow(testDat))
  }
  
  # iterate through chromosomes and get predictions
  predsPerChr = sapply(chroms, function(cr){
    # load  RF
    load(paste0("data/rdata/RFmodel/", predictingTissue,
                "_", cr, ".RData"))
    # subset test Data to chr
    subDat = testDat[testChroms == cr,]
    subTrimers = testTrimers[testChroms == cr]
    # get predictions on chromosome
    yHat = predict(rf,data=subDat,type="response", num.threads=10)
    # return 3 columns: prediction, true value, trimer
    return(data.frame(prediction = yHat$predictions[,2], 
                      trueVal = subDat$mutated, 
                      trimer = subTrimers))
  }, simplify=F)
  predsPerChr = do.call(rbind,predsPerChr)
  
  # get performance per trimer
  ROCperTrimer = apply(compl,1,function(trimer){
    subPreds = predsPerChr[predsPerChr$trimer %in% trimer,]
    temp = prediction(predictions=subPreds$prediction, 
                      labels=subPreds$trueVal)
    auc = performance(temp,"auc")
    return(auc@y.values[[1]])
  })
  return(ROCperTrimer)
})
cat('\n')
save(ROCs, file = paste0("data/rdata/",tissue,"/crossTissue_testTrimer.RData"))
#####

#####
print("subset trainData to trimer")
predictingTissue = tissue
cat("trainData: ", predictingTissue, "\n")
load(paste0("data/procData/traindata/traindata_processed_",
            predictingTissue, ".RData"))
trainDat = dat
trainChroms = datchroms
rm(dat, datchroms)
datTrimer = paste0(trainDat$context_precedingBase, 
                   trainDat$context_ref,
                   trainDat$context_followingBase)
# train a model for only one trimer, excluding chr1 and chr2
dumpTrimers = apply(compl,1,function(trimer){
  # train a model only on one trimer
  subtrainDat = trainDat[datTrimer %in% trimer & 
                           !trainChroms %in% c("chr1", "chr2"),]
  rf = ranger(mutated ~ ., data = subtrainDat, 
              write.forest = T, seed = 1234, num.threads =  10,
              respect.unordered.factors = 'partition',
              probability = T, verbose=F, importance = "impurity_corrected")
  save(rf, file = paste0("data/rdata/RFmodel/", predictingTissue,
                         "_", trimer[1], "_",trimer[2],
                         "_chr1and2.RData"))
  return(NA)
})

# iterate through test tissues and test trimer-specific models on whole data
cat("testData: ")
ROCs = sapply(tissues, function(tissue2predict){
  # load test Data
  cat(tissue2predict,' ')
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue2predict, ".RData"))
  testDat = dat[datchroms %in% c("chr1", "chr2"),]
  rm(dat, datchroms)
  # make sure we have same set of predictors
  predictors = colnames(trainDat)
  testDat = testDat[,colnames(testDat) %in% predictors]
  toAdd = predictors[!predictors %in% colnames(testDat)]
  for(x in toAdd){
    testDat[x] = rep(mean(trainDat[,x]), nrow(testDat))
  }
  # iterate through trimers and get predictions
  ROCsperTrimer = apply(compl,1,function(trimer){
    # train a model only on one trimer
    load(paste0("data/rdata/RFmodel/", predictingTissue,
                "_", trimer[1], "_",trimer[2],
                "_chr1and2.RData"))
    yHat = predict(rf,data=testDat,type="response", num.threads=10)
    temp = prediction(predictions=yHat$predictions[,2], 
                      labels=testDat$mutated)
    auc = performance(temp,"auc")
    return(auc@y.values[[1]])
  })
  return(ROCsperTrimer)
})
cat('\n')
save(ROCs, file = paste0("data/rdata/",tissue,"/crossTissue_trainTrimer.RData"))
######


# both train and test data for specific trimer  #####
predictingTissue = tissue
cat("trainData: ", predictingTissue, "\n")
load(paste0("data/procData/traindata/traindata_processed_",
            predictingTissue, ".RData"))
trainDat = dat
trainChroms = datchroms
rm(dat, datchroms)
datTrimer = paste0(trainDat$context_precedingBase, 
                   trainDat$context_ref,
                   trainDat$context_followingBase)

# iterate through test tissues and test trimer-specific models on whole data
cat("testData: ")
ROCs = sapply(tissues, function(tissue2predict){
  # load test Data
  cat(tissue2predict, " ")
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue2predict, ".RData"))
  testDat = dat[datchroms %in% c("chr1", "chr2"),]
  testTrimers = paste0(testDat$context_precedingBase,
                       testDat$context_ref,
                       testDat$context_followingBase)
  rm(dat, datchroms)
  # make sure we have same set of predictors
  predictors = colnames(trainDat)
  testDat = testDat[,colnames(testDat) %in% predictors]
  toAdd = predictors[!predictors %in% colnames(testDat)]
  for(x in toAdd){
    testDat[x] = rep(mean(trainDat[,x]), nrow(testDat))
  }
  # iterate through trimers and get predictions
  ROCsperTrimer = apply(compl,1,function(trimer){
    subtestDat = testDat[testTrimers %in% trimer,]
    # train a model only on one trimer
    load(paste0("data/rdata/RFmodel/", predictingTissue,
                "_", trimer[1], "_",trimer[2],
                "_chr1and2.RData"))
    yHat = predict(rf,data=subtestDat,type="response", num.threads=10)
    temp = prediction(predictions=yHat$predictions[,2], 
                      labels=subtestDat$mutated)
    auc = performance(temp,"auc")
    return(auc@y.values[[1]])
  })
  return(ROCsperTrimer)
})
cat('\n')
save(ROCs, file = paste0("data/rdata/",tissue,"/crossTissue_bothTrimer.RData"))
######
