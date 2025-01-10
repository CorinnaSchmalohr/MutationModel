#tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
#            "prostate", "skin")
source("scripts/05_analysis/00_NamesAndColors.R")
library(ranger)
library(ROCR)
library(ggplot2)
cr= "chr1"

# only train data subsampled #####
subTrain = sapply(tissues, function(tissue2predict){
   cat("data: ", tissue2predict, "\n")
   load(paste0("data/procData/traindata/traindata_processed_",
               tissue2predict, ".RData"))
   tempData = dat[datchroms == cr,]
   truePreds = tempData$mutated
   cat("model: ")
   perfs = sapply(tissues, function(predictingTissue){
      cat(predictingTissue, " ")
      load(paste0("data/rdata/RFmodel/", predictingTissue, "_subsampled_", cr, ".RData"))
      load(paste0("data/procData/traindata/traindata_processed_",
                  predictingTissue, "_subsampled.RData"))
      predictors = names(rf$variable.importance)
      testData = tempData[,colnames(tempData) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testData)]
      for(x in toAdd){
         testData[x] = rep(mean(dat[,x]), nrow(testData))
      }
      yHat = predict(rf,data=testData,type="response", num.threads=10)
      temp = prediction(yHat$predictions[,2], truePreds)
      roc = performance(temp, "tpr", "fpr")
      pr = performance(temp,"prec", "rec")
      auc = performance(temp,"auc")
      perf = list(roc = roc, pr = pr, auc = auc)
      return(perf)
   }, simplify=F)
   cat("\n")
   return(perfs)
}, simplify = F)
save(subTrain, file = "data/rdata/CrossTissue_subTrain.RData")
#####

# only test data subsampled #####
subTest = sapply(tissues, function(tissue2predict){
   cat("data: ", tissue2predict, "\n")
   load(paste0("data/procData/traindata/traindata_processed_",
               tissue2predict, "_subsampled.RData"))
   tempData = dat[datchroms == cr,]
   truePreds = tempData$mutated
   cat("model: ")
   perfs = sapply(tissues, function(predictingTissue){
      cat(predictingTissue, " ")
      load(paste0("data/procData/traindata/traindata_processed_",
                  predictingTissue, ".RData"))
      load(paste0("data/rdata/RFmodel/", predictingTissue,
                  "_", cr, ".RData"))
      predictors = names(rf$variable.importance)
      testData = tempData[,colnames(tempData) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testData)]
      for(x in toAdd){
         testData[x] = rep(mean(dat[,x]), nrow(testData))
      }
      yHat = predict(rf,data=testData,type="response", num.threads=10)
      temp = prediction(yHat$predictions[,2], truePreds)
      roc = performance(temp, "tpr", "fpr")
      pr = performance(temp,"prec", "rec")
      auc = performance(temp,"auc")
      perf = list(roc = roc, pr = pr, auc = auc)
      return(perf)
   }, simplify=F)
   cat("\n")
   return(perfs)
}, simplify = F)
save(subTest, file = "data/rdata/CrossTissue_subTest.RData")
#####

# both subsampled #####
subBoth = sapply(tissues, function(tissue2predict){
   cat("data: ", tissue2predict, "\n")
   load(paste0("data/procData/traindata/traindata_processed_",
               tissue2predict, "_subsampled.RData"))
   tempData = dat[datchroms == cr,]
   truePreds = tempData$mutated
   cat("model: ")
   perfs = sapply(tissues, function(predictingTissue){
      cat(predictingTissue, " ")
      load(paste0("data/rdata/RFmodel/", predictingTissue,
                  "_subsampled_", cr, ".RData"))
      load(paste0("data/procData/traindata/traindata_processed_",
                  predictingTissue, "_subsampled.RData"))
      predictors = names(rf$variable.importance)
      testData = tempData[,colnames(tempData) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testData)]
      for(x in toAdd){
         testData[x] = rep(mean(dat[,x]), nrow(testData))
      }
      yHat = predict(rf,data=testData,type="response", num.threads=10)
      temp = prediction(yHat$predictions[,2], truePreds)
      roc = performance(temp, "tpr", "fpr")
      pr = performance(temp,"prec", "rec")
      auc = performance(temp,"auc")
      perf = list(roc = roc, pr = pr, auc = auc)
      return(perf)
   }, simplify=F)
   cat("\n")
   return(perfs)
}, simplify = F)
save(subBoth, file = "data/rdata/CrossTissue_subBoth.RData")
#####

# collect ROCs #####
load("data/rdata/CrossTissue_subTrain.RData")
load("data/rdata/CrossTissue_subTest.RData")
load("data/rdata/CrossTissue_subBoth.RData")
ROCs_subTrain = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
   t(sapply(tissues, function(predictingTissue){
      temp = subTrain[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
      return(c("trainedOn" = predictingTissue, 
               "testedOn" = tissue2predict, "AUC" = temp))
   }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subTrain$AUC = as.numeric(ROCs_subTrain$AUC)
ROCs_subTrain$rank = NA

ROCs_subTest = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
   t(sapply(tissues, function(predictingTissue){
      temp = subTest[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
      return(c("trainedOn" = predictingTissue, 
               "testedOn" = tissue2predict, "AUC" = temp))
   }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subTest$AUC = as.numeric(ROCs_subTest$AUC)
ROCs_subTest$rank = NA

ROCS_subBoth = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
   t(sapply(tissues, function(predictingTissue){
      temp = subBoth[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
      return(c("trainedOn" = predictingTissue, 
               "testedOn" = tissue2predict, "AUC" = temp))
   }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCS_subBoth$AUC = as.numeric(ROCS_subBoth$AUC)
ROCS_subBoth$rank = NA

# Add rowrank to each dataset
for(tissue in tissues){
  ROCs_subTrain[ROCs_subTrain$testedOn == tissue,]$rank[order(-ROCs_subTrain[ROCs_subTrain$testedOn == tissue,]$AUC)] = 1:length(tissues)
  ROCs_subTest[ROCs_subTest$testedOn == tissue,]$rank[order(-ROCs_subTest[ROCs_subTest$testedOn == tissue,]$AUC)] = 1:length(tissues)
  ROCS_subBoth[ROCS_subBoth$testedOn == tissue,]$rank[order(-ROCS_subBoth[ROCS_subBoth$testedOn == tissue,]$AUC)] = 1:length(tissues)
}
#####

# visualize ##### 
ggplot(ROCs_subTrain, aes(trainedOn, testedOn, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
   geom_text(aes(label = rank), color = "white", size = 2.5)+
   labs(y = "Tested on ", x = "Trained on") +  
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)
ggsave(filename="fig/CrossTissue_subTrain.png", 
       width = 12, height = 12, units = "cm")
ggplot(ROCs_subTest, aes(trainedOn, testedOn, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
   geom_text(aes(label = rank), color = "white", size = 2.5)+
   labs(y = "Tested on ", x = "Trained on") +  
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)
ggsave(filename="fig/CrossTissue_subTest.png", 
       width = 12, height = 12, units = "cm")
ggplot(ROCS_subBoth, aes(trainedOn, testedOn, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
   geom_text(aes(label = rank), color = "white", size = 2.5)+
   labs(y = "Tested on ", x = "Trained on") +  
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)
ggsave(filename="fig/CrossTissue_subBoth.png", 
       width = 12, height = 12, units = "cm")
#####
