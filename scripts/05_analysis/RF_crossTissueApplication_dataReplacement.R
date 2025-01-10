source("scripts/05_analysis/00_NamesAndColors.R")
library(ranger)
library(ROCR)
library(ggplot2)
cr= "chr1"
load("data/rdata/dataInfos.RData")
mindatsize = min(sapply(dataInfos, function(x)x$nMuts))

# variant 1: e.g. test lung model on lung mutations with kidney predictors #####
print("variant 1")
CrossTissuePreds = sapply(tissues, function(muttissue){
   cat("muts and model from ", muttissue, "\n")
   # load RF
   load(paste0("data/rdata/RFmodel/", muttissue, "_",
               cr, ".RData"))
   predictors = names(rf$variable.importance)

   cat("predictors from ")
   perfs = sapply(tissues, function(predtissue){
      cat(predtissue, " ")
      if(predtissue == muttissue){
         load(paste0("data/procData/traindata/traindata_processed_",
                     muttissue, ".RData"))
      }else{
         load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",
                     muttissue, "With",predtissue,"Preds_processed", ".RData"))
      }


      testData = dat[datchroms == cr,]
      truePreds = testData$mutated
      # add dummy variables for predictors present only in predtissue
      colnames(testData)[colnames(testData) == "aPhased_repeats"] = "aPhased_repeates"
      testData = testData[,colnames(testData) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testData)]
      load(paste0("data/procData/traindata/traindata_processed_",
                  muttissue, ".RData"))
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
save(CrossTissuePreds, file = "data/rdata/CrossTissuePreds.RData")
#####


# variant 2: e.g. test lung model on kidney mutations with lung preds #####
print("variant 2")
CrossTissueMuts = sapply(tissues, function(muttissue){
   cat("preds and model from ", muttissue, "\n")
   # load RF
   load(paste0("data/rdata/RFmodel/", muttissue, "_", 
               cr, ".RData"))
   predictors = names(rf$variable.importance)
   
   cat("muts from ")
   perfs = sapply(tissues, function(predtissue){
      cat(predtissue, " ")
      if(predtissue == muttissue){
         load(paste0("data/procData/traindata/traindata_processed_",
                     predtissue, ".RData"))
      }else{
         load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",
                     predtissue, "With",muttissue,"Preds_processed", ".RData"))
      }
      testData = dat[datchroms == cr,]
      truePreds = testData$mutated
      # add dummy variables for predictors present only in predtissue
      colnames(testData)[colnames(testData) == "aPhased_repeats"] = "aPhased_repeates"
      testData = testData[,colnames(testData) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testData)]
      load(paste0("data/procData/traindata/traindata_processed_",
                  muttissue, ".RData"))
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
save(CrossTissueMuts, file = "data/rdata/RFCrossTissueMuts.RData")
#####


# visualize results #####
# variant 1
load("data/rdata/CrossTissuePreds.RData")
RocsPreds = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
   t(sapply(tissues, function(predtissue){
      temp = CrossTissuePreds[[muttissue]][[predtissue]]$auc@y.values[[1]]
      return(c("predictorSource" = predtissue, 
               "MutationAndModelSource" = muttissue, "AUC" = temp))
   }, USE.NAMES= F))
})), stringsAsFactors=F)
RocsPreds$AUC = as.numeric(RocsPreds$AUC)
RocsPreds$rank <- NA
for(tissue in tissues){
  RocsPreds[RocsPreds$predictorSource == tissue,]$rank[order(-RocsPreds[RocsPreds$predictorSource == tissue,]$AUC)] = 1:length(tissues)
}
ggplot(RocsPreds, aes(y = predictorSource, x = MutationAndModelSource, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
   labs(y = "Predictor Source", x = "Mutation & Model Source") + 
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)+
   geom_text(aes(label = rank), color = "white", size = 2.5)
ggsave(filename="fig/CrossTissuePreds.png", 
       width = 12, height = 12, units = "cm")

# it looks like we have slightly higher performance than for the
# original cross-tissue application. Is that the case?
load("data/rdata/RFperformanceCrossTissue.RData")
ROCs = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
   t(sapply(tissues, function(predictingTissue){
      temp = RFperformanceCrossTissue[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
      return(c("trainedOn" = predictingTissue, 
               "testedOn" = tissue2predict, "AUC" = temp))
   }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs$AUC = as.numeric(ROCs$AUC)
all.equal(ROCs$trainedOn, RocsPreds$predictorSource)
all.equal(ROCs$testedOn, RocsPreds$MutationAndModelSource)

table(ROCs$AUC <= RocsPreds$AUC)
png("fig/crossTissue_origVsRocsPreds.png")
plot(ROCs$AUC, RocsPreds$AUC, xlim = c(0.47,0.7), ylim = c(0.47,0.7), 
     xlab = "original AUCs", ylab = "AUCs from exchanged predictors",
     col = (RocsPreds$AUC > ROCs$AUC)+3)
abline(0,1)
text(x = 0.55, y = 0.68, paste0("n = ",sum(RocsPreds$AUC > ROCs$AUC)), col = 4, font = 2)
text(x = 0.68, y = 0.55, paste0("n = ",sum(RocsPreds$AUC <= ROCs$AUC)), col = 3, font = 2)

dev.off()


# variant 2
load("data/rdata/RFCrossTissueMuts.RData")
RocsMuts = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
   t(sapply(tissues, function(predtissue){
      temp = CrossTissueMuts[[muttissue]][[predtissue]]$auc@y.values[[1]]
      return(c("MutationSource" = predtissue, 
               "PredictorAndModelSource" = muttissue, "AUC" = temp))
   }, USE.NAMES= F))
})), stringsAsFactors=F)
RocsMuts$AUC = as.numeric(RocsMuts$AUC)
RocsMuts$rank <- NA
for(tissue in tissues){
  RocsMuts[RocsMuts$MutationSource == tissue,]$rank[order(-RocsMuts[RocsMuts$MutationSource == tissue,]$AUC)] = 1:length(tissues)
}
ggplot(RocsMuts, aes(y = MutationSource, x = PredictorAndModelSource, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
   labs(y = "Mutation Source", x = "Predictor & Model Source") + 
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)+
   geom_text(aes(label = rank), color = "white", size = 2.5)
ggsave(filename="fig/CrossTissueMuts.png", 
       width = 12, height = 12, units = "cm")



