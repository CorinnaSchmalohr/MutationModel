tissues = c("luad",  "skin", "colon", "ovary", "kidney", "prostate", "breast")
library(ranger)
library(ROCR)
library(ggplot2)
cr = "chr1"
RFperformanceCrossTissue = sapply(tissues, function(tissue2predict){
   cat("data: ", tissue2predict, "\n")
   load(paste0("data/rdata/", tissue2predict, "/completeData_withBinwise.RData"))
   data = data[removed$chr == "chr1",]
   truePreds = data$mutated
   data$context = NULL
   data$pentamer = NULL
   data$trimer = NULL
   data$septamer = NULL
   data$inexon = NULL
   cat("model: ")
   perfs = sapply(tissues, function(predictingTissue){
      cat(predictingTissue, " ")
      load(paste0("data/rdata/", predictingTissue, "/RFmodel/",
                  cr, "_withBinWise.RData"))
      predictors = names(rf$variable.importance)
      tempData = data[,colnames(data) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(tempData)]
      tempData[toAdd] = rep(0, nrow(tempData))
      yHat = predict(rf,data=tempData,type="response", num.threads=10)
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
save(RFperformanceCrossTissue, file = "data/rdata/RFperformanceCrossTissue.RData")
PerfsSelf = sapply(tissues, function(tissue){
   RFperformanceCrossTissue[[tissue]][[tissue]]
}, simplify = F)

tissueCols = setNames(rainbow(length(tissues)), tissues)
# plot ROC
for(predictingTissue in tissues){
   png(paste0("fig/performance/ROC_rf_",
              predictingTissue,"ModelOnOtherTissues.png"), 
       width=800, height=800, pointsize=20)
   plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "FPR", ylab = "TPR", 
        main = paste0("performance of ", predictingTissue,
                      " model on other tissues"), las = 1)
   abline(0,1)
   for(tissue2predict in tissues){
      temp = RFperformanceCrossTissue[[tissue2predict]][[predictingTissue]]
      plot(temp$roc, add = T, col = tissueCols[tissue2predict], lwd = 2)
      if(tissue2predict != predictingTissue){
         temp2 = PerfsSelf[[tissue2predict]]
         plot(temp2$roc, add = T, col = tissueCols[tissue2predict], 
              lty = 2, lwd = 2)
      }
   }
   legend(x = 0.7,y = 0.6,legend=names(tissueCols), col = tissueCols,
          lty = 1, bty = "n", lwd = 2)
   legend(x = 0.6, y = 0.2, bty = "n",
          legend = c(paste0(predictingTissue, " model"), 
                     "tissue-specific model"), lty = 1:2, lwd = 2)
   dev.off()
}

# plot PR
for(predictingTissue in tissues){
   png(paste0("fig/performance/PR_rf_",
              predictingTissue,"ModelOnOtherTissues.png"), 
       width=800, height=800, pointsize=20)
   plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Recall", ylab = "Precision", 
        main = paste0("performance of ", predictingTissue,
                      " model on other tissues"), las = 1)
   for(tissue2predict in tissues){
      temp = RFperformanceCrossTissue[[tissue2predict]][[predictingTissue]]
      plot(temp$pr, add = T, col = tissueCols[tissue2predict], lwd = 2)
      if(tissue2predict != predictingTissue){
         temp2 = PerfsSelf[[tissue2predict]]
         plot(temp2$pr, add = T, col = tissueCols[tissue2predict],
              lty = 2, lwd = 2)
      }
   }
   legend(x = 0.7,y = 0.45,legend=names(tissueCols), col = tissueCols,
          lty = 1, bty = "n", lwd = 2)
   legend(x = 0.55, y = 0.1, bty = "n",
          legend = c(paste0(predictingTissue, " model"), 
                     "tissue-specific model"), lty = 1:2, lwd = 2)
   dev.off()
}

# image of ROCs
ROCs = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
   t(sapply(tissues, function(predictingTissue){
      temp = RFperformanceCrossTissue[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
      return(c("trainedOn" = predictingTissue, 
               "testedOn" = tissue2predict, "AUC" = temp))
   }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs$AUC = as.numeric(ROCs$AUC)
ggplot(ROCs, aes(trainedOn, testedOn, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c() +
   theme_minimal()
ggsave(filename="fig/performance/AUCheatmap_modelCrossTissueApplication.png", 
       width = 12, height = 12, units = "cm")
