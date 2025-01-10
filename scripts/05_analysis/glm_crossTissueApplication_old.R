library(ROCR)
library(ggplot2)
source("scripts/05_analysis/00_NamesAndColors.R")
cr= "chr1"
load("data/rdata/dataInfos.RData")
mindatsize = min(sapply(dataInfos, function(x)x$nMuts))

# apply models to data from other tissues and get performance #####
GLMperformanceCrossTissue = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue2predict, ".RData"))
  tempData = dat[datchroms == cr,]
  truePreds = tempData$mutated
  cat("model: ")
  perfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    load(paste0("data/rdata/GLMmodel/", predictingTissue, "_", 
                cr, "_sig.RData"))
    load(paste0("data/procData/traindata/traindata_processed_",
                predictingTissue, ".RData"))
    
    # Just transfer on the sig. predictors
    load(file = paste0("data/rdata/GLMmodel/", predictingTissue,"_pvals.RData"))
    predictors = names(rowMeans(pvals)[rowMeans(pvals) < 0.05])
    testData = tempData[,colnames(tempData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    for(x in toAdd){ # Add mean value of predictors missing for the tissue to predict on (from predicting tissue)
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    # Predict on new transfered data 
    yHat = predict(logR, newdata = testData, type = "response")
    temp = prediction(yHat, truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
PerfsSelf = sapply(tissues, function(tissue){
  GLMperformanceCrossTissue[[tissue]][[tissue]]
}, simplify = F)
save(GLMperformanceCrossTissue, PerfsSelf, file = "data/rdata/GLMperformanceCrossTissue.RData")
######

# visualize results ######
# load( "data/rdata/GLMperformanceCrossTissue.RData")
# plot ROC
png("fig/modelOnOtherTissues_ROC.png",
    width=1250, height=1000, pointsize=27)
par(mfrow = c(3,4), mar = c(2,2.5,2,1))
for(predictingTissue in tissues){
  plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "FPR", ylab = "TPR", 
       main = paste0(t2T[predictingTissue],
                     " model"), las = 1, mgp = c(3,1,0))
  abline(0,1)
  for(tissue2predict in tissues){
    temp = GLMperformanceCrossTissue[[tissue2predict]][[predictingTissue]]
    plot(temp$roc, add = T, col = tissueCols[tissue2predict], lwd = 2)
    if(tissue2predict != predictingTissue){
      temp2 = PerfsSelf[[tissue2predict]]
      plot(temp2$roc, add = T, col = tissueCols[tissue2predict], 
           lty = 2, lwd = 2)
    }
  }
}
plot.new()
legend("center",legend=t2T[names(tissueCols)], col = tissueCols,cex = 1.2,
       lty = 1, bty = "n", lwd = 2)
plot.new()
legend("center", bty = "n",cex = 1.2,
       legend = c(paste0('"foreign"', " model"), 
                  "tissue-specific\nmodel"), lty = 1:2, lwd = 2)
dev.off()


# plot PR
png("fig/modelOnOtherTissues_PR.png",
    width=1250, height=1000, pointsize=27)
par(mfrow = c(3,4), mar = c(2,2.5,2,1))
for(predictingTissue in tissues){
  plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Recall", ylab = "Precision", 
       main = paste0(t2T[predictingTissue],
                     " model"), las = 1, mgp = c(2.5,1,0))
  for(tissue2predict in tissues){
    temp = GLMperformanceCrossTissue[[tissue2predict]][[predictingTissue]]
    plot(temp$pr, add = T, col = tissueCols[tissue2predict], lwd = 2)
    if(tissue2predict != predictingTissue){
      temp2 = PerfsSelf[[tissue2predict]]
      plot(temp2$pr, add = T, col = tissueCols[tissue2predict],
           lty = 2, lwd = 2)
    }
  }
}
plot.new()
legend("center",legend=t2T[names(tissueCols)], col = tissueCols,cex = 1.2,
       lty = 1, bty = "n", lwd = 2)
plot.new()
legend("center", bty = "n", cex = 1.2,
       legend = c(paste0('"foreign"', " model"),
                  "tissue-specific\nmodel"), lty = 1:2, lwd = 2)
dev.off()


# image of ROCs
ROCs = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = GLMperformanceCrossTissue[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs$AUROC = as.numeric(ROCs$AUROC)
ROCs$rank <- NA
for(tissue in tissues){
  ROCs[ROCs$testedOn == tissue,]$rank[order(-ROCs[ROCs$testedOn == tissue,]$AUROC)] = 1:length(tissues)
}
###
upperL = 0.69
lowerL = 0.48
###

ggplot(ROCs, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
#ggsave(filename="fig/CrossTissue.png", 
#       width = 12, height = 12, units = "cm", bg = "white", dpi = 200)
ggsave(filename="fig/CrossTissue.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)
#####