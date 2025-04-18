.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ggplot2)
library(ROCR)
library(ranger) 
library(plotrix)
source("lib/general_function.R")
source("scripts/05_analysis/00_NamesAndColors.R")
dir.create("data/Modeling/WholeGenomeData/CrossTissue", showWarnings = F)
plotEnding = "_WGS_20241209"
tissues = c("brain", "breast", "esophagus", "kidney", "liver",
            "ovary", "prostate", "skin")
tissueCols = c(tissueCols[tissues], 
               "generalModel" = "darkgreen")
t2T = c(t2T[tissues],"generalModel" = "All-tissue general model")
nThread = 14

# Mutation overview #####
print("mutation overview")
# prepare counts of mutation types
bases = c("A", "C", "G", "T")
trimClass = cbind(paste0(paste0(rep(rep(bases, each = 4), 6),"[",
                                c(rep("C", 48), rep("T", 48)), ">",
                                c(rep(c("A", "G", "T"), each =16),
                                  rep(c("A", "C", "G"), each =16)), "]",
                                rep(bases, 24))),
                  paste0(paste0(rep(rep(rev(bases), each = 4),  6),"[",
                                c(rep("G", 48), rep("A", 48)), ">",
                                c(rep(c("T", "C", "A"), each =16),
                                  rep(c("T", "G", "C"), each =16)), "]",
                                rep(rev(bases), 24))))
rownames(trimClass) = paste(trimClass[,1], trimClass[,2], sep = " and ")
trimTypes = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
              tissue, "_mapped.RData"))
  muts = data$muts[data$muts$mutated == 1,]
  temp = paste0(substr(muts$context,2,2),"[",muts$ref,">", muts$alt,"]",
                substr(muts$context,4,4))
  apply(trimClass,1, function(x){
    sum(temp %in% x)
  })
})
trimTypes = cbind(trimTypes, "generalModel" = rowSums(trimTypes))
# plot
png(paste0("fig/MutationOverview_generalModel", plotEnding, ".png"),
    width = 2900, height = 1500, res=200, pointsize=18)
layout(rbind(1:9))
par(mar = c(0,1,2,0), oma = c(3.3,5,0,1))
temp =apply(trimTypes, 2,function(x){x/sum(x)})
dumpVar = sapply(1:ncol(temp), function(i){
  barplot(rev(temp[,i]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
          names.arg="",cex.names=0.5, cex.axis = 1.1,main  = c(t2T[tissues],"Combined")[i],
          col = rep(rev(basesCol), each=16), xlim = c(0,0.16), xaxp = c(0, 0.1, 1))
  if(i==1){
    classes = trimClass[rownames(temp),1]
    names = rev(paste0(substr(classes,1,1), substr(classes,3,3), substr(classes,7,7)))
    names2 = rev(paste0(" ", substr(classes,3,3), " "))
    text(x = 0, y = 0.99:96-0.5,  labels=names,family = "mono",
         xpd = NA, adj=1.1,  cex = 0.5)
    text(x = 0, y = 0.99:96-0.5,  labels=names2,family = "mono",
         xpd = NA, adj=1.1, col = rep(rev(basesCol), each=16), font = 2, cex = 0.5)
    text(x = -0.055, y = 1:6*16-8,
         labels = c("T>G", "T>C", "T>A" ,"C>T", "C>G", "C>A"),
         xpd = NA, adj=1, col = rev(basesCol), font = 2)
    segments(y0=1:6*16-16, y1 = 1:6*16, x0=-0.045, lwd = 4,
             col = rev(basesCol), xpd = NA, lend=1, cex = 1.1)
    title(ylab = "Mutation type", mgp = c(5,1,0), xpd = NA, cex.lab = 1.1)
  }
  abline(h=seq(16,80,length.out = 5), lty = 2, lwd = 1.5)
  abline(h=0, lty = 1, lwd = 2)
  #abline(h=96, lty = 1, lwd = 1.8)
})
axis(3,at =1:length(tissueCols)*2-0.5, labels= t2T[names(tissueCols)], tick=F,
     mgp = c(2,0.5,0), cex.axis  = 1.1)
title(xlab = "Proportion of variants", outer = T, mgp = c(2,1,0), cex.lab = 1.1)
dev.off()
rm(trimClass, bases, trimTypes, temp, dumpVar)
#####

# # compute performances #####
# load("data/Modeling/WholeGenomeData/RF/TissueCombination_predictions.RData") # predictions
# PerfPerChr  = lapply(predictions, function(x){ #iterate through chromosomes
#   perf = prediction(x$pred, x$label)
#   roc = performance(perf, "tpr", "fpr")
#   pr = performance(perf,"prec", "rec")
#   auc = performance(perf,"auc")@y.values[[1]]
#   return(list(roc = roc, pr = pr, auc = auc))
# })
# predConcat = do.call(rbind,predictions)
# temp = prediction(predConcat$pred, predConcat$label)
# PerfConcat = list(roc = performance(temp, "tpr", "fpr"),
#                   pr = performance(temp,"prec", "rec"),
#                   auc = performance(temp,"auc")@y.values[[1]])
# #####

# # ROC and PR of CWCV ####
# print("ROC and PR")
# 
# # plot
# png(paste0("fig/modelEvaluationWGS/ROC_PR_CWCV_generalModel", plotEnding, ".png"),
#     width = 1800, height = 1200, res=200, pointsize=12)
# layout(cbind(1,2))
# # ROC
# plot(NA, xlim = c(0,1), ylim = c(0,1),
#      mgp = c(2,0.7,0), las = 1,
#      xlab ="FPR", ylab = "TPR")
# abline(0,1, col = "grey", lty = 2)
# dump = sapply(names(chrCols), function(cr){
#   lines(PerfPerChr[[cr]]$roc@x.values[[1]],
#         PerfPerChr[[cr]]$roc@y.values[[1]],
#         col = rgb(0,0,0,0.3))
# })
# lines(PerfConcat$roc@x.values[[1]],
#       PerfConcat$roc@y.values[[1]],
#       col = tissueCols["generalModel"], lwd = 3)
# legend("bottomright", lty = 1, lwd = c(1,3), bty = "n",
#        col = c(rgb(0,0,0,0.5), tissueCols["generalModel"]),
#        legend=c("CWCV", "Total"))
# #  PR
# plot(NA, xlim = c(0,1), ylim = c(0,1),
#      mgp = c(2,0.8,0), las = 1,xlab = "Recall", ylab = "Precision")
# dump = sapply(names(chrCols), function(cr){
#   lines(PerfPerChr[[cr]]$pr@x.values[[1]],
#         PerfPerChr[[cr]]$pr@y.values[[1]],
#         col = rgb(0,0,0,0.3))
# })
# lines(PerfConcat$pr@x.values[[1]],
#       PerfConcat$pr@y.values[[1]],
#       col = tissueCols["generalModel"], lwd = 3)
# dev.off()
# #####


# # comparison of tissue performance #####
# print("comparison of tissues")
# load("data/Modeling/WholeGenomeData/RF/ROC_PR_RF_concat.RData")
# 
# png(paste0("fig/modelEvaluationWGS/compareRFperformance_generalModel",
#            plotEnding, ".png"),
#     width=2400, height=800, pointsize=35)
# par(mfrow = c(1,3), mar = c(5,5,0.1,0.1), mgp = c(2.5,1,0))
# # roc
# plot(NA, xlim = c(0,1), ylim = c(0,1),
#      xlab = "FPR",
#      ylab = "TPR", las = 1)
# abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
# plotDump = sapply(tissues, function(tissue){
#   lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
#         ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]],
#         col = tissueCols[tissue], lwd = 2)
# })
# lines(PerfConcat$roc@x.values[[1]],
#       PerfConcat$roc@y.values[[1]],
#       col = tissueCols["generalModel"], lwd = 2)
# # pr
# plot(NA, xlim = c(0,1), ylim = c(0,1),
#      xlab = "Recall",
#      ylab = "Precision", las = 1)
# plotDump = sapply(tissues, function(tissue){
#   lines(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
#         ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]],
#         col = tissueCols[tissue], lwd = 2)
# })
# lines(PerfConcat$pr@x.values[[1]],
#       PerfConcat$pr@y.values[[1]],
#       col = tissueCols["generalModel"], lwd = 2)
# legend("bottomright", col=tissueCols, legend=t2T[names(tissueCols)], lwd = 2)
# # AUCs
# AUCs = c(sapply(tissues, function(tissue){
#   ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
# }),   PerfConcat$auc)
# par(mar = c(5,10,0.1,0.1))
# temp = barplot(rev(AUCs), las = 1, xlab = "AUC", col = rev(tissueCols), names.arg = rev(t2T),
#                horiz = T,)
# # text(y = temp, x = -0.02, labels = t2T,srt = 45, c(1,0),
# #      xpd = T)
# dev.off()
# 
# #####



# # comparison of predictor importance #####
# print("comparison of predictor importances")
# rf_imps = do.call(rbind,sapply(names(tissueCols), function(tissue){
#   print(tissue)
#   if(tissue=="generalModel"){
#     load(paste0("data/Modeling/WholeGenomeData/RF/TissueCombination_", 
#                 "finalModel.RData"))
#   } else{
#     load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_",
#                 "finalModel.RData")) # rf, importance
#   }
#   res = data.frame("tissue" = tissue,
#                    "predictor" = names(p2P),
#                    "gini" = importance[names(p2P)],
#                    "gini_scaled" = scale(importance[names(p2P)], center = F))
#   return(res)
# }, simplify=F))
# rf_imps$predictor = factor(rf_imps$predictor, levels = names(p2P))
# rf_imps$tissue = factor(rf_imps$tissue, levels = names(t2T))
# save(rf_imps, file = "temp/rf_imps.RData")
# 
# ggplot(rf_imps, aes(x = tissue, y = predictor, fill = gini_scaled)) +
#   geom_raster() +
#   scale_fill_gradient(low="grey90", high="red",na.value="grey") +
#   labs(y = "Predictor", x = "Tissue") +
#   labs(fill = "Scaled Gini\nImportance")+
#   scale_x_discrete(labels=t2T) +
#   scale_y_discrete(labels = p2P)+
#   theme(axis.text.y=element_text(size=5.5),
#         axis.text.x=element_text(size=8 , angle = 45,vjust = 1, hjust=1),
#         axis.title=element_text(size=14,face="bold"))
# ggsave(paste0("fig/modelEvaluationWGS/RFginiScaled_generalModel", plotEnding, ".png"),
#        height=8, width=6)
# rm(rf_imps); gc()
# #####


# apply to other tissues #####
print("apply to other tissues")
crossTissuePerformance =  sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_", 
              tissue, "_mapped_processed.RData")) # dat, datchroms
  load("data/Modeling/WholeGenomeData/RF/TissueCombination_finalModel.RData")
  yHat = predict(rf,data=dat,type="response", num.threads=nThread)
  # compute performance
  temp = prediction(predictions = yHat$predictions[,2], labels = as.integer(dat$mutated)-1)
  auc = performance(temp,"auc")
  save(auc@y.values[[1]], file = paste0("data/Modeling/WholeExomeData/contextModel/generalModel_on", tissue, ".RData"))
  return(auc@y.values[[1]])
})
save(crossTissuePerformance, 
     file = "data/Modeling/WholeGenomeData/CrossTissue/generalModel_crossTissuePerformance.RData")
######


#visualize cross-tissue performance #####
load("data/Modeling/WholeGenomeData/CrossTissue/generalModel_crossTissuePerformance.RData")
load("data/Modeling/WholeGenomeData/RF/ROC_PR_RF_concat.RData")
tissueSpecificPerformance = sapply(tissues, function(tissue){
  ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
})
AUCs = rbind(tissueSpecificPerformance,crossTissuePerformance)
AUCs = AUCs[,ncol(AUCs):1]
png(paste0("fig/modelEvaluationWGS/generalModelCrossTissuePerformance", plotEnding, ".png"),
    width=400, height=800, pointsize=12)
par(mar = c(4,6,6,0.1))
temp = barplot(AUCs-0.5, beside = T,  las = 1, xlab = "AUC", horiz = T,mgp = c(2,1,0),
               names.arg = t2T[colnames(AUCs)], col = rep(tissueCols[colnames(AUCs)], each = 2),xaxt = "n",
               density = c(-1, 10), angle = c(-1,45), legend.text = F, bg = "grey")
axis(1, at = axTicks(1), labels = c(0,axTicks(1)[-1]+0.5))
axis.break(axis = 1, breakpos = 0.01, style = "gap")
legend(y = max(temp)*1.15, x =0, inset = 0.01,xpd = NA,
       legend = c("Tissue-specific model","All-tissue general model"),
       density = c(-1, 30), angle = c(-1,45),fill = "black")
dev.off()
#####


