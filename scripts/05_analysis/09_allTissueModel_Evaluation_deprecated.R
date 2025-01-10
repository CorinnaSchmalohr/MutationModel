library(ROCR)
library(sinaplot)
library(corrplot)
source("scripts/05_analysis/00_NamesAndColors.R")
predictorOrder = read.table("scripts/04_trainmodels/predictorOrder_glm.txt")[,1]

## Load feature p-values, features coefficients and predictions
load(paste0("data/rdata/GLMmodel/allTissues_pvals_subsampled.RData"))
pvals_sub = pvals
load(paste0("data/rdata/GLMmodel/allTissues_pvals.RData"))

load(paste0("data/rdata/GLMmodel/allTissues_predictions_subsampled_sig.RData"))
predictions_sub = predictions
load(paste0("data/rdata/GLMmodel/allTissues_predictions_sig.RData"))

load(file = paste0("data/rdata/GLMmodel/allTissues_importances_subsampled_sig.RData"))
imp_sub = imp
load(file = paste0("data/rdata/GLMmodel/allTissues_importances_sig.RData"))


## Comparison to tissue specific predictions 
predPerTissue = sapply(tissues, function(tissue){
  # Load predictions from tissue specific model
  load(paste0("data/rdata/GLMmodel/", tissue,"_predictions_sig.RData"))
  predictions_org = predictions
  predictions_org = do.call(rbind,predictions_org)

  # Load predictions from reduced generalized model on tissue specific data 
  load(file = paste0("data/rdata/GLMmodel/allTissues_predictions_on_", tissue, "_subsampled.RData"))
  predictions_sub = predictions
  predictions_sub = do.call(rbind,predictions)
  
  # Load predictions from full generalized model on tissue specific data 
  load(paste0("data/rdata/GLMmodel/allTissues_predictions_on_", tissue, ".RData"))
  predictions = do.call(rbind,predictions)
  
  return(cbind(pred_org = predictions_org$pred, pred = predictions, pred_sub = predictions_sub$pred))
}, simplify=F)


## CVCW and genome average performance of generalized model #####
ROC_PR_glm_perChr_allTissues = lapply(predictions, function(x){ 
  perf = prediction(x$pred, x$label)
  roc = performance(perf, "tpr", "fpr")
  pr = performance(perf,"prec", "rec")
  auc = performance(perf,"auc")@y.values[[1]]
  return(list(roc = roc, pr = pr, auc = auc))
})
save(ROC_PR_glm_perChr_allTissues, file = "data/rdata/GLMmodel/ROC_PR_glm_perChr_allTissues.RData")

ROC_PR_glm_perChr_allTissues_subsampled = lapply(predictions_sub, function(x){
  perf = prediction(x$pred, x$label)
  roc = performance(perf, "tpr", "fpr")
  pr = performance(perf,"prec", "rec")
  auc = performance(perf,"auc")@y.values[[1]]
  return(list(roc = roc, pr = pr, auc = auc))
})
save(ROC_PR_glm_perChr_allTissues, file = "data/rdata/GLMmodel/ROC_PR_glm_perChr_allTissues_subsampled.RData")

# genome average training performance of generalized model 
predConcat = do.call(rbind,predictions)
predConcat_subsampled = do.call(rbind,predictions_sub)

ROC_concat = performance(prediction(predConcat$pred, 
                                    predConcat$label), "tpr", "fpr")
AUC_concat = performance(prediction(predConcat$pred, 
                                    predConcat$label), "auc")
PR_concat = performance(prediction(predConcat$pred, 
                                   predConcat$label), "prec", "rec")

ROC_concat_subsampled = performance(prediction(predConcat_subsampled$pred, 
                                    predConcat_subsampled$label), "tpr", "fpr")
AUC_concat_subsampled = performance(prediction(predConcat_subsampled$pred, 
                                    predConcat_subsampled$label), "auc")
PR_concat_subsampled = performance(prediction(predConcat_subsampled$pred, 
                                   predConcat_subsampled$label), "prec", "rec")
#####

# Performance measure of tissue specific predictions #####
AUC_perTissue = sapply(tissues, function(tissue){
  predConcat = predPerTissue[[tissue]]
  # get performances
  AUC = performance(prediction(predConcat$pred.pred, 
                               predConcat$pred.label), "auc")
  AUC_sub = performance(prediction(predConcat$pred_sub, 
                               predConcat$pred.label), "auc")
  AUC_org = performance(prediction(predConcat$pred_org, 
                                   predConcat$pred.label), "auc")
  out = rbind(auc = AUC@y.values[[1]], auc_sub = AUC_sub@y.values[[1]], auc_org = AUC_org@y.values[[1]])
  colnames(out) = tissue
  
  return(out)
}, simplify=F)
AUC_perTissue = do.call(cbind,AUC_perTissue)
#####



# Plots #####
png(paste0("fig/modelEvaluation/summary_allTissues.png"),
    width=1200, height=1000, pointsize = 28)
m = rbind(c(1,1,1,3,3,3,3), 
          c(1,1,1,3,3,3,3), 
          c(1,1,1,3,3,3,3), 
          c(2,2,2,3,3,3,3), 
          c(2,2,2,3,3,3,3),
          c(2,2,2,3,3,3,3))
par(oma = c(0.1,0.1,1.3,0.1), mar = c(4,4,0.5,0))
layout(m)
# ROC over chr
plot(NA, xlim = c(0,1), ylim = c(0,1),
     mgp = c(2.7,1,0), yaxt="n", 
     xaxt="n", xlab = "", ylab = "")
abline(0,1, col = "grey", lty = 2, lwd = 2)
axis(1,cex.axis=1.3, las = 1)
axis(2,cex.axis=1.3, las = 1)
mtext("False positive rate", side = 1, cex= 0.9, line = 2.5)
mtext("True positive rate", side = 2, cex= 0.9, line = 3)
dump = sapply(seq_len(length(chrCols)), function(i){
  plot(ROC_PR_glm_perChr_allTissues[[i]]$roc, col = "grey70", add = T)
})
plot(ROC_concat, col = "black", add = T, lwd = 2, main = "ROC")

# PR over chr
plot(NA, xlim = c(0,1), ylim = c(0,1),
     mgp = c(2.7,1,0), yaxt="n", 
     xaxt="n", xlab = "", ylab = "")
axis(1,cex.axis=1.3, las = 1)
axis(2,cex.axis=1.3, las = 1)
mtext("Recall", side = 1, cex= 0.9, line = 2.5)
mtext("Precision", side = 2, cex= 0.9, line = 3)
dump = sapply(seq_len(length(chrCols)), function(i){
  plot(ROC_PR_glm_perChr_allTissues[[i]]$pr, col = "grey70", add = T)
})
plot(PR_concat, col = "black", add = T, lwd = 2, main = "PR")

# GLM coefficients
par(mar = c(4,17,0,0.7))
bxplt=boxplot(t(imp), las = 1, drop =T,
              horizontal = T, cex.axis = 0.5, xaxt = "n", yaxt = "n")
axis(1,cex.axis=1.3)
axis(2, at=seq(1, length(rownames(imp))), labels= p2P[rownames(imp)], las = 2)
mtext("GLM coefficients", side = 1, cex= 0.9, line = 2.5)
segments(y0= 1:nrow(imp), y1 = 1:nrow(imp),
         x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
abline(v = 0, lwd = 2)
dev.off()

# Correlation plot 
#load("data/rdata/dataInfos_allTissues.RData")
#png(paste0("fig/modelEvaluation/corrplot_allTissues.png"),
#    width=1550, height=1300, pointsize = 25)
#cors = dataInfos$cors
#rownames(cors) = p2P[rownames(cors)]
#colnames(cors) = p2P[colnames(cors)]
#corrplot(corr = cors, tl.cex=0.8, tl.col="black")
#dev.off()


# Comparison 
png(paste0("fig/modelEvaluation/AUC_comp_AllTissues.png"),
    width=1200, height=600, pointsize = 20)
par(oma = c(0.1,0.1,1.3,0.1), mar = c(4,4,0.5,0))
barplot(height = AUC_perTissue, beside = TRUE, ylab = "AUROC", las = 1,
        names.arg = t2T, col = c("grey70", "grey40", "gray10"), 
        legend.text = c("Generalized model", "Generalized model on subsampled data", "Tissue-specific model"),
        args.legend = list(x = "bottomright"))
dev.off()

