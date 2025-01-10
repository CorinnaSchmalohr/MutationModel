library(ggplot2)
library(data.table) 
library(ROCR)
library(ranger) 
library(zoo)
library(berryFunctions) # for smallPlot
library(corrplot)
library(sinaplot)
library(plotrix) # for point labels                            
library(scales)
library(RColorBrewer)
source("lib/general_function.R")
source("scripts/05_analysis/00_NamesAndColors.R")

dir.create("fig/modelEvaluation", showWarnings=F)

predictorOrder = read.table("scripts/04_trainmodels/predictorOrder_glm.txt")[,1]
rarePredictors = c("ZBTB33_100bp", "YY1_100bp", "TAF1_100bp", "SP1_100bp",
                   "RXRA_100bp", "REST_100bp", "RAD21_100bp",
                   "NR2F2_100bp", "MAX_100bp", "JUND_100bp", "HNF4G_100bp",
                   "HNF4A_100bp", "GABPA_100bp", "FOXA2_100bp", 
                   "FOXA1_100bp", "EGR1_100bp", "ATF3_100bp")

## Data gathering ##
# load CWCV GLM predictions for each tissue #####
predPerTissueGLM = sapply(tissues, function(tissue){
  load(file = paste0("data/rdata/GLMmodel/", tissue, "_predictions_sig.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueGLM, file = "data/rdata/GLMmodel/predPerTissueGLM.RData")
#####

# get ROC, PR, and AUROC for each chr from CWCV glm #####
ROC_PR_glm_perChr = sapply(tissues, function(tissue){
  pred = predPerTissueGLM[[tissue]]
  # get performances
  res = lapply(pred, function(x){ #iterate through chromosomes
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  })
  return(res)
}, simplify=F)
save(ROC_PR_glm_perChr, file = "data/rdata/GLMmodel/ROC_PR_glm_perChr.RData")
#####

#####
# get ROC, PR, and AUROC concatenated from CWCV glm and the windowed preformances #####
ROC_PR_glm_concat = sapply(tissues, function(tissue){
  predConcat = do.call(rbind,predPerTissueGLM[[tissue]])
  ROC = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC, pr = PR, auc = AUC))
}, simplify=F)
save(ROC_PR_glm_concat, file = "data/rdata/GLMmodel/ROC_PR_glm_concat.RData")
#####

# get ROC, PR, and AUROC concatenated from windowed CWCV glm predictions #####
windows = c("10kb" = 10000, "1kb" = 1000, "100bp" = 100)

ROC_PR_glm_concat_window = lapply(tissues, function(tissue){
  predConcat = do.call(rbind,predPerTissueGLM[[tissue]])
  
  temp = lapply(windows, function(w){
    # partial = TRUE, then the subset of indexes that are in range are passed to FUN
    windowed_preformance = rollapply(predConcat$pred, width = w, by = 1, FUN = mean, align = "center", partial = TRUE)
    
    ROC = performance(prediction(windowed_preformance, predConcat$label), "tpr", "fpr")
    AUC = performance(prediction(windowed_preformance, predConcat$label), "auc")
    PR = performance(prediction(windowed_preformance, predConcat$label), "prec", "rec")
    return(list(roc = ROC, pr = PR, auc = AUC))
  })
  return(temp)
})
names(ROC_PR_glm_concat_window) = tissues
save(ROC_PR_glm_concat_window, file = "data/rdata/GLMmodel/ROC_PR_glm_concat_window.RData")
#####

# get predictor coefficients and pvals from final glm #####
glm_imps = lapply(tissues, function(tissue){
  load(file = paste0("data/rdata/GLMmodel/", tissue, "_sig.RData"))
  imp_sig = logR$coefficients
  load(file = paste0("data/rdata/GLMmodel/", tissue, ".RData"))
  imp = logR$coefficients
  
  temp = sapply(predictorOrder, function(x){
    if(x %in% names(imp_sig)){
      return(c(imp[x], imp_sig[x]))
    } else if(x %in% names(imp)[!names(imp) %in% names(imp_sig)]){
      return(c(imp[x], 0))
    } else {
      return(c(NA, NA))
    }
  })
  rownames(temp) = c("coefficients", "coefficients_sig")
  return(t(as.data.frame(temp)))
})
names(glm_imps) = tissues

glm_pvals = lapply(tissues, function(tissue){
  load(file = paste0("data/rdata/GLMmodel/", tissue, "_sig.RData"))
  pvals_sig = coef(summary(logR))[,4][-1]
  load(file = paste0("data/rdata/GLMmodel/", tissue, ".RData"))
  pvals = coef(summary(logR))[,4][-1]
  
  temp = sapply(predictorOrder, function(x){
    if(x %in% names(pvals_sig)){
      return(c(pvals[x], pvals_sig[x]))
    } else if(x %in% names(pvals)[!names(pvals) %in% names(pvals_sig)]){
      return(c(pvals[x], NA))
    } else {
      return(c(NA, NA))
    }
  })
  rownames(temp) = c("coefficients", "coefficients_sig")
  temp[temp == 0] = 2e-16
  return(t(as.data.frame(temp)))
})
names(glm_pvals) = tissues
save(glm_imps, glm_pvals, file = "data/rdata/GLMmodel/GLM_impsAndPvals.RData")
#####

#####

# load data #####
load("data/rdata/dataInfos.RData")
load("data/rdata/GLMmodel/predPerTissueGLM.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_perChr.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat_window.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals.RData")
load("data/rdata/nPerTissue.RData")
#####


# collection of plots for every tissue #####
plotDump = sapply(tissues, function(tissue){
  #pdf(paste0("fig/modelEvaluation/perfGLMsummary_", tissue, ".pdf"),
  #    width=12, height=10, pointsize = 20)
  png(paste0("fig/modelEvaluation/perfGLMsummary_", tissue, ".png"),
      width=1670, height=1300, pointsize = 38)
  m = rbind(c(1,1,2,2,5,5,5,5), 
            c(1,1,2,2,5,5,5,5), 
            c(3,3,4,4,5,5,5,5), 
            c(3,3,4,4,5,5,5,5), 
            c(6,6,6,6,5,5,5,5),
            c(6,6,6,6,5,5,5,5))
  par(oma = c(0.1,0.1,0.1,0.1), mar = c(4,4,0.5,0))
  layout(m)
  
  # 1. roc over Chrs 
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2.7,1,0), las = 1, 
       xlab = "False positive rate", ylab = "True positive rate")
  abline(0,1, col = "grey", lty = 2, lwd = 2)
  dump = sapply(seq_len(length(chrCols)), function(i){
    plot(ROC_PR_glm_perChr[[tissue]][[i]]$roc, 
         col = "grey70", add = T)
  })
  plot(ROC_PR_glm_concat[[tissue]]$roc, 
       col = "black", add = T, lwd = 2, main = "ROC")
  
  # 2. pr over Chrs
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2.7,1,0), las = 1,
       xlab = "Recall", ylab = "Precision")
  dump = sapply(seq_len(length(chrCols)), function(i){
    plot(ROC_PR_glm_perChr[[tissue]][[i]]$pr, 
         col = "grey70", add = T)
  })
  plot(ROC_PR_glm_concat[[tissue]]$pr, 
       col = "black", add = T, lwd = 2, main = "PR")
  
  # 3.  AUC vs mutrate
  AUCs = sapply(ROC_PR_glm_perChr[[tissue]], function(y){y$auc})
  
  # 3. AUC vs nPos 
  nMuts = dataInfos[[tissue]]$nMutsPerChr
  plot(nMuts, AUCs, las = 1, mgp = c(2.7,1,0),
       ylab = "AUROC", xlab = "*1000 positions",yaxt = "n")
  axis(side = 2, at = round(seq(min(AUCs), max(AUCs), length = 6), digits = 2), las = 1)#, labels = round(seq(min(AUCs), max(AUCs), length = 6), digits = 2))
  mod = lm(AUCs~nMuts)
  abline(mod, col = "grey", lty = 2, lwd = 2)
  p = summary(mod)$coefficients[2,4]
  r2 = summary(mod)$adj.r.squared
  rp = vector('expression',2)
  rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                     list(MYVALUE = format(r2,dig=3)))[2]
  rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                     list(MYOTHERVALUE = format(p, digits = 2)))[2]
  legend('bottomright', legend = rp, bty = 'n')
  
  # 4. TP vs TN preds violin plot
  par(mar = c(4,4,0.5,0))
  predsConcat = do.call(rbind,predPerTissueGLM[[tissue]])
  sinaplot(pred ~ label, predsConcat, las = 1, 
           xlab = "True value", ylab = "Predicted value", 
           col = c("light grey", "dark grey"), maxwidth = 0.9,mgp = c(2.7,1,0))
  boxplot(pred ~ label, predsConcat, add = T, yaxt = "n",xaxt = "n",
          col = "white", boxwex = c(0.5,0.5), outline = F)
  
  # 5. plot of all significant glm coefficients of the final model
  par(mar = c(4,17.5,0.5,0.7))
  imps = as.data.frame(glm_imps[[tissue]][,2])
  imps[imps==0] = NA
  imps = na.omit(imps)
  imps = imps[!names(imps) %in% rarePredictors, ,drop = FALSE]
  imps = cbind(imps, c(1:nrow(imps)))
  bxplt=boxplot(t(imps), las = 1, drop =T,horizontal = T, cex.axis = 0.5,xaxt = "n", yaxt = "n", plot = FALSE)
  plt = plot(imps, las = 1, cex.axis = 0.5, col = "white",
             xaxt = "n", yaxt = "n", xlab = " ", ylab = " ")
  axis(1)
  axis(2, at=seq(1, length(rownames(imps))), labels= p2P[rownames(imps)], las = 2)
  mtext("Final GLM coefficients", side = 1, cex= 2/3, line = 2.5)
  segments(y0= 1:nrow(imps), y1 = 1:nrow(imps),
           x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
  points(imps, pch = 19)
  abline(v = 0, lwd = 2)
  
  
  # 6. mutation type
  load(paste0("data/procData/traindata/traindata_", tissue, ".RData"))
  muts = data$muts
  muts = data$muts[ data$muts$mutated == 1,]
  muttypes = paste0(muts$ref, ">", muts$alt)
  muttypes = table(muttypes)
  mutDat = rbind(cbind(muttypes[4:6], -muttypes[9:7]), cbind(muttypes[10:12], -muttypes[3:1]))
  labels = paste(c(names(muttypes)[4:6], names(muttypes)[10:12]), 
                 c(names(muttypes)[9:7], names(muttypes)[3:1]), sep=" / ")
  par(mar = c(4,5.5,0.5,0))
  barplot(rev(mutDat[,1]), horiz=T, col = rev(basesCol),
          las = 1, xlim = c(-max(mutDat),max(mutDat)),names.arg=rev(labels),
          cex.lab=1,cex.names=1, xaxt = "n", xlab = "")
  barplot(rev(mutDat[,2]), add = T, horiz = T, col = rev(basesCol),xaxt = "n", yaxt = "n")
  mtext("Count", side = 1, cex= 2/3, line = 2.5)
  temp = axTicks(1)
  axis(1, at = temp, labels=abs(temp))
  dev.off()
  
  # separately: corrplot
  png(paste0("fig/modelEvaluation/corrplot_", tissue, ".png"),
      width=1600, height=1400, pointsize = 28)
  cors = dataInfos[[tissue]]$cors
  rownames(cors) = p2P[rownames(cors)]
  colnames(cors) = p2P[colnames(cors)]
  corrplot(corr = cors, tl.cex=0.8, tl.col="black", type = 'lower', tl.srt = 45)
  dev.off()
})
#####



# collection of plots for every tissue improved #####
png(paste0("fig/modelEvaluation/GLM_ROCsOverChrsAllTissues.png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.7,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  abline(0,1, col = "grey", lty = 2, lwd = 2)
  dump = sapply(names(chrCols), function(cr){
    plot(ROC_PR_glm_perChr[[tissue]][[cr]]$roc, 
         col = rgb(0,0,0,0.5), add = T)
  })
  lines(ROC_PR_glm_concat[[tissue]]$roc@x.values[[1]], 
        ROC_PR_glm_concat[[tissue]]$roc@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2.5)
  text(x=0.05,y=0.9,  labels=t2T[tissue],font= 2, adj= c(0,0), cex= 1.3)
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
})
title(ylab = "TPR", mgp = c(2,0.7,0), xpd = NA, outer=T, cex.lab = 1.3)
title(xlab = "FPR", mgp = c(1.8,0.7,0), xpd = NA, outer=T, cex.lab = 1.3)
legend("bottomright", lty = 1, lwd = c(1.5,3), bty = "n",
       col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"), cex = 1.2)
dev.off()

png(paste0("fig/modelEvaluation/GLM_PRsOverChrsAllTissues.png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.8,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  dump = sapply(names(chrCols), function(cr){
    plot(ROC_PR_glm_perChr[[tissue]][[cr]]$pr, 
         col = rgb(0,0,0,0.5), add = T)
  })
  plot(ROC_PR_glm_concat[[tissue]]$pr, 
       col = tissueCols[tissue], add = T, lwd = 2.5, main = "")
  text(x=0.5,y=0.9,  labels=t2T[tissue],font= 2, adj= c(0.5,0), cex= 1.3)
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
})
title(xlab = "Recall", mgp = c(1.8,0.7,0), xpd = NA, outer = T, cex.lab = 1.3)
title(ylab = "Precision", mgp = c(2,0.7,0), xpd = NA, outer = T, cex.lab = 1.3)
legend("bottomright", lty = 1, lwd = c(1,2.5), bty = "n",
       col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"), cex = 1.2)
dev.off()

# comparison of tissues performance #####
AUCs = sapply(tissues, function(tissue){
  ROC_PR_glm_concat[[tissue]]$auc@y.values[[1]]
})
AUCsPerChr = sapply(tissues, function(tissue){
  sapply(ROC_PR_glm_perChr[[tissue]], function(x){
    x$auc
  })
})
AUCsSD = apply(AUCsPerChr,2,sd)
nMuts = sapply(dataInfos, function(x){x$nMuts})

png("fig/modelEvaluation/compareGLMperformanceTissues.png", 
    width=2400, height=800, pointsize=35)
par(mfrow = c(1,3), mar = c(4,4,0.1,0.1), mgp = c(2.5,1,0))
# roc
plot(NA, xlim = c(0,1), ylim = c(0,1), 
     xlab = "False positive rate", 
     ylab = "True positive rate", las = 1)
abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
plotDump = sapply(tissues, function(tissue){
  lines(ROC_PR_glm_concat[[tissue]]$roc@x.values[[1]], 
        ROC_PR_glm_concat[[tissue]]$roc@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2)
})
legend("bottomright", col=tissueCols, legend=t2T[names(tissueCols)], lwd = 2)
# pr
plot(NA, xlim = c(0,1), ylim = c(0,1), 
     xlab = "Recall", 
     ylab = "Precision", las = 1)
plotDump = sapply(tissues, function(tissue){
  lines(ROC_PR_glm_concat[[tissue]]$pr@x.values[[1]], 
        ROC_PR_glm_concat[[tissue]]$pr@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2)
})
## AUCs vs nMuts ##
plot(nMuts, colMeans(AUCsPerChr), las = 1,
     bg = tissueCols, pch = 21, cex = 1.5, 
     ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
              max(colMeans(AUCsPerChr)+AUCsSD)), 
     xlim = c(min(nMuts), max(nMuts)+50000), 
     xlab = "n positions", ylab = "AUROC")
abline(lm(colMeans(AUCsPerChr) ~ nMuts), col = "grey", lty = 2)
arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nMuts,
       y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nMuts, code = 3,
       angle = 90, length = 0.1)
points(nMuts, colMeans(AUCsPerChr),
       bg = tissueCols, pch = 21, cex = 1.5)
points(x=nMuts, y = AUCs, pch = 4, bg = tissueCols, cex = 1.5)
text(x=nMuts, y = colMeans(AUCsPerChr), labels=t2T[tissues], pos=4)
dev.off()


# N muts overview plots 
pdf(paste0("fig/modelEvaluation/GLM_AUROCvsNmutsComparison.pdf"),
    width=17, height=6, pointsize = 26)
#png("fig/modelEvaluation/GLM_AUROCvsNmutsComparison_120822.png",
#    width=2400, height=800, pointsize=35)
par(mfrow = c(1,3), oma = c(3.5,3.8,0.2,0.8), mar = c(1,1,0,00))
# plot again: AUCs vs nMuts per tissue
plot(nMuts, colMeans(AUCsPerChr), las = 1,
     bg = tissueCols, pch = 21, cex = 1.5, 
     ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
              max(colMeans(AUCsPerChr)+AUCsSD)), 
     xlim = c(min(nMuts), max(nMuts)+65000))#, cex.axis = 1.3)
mod = lm(colMeans(AUCsPerChr) ~ nMuts)
abline(mod, col = "grey", lty = 2, lwd = 2)
p = summary(mod)$coefficients[2,4]
r2 = summary(mod)$adj.r.squared
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(p, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n')
arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nMuts,
       y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nMuts, code = 3,
       angle = 90, length = 0.1)
points(nMuts, colMeans(AUCsPerChr),
       bg = tissueCols, pch = 21, cex = 1.5)
text(x=nMuts, y = colMeans(AUCsPerChr), labels=t2T[tissues], offset=0.5, pos=4)
title(ylab = "AUROC", mgp = c(2.5,0.7,0), xpd = NA, outer = T)#, cex.lab = 1.6)
title(xlab = "Number of positions", mgp = c(2.7,0.7,0), xpd = NA, outer=F)#, cex.lab = 1.6)

# AUC vs nMuts per sample
medianNmuts =  sapply(nPerTissue, median)
quantNmuts = sapply(nPerTissue, quantile,probs=c(0.25, 0.75))
plot(medianNmuts, colMeans(AUCsPerChr), las = 1,
     bg = tissueCols, pch = 21, cex = 1.5, 
     xlim = c(min(quantNmuts), max(quantNmuts)),
     ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
              max(colMeans(AUCsPerChr)+AUCsSD)), 
     xlab = "", ylab = "",yaxt="n")#, cex.lab = 1.3, cex.axis = 1.3)
mod = lm(colMeans(AUCsPerChr) ~ medianNmuts)
abline(mod, col = "grey", lty = 2, lwd = 2)
p = summary(mod)$coefficients[2,4]
r2 = summary(mod)$adj.r.squared
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(p, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n')
arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = medianNmuts,
       y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = medianNmuts, code = 3,
       angle = 90, length = 0.1)
arrows(y0 = colMeans(AUCsPerChr), x0 = quantNmuts[2,],
       y1 = colMeans(AUCsPerChr), x1 = quantNmuts[1,],
       code = 3, angle = 90, length = 0.1)
points(medianNmuts, colMeans(AUCsPerChr),
       bg = tissueCols, pch = 21, cex = 1.5)
text(x=medianNmuts, y = colMeans(AUCsPerChr), labels=t2T[tissues], pos=4, offset=0.5)
title(xlab = "Number of mutations per sample", mgp = c(2.7,0.7,0), xpd = NA, outer=F)#, cex.lab = 1.6)

# AUC vs nSamples per tissue
nSamples = sapply(nPerTissue, length)
plot(nSamples, colMeans(AUCsPerChr), las = 1,
     bg = tissueCols, pch = 21, cex = 1.5, 
     xlim = c(min(nSamples), max(nSamples)+150),
     ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
              max(colMeans(AUCsPerChr)+AUCsSD)), 
     xlab = "", ylab = "",yaxt="n")#, cex.lab = 1.3, cex.axis = 1.3)
mod = lm(colMeans(AUCsPerChr) ~ nSamples)
abline(mod, col = "grey", lty = 2, lwd = 2)
p = summary(mod)$coefficients[2,4]
r2 = summary(mod)$adj.r.squared
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2,dig=3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(p, digits = 2)))[2]
legend('bottomright', legend = rp, bty = 'n')
arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nSamples,
       y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nSamples, code = 3,
       angle = 90, length = 0.1)
points(nSamples, colMeans(AUCsPerChr),
       bg = tissueCols, pch = 21, cex = 1.5)
text(x=nSamples, y = colMeans(AUCsPerChr), labels=t2T[tissues], pos=4, offset=0.5)
title(xlab = "Number of samples per tissue", mgp = c(2.7,0.7,0), xpd = NA, outer=F)#, cex.lab = 1.6, )
dev.off()
#####


# glm p-value heatmap #####
pvals_table = do.call(rbind,lapply(tissues, function(tissue){
  d = data.frame(pvalue = glm_pvals[[tissue]][,1], logpvalue = -log10(glm_pvals[[tissue]][,1]),tissue, direction = -log10(glm_pvals[[tissue]][,1]))
  d$predictor = names(glm_pvals[[tissue]][,1])
  for(f in d$predictor){
    f_dir = glm_imps[[tissue]][,1][f] # check effect direction
    if(f_dir < 0 && !is.na(f_dir)){
      d[f,]$direction = -d[f,]$direction # change sign of pvalue based on direction
    }
  }
  return(d)
}))
pvals_table$predictor = factor(pvals_table$predictor,levels=predictorOrder)
pvals_table["direction_cutoffs"] = pvals_table$direction
pvals_table$direction_cutoffs[ifelse(pvals_table$direction < -20, T, F)] = -20
pvals_table$direction_cutoffs[ifelse(pvals_table$direction > 20, T, F)] = 20

# removing TFs only present in liver
subMeanPvals = pvals_table[!pvals_table$predictor %in% c(rarePredictors),]

# GLM P-values
ggplot(pvals_table, aes(x = tissue, y = predictor, colour="")) + 
  geom_raster(aes(fill=direction_cutoffs)) +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(na.omit(pvals_table$direction_cutoffs)),0,max(na.omit(pvals_table$direction_cutoffs))),
    labels = c(round(min(na.omit(pvals_table$direction_cutoffs), digits = -1)), 0, round(max(na.omit(pvals_table$direction_cutoffs)), digits = -1)))+
  geom_text(aes(label = ifelse(pvalue < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Effect score")+ # labs(fill = "Log p-value &\neffect direction")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = p2P[as.character(pvals_table$predictor)])+
  scale_colour_manual(name = element_blank(), values=c("grey50", "grey50"))+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
ggsave("fig/GLMpvals_heatmap.pdf", height=8, width=7, bg = "white", dpi = 200)

ggplot(subMeanPvals, aes(x = tissue, y = predictor, colour="")) + 
  geom_raster(aes(fill=direction_cutoffs)) +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(na.omit(pvals_table$direction_cutoffs)),0,max(na.omit(pvals_table$direction_cutoffs))),
    labels = c(round(min(na.omit(pvals_table$direction_cutoffs), digits = -1)), 0, round(max(na.omit(pvals_table$direction_cutoffs)), digits = -1)))+
  geom_text(aes(label = ifelse(pvalue < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Effect score")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = p2P[as.character(subMeanPvals$predictor)])+
  scale_colour_manual(name = element_blank(), values=c("grey50", "grey50"))+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
#ggsave("fig/GLMpvals_heatmap_noTFs.png", height=8, width=7,, bg = "white", dpi = 200)
ggsave("fig/GLMpvals_heatmap_noTFs_V3.pdf", height=8, width=7, bg = "white", dpi = 200)


# glm feature coefficients #####
imps_table$predictor = factor(imps_table$predictor,levels=predictorOrder)
# removing TFs only present in liver
subMeanImps = imps_table[!imps_table$predictor %in% c(rarePredictors),]

# GLM feature importances 
ggplot(imps_table, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=importance)) +
  scale_fill_gradient2(low="blue", mid="grey90", high = "red", na.value="grey")+
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Feature\ncoefficients")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = p2P[as.character(imps_table$predictor)])+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
ggsave("fig/GLMimportances_heatmap.pdf", height=8, width=7, bg = "white", dpi = 200)


ggplot(subMeanImps, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=importance)) +
  scale_fill_gradient2(low="blue", mid="grey90", high = "red", na.value="grey")+
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Feature\ncoefficients")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = p2P[as.character(imps_table$predictor)])+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length= unit(0, "cm"),
        axis.line= element_line(linewidth = NA))
ggsave("fig/GLMimportances_heatmap_noTFs.pdf", height=8, width=7, bg = "white", dpi = 200)
#####

# Windowed prediction performance #####
window_cols = c("grey50", "grey40", "grey30")

# roc
png(paste0("fig/modelEvaluation/GLM_ROCs_WindowedPredictions.png"),
    width=1500, height=500, res=150)
par(oma = c(3,3.5,0.2,1), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.7,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  abline(0,1, col = "grey", lty = 2, lwd = 1)
  dump = sapply(1:length(windows), function(i){
    plot(ROC_PR_glm_concat_window[[tissue]][[names(windows[i])]]$roc, 
         col = window_cols[i], add = T, lwd = 1.5)
  })
  lines(ROC_PR_glm_concat[[tissue]]$roc@x.values[[1]], 
        ROC_PR_glm_concat[[tissue]]$roc@y.values[[1]], 
        col = tissueCols[tissue], lwd = 1.5)
  text(x=0.05,y=0.9,  labels=t2T[tissue],font= 2, adj= c(0,0), cex= 1.3)
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
})
title(ylab = "TPR", mgp = c(2,0.7,0), xpd = NA, outer=T, cex.lab = 1.3)
title(xlab = "FPR", mgp = c(1.8,0.7,0), xpd = NA, outer=T, cex.lab = 1.3)
legend("bottomright", lty = 1, lwd = 3, bty = "n",
       col = c(window_cols, "black"), legend=c(names(windows), "1bp"), cex = 1.2)
dev.off()

# pr
png(paste0("fig/modelEvaluation/GLM_PRs_WindowedPredictions.png"),
    width=1500, height=500, res=150)
par(oma = c(3,3.5,0.2,1), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.8,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  dump = sapply(1:length(windows), function(i){
    plot(ROC_PR_glm_concat_window[[tissue]][[names(windows[i])]]$pr, 
         col = window_cols[i], add = T, lwd = 1.5)
  })
  plot(ROC_PR_glm_concat[[tissue]]$pr, 
       col = tissueCols[tissue], add = T, lwd = 1.5, main = "")
  text(x=0.5,y=0.9,  labels=t2T[tissue],font= 2, adj= c(0.5,0), cex= 1.3)
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0), cex.axis=1.3)
  }
})
title(ylab = "TPR", mgp = c(2,0.7,0), xpd = NA, outer=T, cex.lab = 1.3)
title(xlab = "FPR", mgp = c(1.8,0.7,0), xpd = NA, outer=T, cex.lab = 1.3)
legend("bottomright", lty = 1, lwd = 3, bty = "n",
       col = c(window_cols, "black"), legend=c(names(windows), "1bp"), cex = 1.2)
dev.off()
#####
