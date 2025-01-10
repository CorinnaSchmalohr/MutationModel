library(ggplot2)
library(data.table) 
library(ROCR)
library(ranger) 
library(berryFunctions) # for smallPlot
library(corrplot)
library(sinaplot)
library(plotrix) # for point labels                            # Install & load scales
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
# load GLM predictions for each tissue #####
predPerTissueGLM = sapply(tissues, function(tissue){
  load(paste0("data/rdata/GLMmodel/", tissue,
              "_predictions.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueGLM, file = "data/rdata/GLMmodel/predPerTissueGLM.RData")


# Load predictions of significant model #####
predPerTissueGLMsig = sapply(tissues, function(tissue){
  load(file = paste0("data/rdata/GLMmodel/", tissue,
                     "_predictions_sig.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueGLMsig, file = "data/rdata/GLMmodel/predPerTissueGLM_sig.RData")
#####

# get ROC, PR, and AUROC for each chr from glm #####
ROC_PR_glm_perChr = sapply(names(predPerTissueGLM), function(tissue){
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

# and for the only sig. glm
ROC_PR_glm_perChr_sig = sapply(names(predPerTissueGLMsig), function(tissue){
  pred = predPerTissueGLMsig[[tissue]]
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
save(ROC_PR_glm_perChr_sig, file = "data/rdata/GLMmodel/ROC_PR_glm_perChr_sig.RData")

#####

#####
# get ROC, PR, and AUROC concatenated from glm #####
ROC_PR_glm_concat = sapply(names(predPerTissueGLM), function(tissue){
  predConcat = do.call(rbind,predPerTissueGLM[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_glm_concat, file = "data/rdata/GLMmodel/ROC_PR_glm_concat.RData")

# and for the sig. glm
ROC_PR_glm_concat_sig = sapply(names(predPerTissueGLMsig), function(tissue){
  predConcat = do.call(rbind,predPerTissueGLMsig[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_glm_concat_sig, file = "data/rdata/GLMmodel/ROC_PR_glm_concat_sig.RData")
#####

# get predictor importances and pvals from glm #####
glm_imps = sapply(tissues, function(tissue){
  load(paste0("data/rdata/GLMmodel/", tissue,
              "_importances.RData"))
  temp = t(sapply(predictorOrder, function(x){
    if(x %in% rownames(imp)){
      return(imp[x,])
    } else{
      return(rep(NA, ncol(imp)))
    }
  }))
  return(temp)
}, simplify=F)
glm_pvals = sapply(tissues, function(tissue){
  load(paste0("data/rdata/GLMmodel/", tissue,
              "_pvals.RData"))
  pvals = as.matrix(pvals)
  temp = t(sapply(predictorOrder, function(x){
    if(x %in% rownames(pvals)){
      return(pvals[x,])
    } else{
      return(rep(NA, ncol(pvals)))
    }
  }))
  temp[temp == 0] = 2e-16
  return(temp)
}, simplify = F)
save(glm_imps, glm_pvals, file = "data/rdata/GLMmodel/GLM_impsAndPvals.RData")

# and for the sig. glm
glm_imps_sig = sapply(tissues, function(tissue){
  load(paste0("data/rdata/GLMmodel/", tissue,
              "_importances_sig.RData"))
  temp = t(sapply(predictorOrder, function(x){
    if(x %in% rownames(imp)){
      return(imp[x,])
    } else{
      return(rep(NA, ncol(imp)))
    }
  }))
  return(temp)
}, simplify=F)
glm_pvals_sig = sapply(tissues, function(tissue){
  load(paste0("data/rdata/GLMmodel/", tissue,
              "_pvals_sig.RData"))
  pvals = as.matrix(pvals)
  temp = t(sapply(predictorOrder, function(x){
    if(x %in% rownames(pvals)){
      return(pvals[x,])
    } else{
      return(rep(NA, ncol(pvals)))
    }
  }))
  temp[temp == 0] = 2e-16
  return(temp)
}, simplify = F)
save(glm_imps_sig, glm_pvals_sig, file = "data/rdata/GLMmodel/GLM_impsAndPvals_sig.RData")
#####

#####

# load data #####
load("data/rdata/dataInfos.RData")
load("data/rdata/GLMmodel/predPerTissueGLM.RData")
load("data/rdata/GLMmodel/predPerTissueGLM_sig.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_perChr.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_perChr_sig.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat_sig.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals_sig.RData")
load("data/rdata/nPerTissue.RData")
#####

# collection of plots comparison with glm and the sig. glm #####
png("fig/modelEvaluation/compareGLMandSigGLMperformance.png", 
    height=1500, width=1300, pointsize=30)
par(mfrow = c(length(tissues),3), mar = c(0.1,4,0.1,0.1), 
    mgp = c(2.5,1,0), oma = c(3,1,0.05,0.05))
#predictorCols = setNames(rainbow(length(predictorOrder)),predictorOrder )
plotDump = sapply(tissues, function(tissue){
  print(tissue)
  # roc and auroc
  plot(ROC_PR_glm_concat[[tissue]]$roc@x.values[[1]],
       ROC_PR_glm_concat[[tissue]]$roc@y.values[[1]],lwd = 2, #
       las = 1,xlab = "",  ylab = "", type = "l", xaxt = "n",yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7)
  abline(0,1, col = "grey", lty = 2)
  plot(ROC_PR_glm_concat_sig[[tissue]]$roc, add = T,
       lwd = 2, col = "red", lty = "4414")
  mtext(side=2, text=t2T[tissue], line = 3.7)
  if(tissue == tail(tissues,1)){   
    axis(1, mgp = c(2,0.7,0), las = 1,
         at = c(0,0.5), lwd = 0, lwd.ticks=1)
    axis(1, mgp = c(2,0.7,0), las = 1,at = 1)
    title(xlab = "FPR", line = 2, xpd = NA, cex.lab = 1.2)
  }
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    title(ylab = "TPR", xpd = NA, 
          mgp = c(2,1,0), cex.lab = 1.2, line = 2)
  }
  AUCs = data.frame(GLM = sapply(ROC_PR_glm_perChr[[tissue]],function(x){x$auc}),
                    GLM_sig = sapply(ROC_PR_glm_perChr_sig[[tissue]],function(x){x$auc}))
  lim = par("plt")
  smallPlot(expr={
    boxplot(AUCs, las = 1, 
            border = c("black", "red"), xaxt = "n", yaxt = "n")
    title(ylab = "AUROC", line=2.7, cex.lab = 0.8)
    mtext(text=c("GLM", "GLM*"),side=1,  at=1:2, cex = 0.6, line=0.2,
          col=c("black", "red"), font = 1)
    axis(2, mgp = c(3,0.7,0), las = 1,cex=0.7,
         at = par("yaxp")[1:2], lwd = 0, lwd.ticks=1)
  }, 
  x1 = lim[1]+0.65*(lim[2]-lim[1]), x2 = lim[2]-0.04*(lim[2]-lim[1]),
  y1 = lim[3]+0.2*(lim[4]-lim[3]), y2 = lim[4]-0.45*(lim[4]-lim[3]), xpd = F,
  mar = c(0,0,0,0), border = "transparent")
  
  #pr
  plot(ROC_PR_glm_concat[[tissue]]$pr@x.values[[1]],
       ROC_PR_glm_concat[[tissue]]$pr@y.values[[1]], lwd = 2, 
       las = 1,xlab = "", xaxt = "n", type = "l", 
       ylab = "",  ylim = c(0,1), mgp =  c(2,0.7,0),  xaxt = "n",yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7)
  plot(ROC_PR_glm_concat_sig[[tissue]]$pr, lwd = 2, 
       col = "red", add = T, lty = "4414")
  rect(xleft = -0.02, xright=0.1, ybottom=0.5, ytop=1, lwd = 1.5)
  if(tissue == tail(tissues,1)){   
    axis(1, mgp = c(2,0.7,0), las = 1,
         at = c(0,0.5), lwd = 0, lwd.ticks=1)
    axis(1, mgp = c(2,0.7,0), las = 1,at = 1)
    title(xlab = "Recall", line = 2, xpd = NA, cex.lab = 1.2)
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    title(ylab = "Precision", xpd = NA, 
          mgp =  c(2,0.7,0), cex.lab = 1.2)
  }
  lim = par("plt")
  xlims = lim[2]-lim[1]
  ylims = lim[4]-lim[3]
  inset = 0.05
  smallPlot(expr={
    plot(x = ROC_PR_glm_concat[[tissue]]$pr@x.values[[1]],
         y = ROC_PR_glm_concat[[tissue]]$pr@y.values[[1]],
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         xlim = c(0,0.1),ylim = c(0.5,1), 
         lwd = 2, type = "l")
    lines(x = ROC_PR_glm_concat_sig[[tissue]]$pr@x.values[[1]],
          y = ROC_PR_glm_concat_sig[[tissue]]$pr@y.values[[1]],
          col = "red", lty = 2, lwd = 2)
    box(lwd = 2)
  }, 
  x1 = lim[1]+0.3*xlims, x2 = lim[2]-0.05*xlims,
  y1 = lim[3]+0.05*ylims, y2 = lim[4]-0.3*ylims, xpd = F,
  mar = c(0,0,0,0), border = "transparent")
  
  
  # violin of predictions
  glmpreds = do.call(rbind,predPerTissueGLM[[tissue]])
  glmSigpreds = do.call(rbind,predPerTissueGLMsig[[tissue]])
  glmpredssplit = split(glmSigpreds$pred, glmpreds$label)
  glmSigpredssplit = split(glmSigpreds$pred, glmpreds$label)
  plotDat = list("TN.GLM" = glmpredssplit$`0`,
                 "TN.GLMsig" = glmSigpredssplit$`0`,
                 # "gap" = c(0.5,0.6),
                 "TP.GLM" = glmpredssplit$`1`,
                 "TP.GLMsig" = glmSigpredssplit$`1`)
  sinaplot(plotDat, col = c("black", "red", "black", "red"), 
           xaxt = "n",yaxt = "n",las = 1,  ylab = "", cex = 0.8, ylim= c(0,1))
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7)
  abline(h=0.5, col = "grey", lty = 2, lwd = 2)
  abline(v=2.5)
  boxplot(plotDat, border = c("dark grey", "dark red", "dark grey", "dark red"),
          add = T, boxwex = 0.2, outline = F, col = rgb(0,0,0,alpha = 0),
          ann = F, yaxt = "n", xaxt = "n", mgp = c(2,0.7,0))
  if(tissue == tail(tissues,1)){   
    axis(1,at = c(1.5,3.5), labels=c("0", "1"), mgp = c(2,0.7,0))
    title(xlab = "True labels", line = 2, xpd = NA, cex.lab = 1.2, 
          mgp = c(2,0.7,0))
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    title(ylab = "Prediction", xpd = NA, 
          mgp =  c(2,0.7,0), cex.lab = 1.2)
  }

})
dev.off()

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
    plot(ROC_PR_glm_perChr_sig[[tissue]][[i]]$roc, 
         col = "grey70", add = T)
  })
  plot(ROC_PR_glm_concat_sig[[tissue]]$roc, 
       col = "black", add = T, lwd = 2, main = "ROC")
  
  # 2. pr over Chrs
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2.7,1,0), las = 1,
       xlab = "Recall", ylab = "Precision")
  dump = sapply(seq_len(length(chrCols)), function(i){
    plot(ROC_PR_glm_perChr_sig[[tissue]][[i]]$pr, 
         col = "grey70", add = T)
  })
  plot(ROC_PR_glm_concat_sig[[tissue]]$pr, 
       col = "black", add = T, lwd = 2, main = "PR")
  
  # 3.  AUC vs mutrate
  AUCs = sapply(ROC_PR_glm_perChr_sig[[tissue]], function(y){y$auc})
  
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
  predsConcat = do.call(rbind,predPerTissueGLMsig[[tissue]])
  sinaplot(pred ~ label, predsConcat, las = 1, 
           xlab = "True value", ylab = "Predicted value", 
           col = c("light grey", "dark grey"), maxwidth = 0.9,mgp = c(2.7,1,0))
  boxplot(pred ~ label, predsConcat, add = T, yaxt = "n",xaxt = "n",
          col = "white", boxwex = c(0.5,0.5), outline = F)
  
  # 5. plot of glm coefficients
  par(mar = c(4,17.5,0.5,0.7))
  imps = glm_imps_sig[[tissue]]
  imps = imps[apply(imps,1,function(x){!all(is.na(x))}),]
  imps = imps[!rownames(imps) %in% rarePredictors,]
  # Add * to sig. features
  #y_lab = sapply(p2P[rownames(imps)], function(name){
  #  if(name %in% p2P[names(pvals)]){
  #    return(paste0(name, " *"))
  #  } else {
  #    return(name)
  #  }
  #})
  #pvals = rowMeans(glm_pvals[[tissue]])
  #pvals = pvals[!is.na(pvals)]
  #pvals[!names(pvals) %in% rarePredictors]
  #pvals = pvals[pvals < 0.05]´
  bxplt=boxplot(t(imps), las = 1, drop =T,
                horizontal = T, cex.axis = 0.5, 
                xaxt = "n", yaxt = "n")
  axis(1)
  axis(2, at=seq(1, length(rownames(imps))), labels= p2P[rownames(imps)], las = 2)
  mtext("GLM coefficients", side = 1, cex= 2/3, line = 2.5)
  segments(y0= 1:nrow(imps), y1 = 1:nrow(imps),
           x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
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
dev.off()
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
    plot(ROC_PR_glm_perChr_sig[[tissue]][[cr]]$roc, 
         col = rgb(0,0,0,0.5), add = T)
  })
  lines(ROC_PR_glm_concat_sig[[tissue]]$roc@x.values[[1]], 
        ROC_PR_glm_concat_sig[[tissue]]$roc@y.values[[1]], 
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
    plot(ROC_PR_glm_perChr_sig[[tissue]][[cr]]$pr, 
         col = rgb(0,0,0,0.5), add = T)
  })
  plot(ROC_PR_glm_concat_sig[[tissue]]$pr, 
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

### 7. plot of coefficients --> similar to heatmaps
#png(paste0("fig/modelEvaluation/GLMCoefficientsOverChrsAllTissues_120822.png"),
#    width=1200, height=900, pointsize=20)
#par(mar = c(3,0,1.5,0.1), oma = c(0,11,0.1,0.1),mfrow = c(1,10))
#plotDump = sapply(tissues, function(tissue){
#  coeffs = glm_imps_sig[[tissue]]
#  bxplt=boxplot(t(coeffs), las = 1, drop =T,yaxt = "n", xaxt = "n",
#                horizontal = T, cex.axis = 0.75, main = t2T[tissue], 
#                xlim = c(3,nrow(coeffs)-2), names = NA, xaxs = "i")
#  axis(1,labels = F)
#  axis(1, mgp = c(2,0.5,0), cex.axis = 0.75, 
#       at = par("xaxp")[1], tick=F)
#  axis(1, mgp = c(2,0.5,0), cex.axis = 0.75, 
#       at = par("xaxp")[2], hadj = 0.8,  tick=F)
#  if(tissue == tissues[1]){
#    axis(2,at = 1:nrow(coeffs), labels=p2P[rownames(coeffs)], las = 1,
#         cex.axis = 0.75)
#  }
#  segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
#           x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
#})
#title(xlab = "Significant coefficients", 
#      xpd = NA, outer = T, line = -1.5)
#title(ylab = "Predictors", outer = T, line = 10)
dev.off()





# comparison of tissues performance #####
AUCs = sapply(tissues, function(tissue){
  ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]]
})
AUCsPerChr = sapply(tissues, function(tissue){
  sapply(ROC_PR_glm_perChr_sig[[tissue]], function(x){
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
  lines(ROC_PR_glm_concat_sig[[tissue]]$roc@x.values[[1]], 
        ROC_PR_glm_concat_sig[[tissue]]$roc@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2)
})
legend("bottomright", col=tissueCols, legend=t2T[names(tissueCols)], lwd = 2)
# pr
plot(NA, xlim = c(0,1), ylim = c(0,1), 
     xlab = "Recall", 
     ylab = "Precision", las = 1)
plotDump = sapply(tissues, function(tissue){
  lines(ROC_PR_glm_concat_sig[[tissue]]$pr@x.values[[1]], 
        ROC_PR_glm_concat_sig[[tissue]]$pr@y.values[[1]], 
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
meanPvals = do.call(rbind,lapply(tissues, function(tissue){
  imp = rowMeans(glm_pvals[[tissue]])
  d = data.frame(pvalue = imp, logpvalue = -log10(imp),tissue, direction = -log10(imp))
  d$predictor = rownames(d)
  for(f in rownames(d)){
    f_dir = mean(glm_imps[[tissue]][f,]) # check effect direction
    if(f_dir < 0 && !is.na(f_dir)){
      d[f,]$direction = -d[f,]$direction # change sign  of pvalue based on direction
    }
  }
  return(d)
}))
meanPvals$predictor = factor(meanPvals$predictor,levels=predictorOrder)
meanPvals = meanPvals[!meanPvals$predictor %in% c("followingBase_CG1", "ref_CG1", "precedingBase_CG1"),]
meanPvals["direction_cutoffs"] = meanPvals$direction
meanPvals$direction_cutoffs[ifelse(meanPvals$direction < -20, T, F)] = -20
meanPvals$direction_cutoffs[ifelse(meanPvals$direction > 20, T, F)] = 20

# removing TFs only present in liver
subMeanPvals = meanPvals[!meanPvals$predictor %in% c(rarePredictors),]

# GLM P-values
ggplot(meanPvals, aes(x = tissue, y = predictor, colour="")) + 
  geom_raster(aes(fill=direction_cutoffs)) +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(na.omit(meanPvals$direction_cutoffs)),0,max(na.omit(meanPvals$direction_cutoffs))),
    labels = c(round(min(na.omit(meanPvals$direction_cutoffs), digits = -1)), 0, round(max(na.omit(meanPvals$direction_cutoffs)), digits = -1)))+
  geom_text(aes(label = ifelse(pvalue < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Effect score")+ # labs(fill = "Log p-value &\neffect direction")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(meanPvals$predictor)])+
  scale_colour_manual(name = element_blank(), values=c("grey50", "grey50"))+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.y = element_line(colour = "grey70"),
        axis.ticks.length.x.bottom= unit(0, "cm"),
        axis.ticks.length.y.left= unit(0.1, "cm"),
        #axis.ticks.y.left = element_line(linewidth = NA),
        axis.line= element_line(linewidth = NA))
ggsave("fig/GLMpvals_heatmap.pdf", height=8, width=7, bg = "white", dpi = 200)

ggplot(subMeanPvals, aes(x = tissue, y = predictor, colour="")) + 
  geom_raster(aes(fill=direction_cutoffs)) +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(na.omit(meanPvals$direction_cutoffs)),0,max(na.omit(meanPvals$direction_cutoffs))),
    labels = c(round(min(na.omit(meanPvals$direction_cutoffs), digits = -1)), 0, round(max(na.omit(meanPvals$direction_cutoffs)), digits = -1)))+
  geom_text(aes(label = ifelse(pvalue < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Effect score")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 1, angle = 40)) + 
  scale_y_discrete(labels = p2P[as.character(subMeanPvals$predictor)])+
  scale_colour_manual(name = element_blank(), values=c("grey50", "grey50"))+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length.x.bottom= unit(.2, "cm"),
        axis.ticks = element_line(colour = "grey70"),
        axis.ticks.length.y.left= unit(0.1, "cm"),
        #axis.ticks.y.left = element_line(linewidth = NA),
        axis.line= element_line(linewidth = NA))
#ggsave("fig/GLMpvals_heatmap_noTFs.png", height=8, width=7,, bg = "white", dpi = 200)
ggsave("fig/GLMpvals_heatmap_noTFs_V3.pdf", height=8, width=7, bg = "white", dpi = 200)


# glm feature importance #####
meanImps = do.call(rbind,lapply(tissues, function(tissue){
  imp = rowMeans(glm_imps[[tissue]])
  d = data.frame(importance = imp,tissue, pval = rowMeans(glm_pvals[[tissue]]))
  d$predictor = rownames(d)
  return(d)
}))
meanImps$predictor = factor(meanImps$predictor,levels=predictorOrder)
pValsBases = meanImps[meanImps$predictor %in% c("followingBase_CG", "ref_CG", "precedingBase_CG"),]$pval
meanImps[meanImps$predictor %in% c("followingBase_CG1", "ref_CG1", "precedingBase_CG1"),]$pval = pValsBases
meanImps = meanImps[!meanImps$predictor %in% c("followingBase_CG", "ref_CG", "precedingBase_CG"),]
# removing TFs only present in liver
subMeanImps = meanImps[!meanImps$predictor %in% c(rarePredictors),]

# GLM feature importances 
ggplot(meanImps, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=importance)) +
  scale_fill_gradient2(low="blue", mid="grey90", high = "red", na.value="grey")+
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Feature\ncoefficients")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(meanImps$predictor)])+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length.x.bottom= unit(.2, "cm"),
        axis.ticks.y.left = element_line(linewidth = NA),
        axis.line= element_line(linewidth = NA))
ggsave("fig/GLMimportances_heatmap.pdf", height=8, width=7, bg = "white", dpi = 200)


ggplot(subMeanImps, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=importance)) +
  scale_fill_gradient2(low="blue", mid="grey90", high = "red", na.value="grey")+
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Feature\ncoefficients")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(meanImps$predictor)])+
  theme(plot.margin = margin(t = 0, r = 0, b = 0, l = 0), 
        panel.background = element_rect(fill = NA),
        axis.ticks.length.x.bottom= unit(.2, "cm"),
        axis.ticks.y.left = element_line(linewidth = NA),
        axis.line= element_line(linewidth = NA))
ggsave("fig/GLMimportances_heatmap_noTFs.pdf", height=8, width=7, bg = "white", dpi = 200)
#####