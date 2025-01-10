library(ggplot2)
library(data.table) 
library(ROCR)
library(ranger) 
library(berryFunctions) # for smallPlot
library(corrplot)
library(sinaplot)
library(scales)
library(RColorBrewer)
library(plotrix) # for point labels
source("scripts/05_analysis/00_NamesAndColors.R")
basesCol = c("#f18aad","#ea6759","#f88f58","#f3c65f","#8bc28c","#6667ab") #brewer.pal(n = 6, name = "Dark2")

rarePredictors = c("ZBTB33_100bp", "YY1_100bp", "TAF1_100bp", "SP1_100bp",
                   "RXRA_100bp", "REST_100bp", "RAD21_100bp",
                   "NR2F2_100bp", "MAX_100bp", "JUND_100bp", "HNF4G_100bp",
                   "HNF4A_100bp", "GABPA_100bp", "FOXA2_100bp", 
                   "FOXA1_100bp", "EGR1_100bp", "ATF3_100bp")

load("data/rdata/RFmodel/predPerTissueRF.RData")
load("data/rdata/GLMmodel/predPerTissueGLM.RData")
load("data/rdata/GLMmodel/predPerTissueGLM_sig.RData")
load("data/rdata/dataInfos.RData")
load("data/rdata/RFmodel/ROC_PR_RF_perChr.RData")
load("data/rdata/RFmodel/ROC_PR_RF_concat.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_perChr.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat.RData")
load("data/rdata/RFmodel/RF_imps.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals_sig.RData")
load("data/rdata/nPerTissue.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_perChr_sig.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat_sig.RData")
load("data/rdata/GLMperformanceCrossTissue.RData")
load("data/rdata/GLMCrossTissue_subTrain.RData")
load("data/rdata/GLMCrossTissue_subTest.RData")
load("data/rdata/GLMCrossTissue_subBoth.RData")
load("data/rdata/GLMCrossTissuePreds.RData")
load("data/rdata/GLMCrossTissueMuts.RData")
predictorOrder = read.table("scripts/04_trainmodels/predictorOrder_glm.txt")[,1]


# create overview figure nMuts per sample and mutation types #####
sortByMedian = nPerTissue[order(sapply(nPerTissue,median))]
mutClass = list(c("T>A", "A>T"),
                c("T>C", "A>G"),
                c("T>G", "A>C"),
                c("C>A", "G>T"),
                c("C>G", "G>C"),
                c("C>T", "G>A"))
names(mutClass) = sapply(mutClass, paste, collapse = "/")
mutTypes = sapply(names(sortByMedian), function(tissue){
  load(paste0("data/procData/traindata/traindata_",
              tissue, ".RData"))
  muts = data$muts[data$muts$mutated == 1,]
  temp = paste(muts$ref, muts$alt, sep=">")
  sapply(mutClass, function(x){
    sum(temp %in% x)
  })
  # table(paste(muts$ref, muts$alt, sep=">"))
})
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
trimTypes = sapply(names(sortByMedian), function(tissue){
  load(paste0("data/procData/traindata/traindata_",
              tissue, ".RData"))
  muts = data$muts[data$muts$mutated == 1,]
  temp = paste0(substr(muts$context,2,2),"[",muts$ref,">", muts$alt,"]",
                substr(muts$context,4,4))
  apply(trimClass,1, function(x){
    sum(temp %in% x)
  })
})

pdf("fig/MutationOverview.pdf",width=19, height=15, pointsize = 35)
#png("fig/MutationOverview.png", width=2000, height=1800, res=200)
layout(rbind(c(1,1,1,1,1,1,1,1,1,1),
             c(2,3,4,5,6,7,8,9,10,11),
             c(2,3,4,5,6,7,8,9,10,11),
             c(2,3,4,5,6,7,8,9,10,11)))
par(mar = c(0.1,0.15,0,0), oma = c(3.3,5,2,0))
plot(NA, xlim = c(1,20), ylim = c(1,max(unlist(nPerTissue))), log = "y",bty="l", 
     ylab = "Mutations per tumor", xaxt = "n", las = 1, xlab = "", xpd = NA,
     mgp = c(4,1,0), cex.lab = 1.1,cex.axis = 1.1)
dumpVar = sapply(1:length(sortByMedian), function(i){
  tissue= names(sortByMedian)[i]
  x = sortByMedian[[i]]
  x = sort(x)
  pos = seq(0,1,length.out = length(x))
  points(i*2-1+pos,x, col = tissueCols[tissue], pch = 1, cex=0.7)
  segments(x0=i*2-1,x1=i*2, y0 = median(x), lwd = 2, lty = 2)
  # return(cbind(i*2-1+pos,x))
})
axis(3,at =1:length(sortByMedian)*2-0.5, labels= t2T[names(sortByMedian)], tick=F,
     mgp = c(2,0.5,0), cex.axis  = 1.1)

par(mar = c(0,0.15,0,0))
temp =apply(trimTypes, 2,function(x){x/sum(x)})
sapply(1:ncol(temp), function(i){
  barplot(rev(temp[,i]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
          names.arg="",cex.names=0.5, cex.axis = 1.1,
          col = rep(basesCol, each=16), xlim = c(0,0.16), xaxp = c(0, 0.1, 1))
  if(i==1){
    classes = trimClass[rownames(temp),1]
    names = rev(paste0(substr(classes,1,1), substr(classes,3,3), substr(classes,7,7)))
    names2 = rev(paste0(" ", substr(classes,3,3), " "))
    text(x = 0, y = 0.99:96-0.5,  labels=names,family = "mono",
         xpd = NA, adj=1.1,  cex = 0.5)
    text(x = 0, y = 0.99:96-0.5,  labels=names2,family = "mono",
         xpd = NA, adj=1.1, col = rep(basesCol, each=16), font = 2, cex = 0.5)
    text(x = -0.055, y = 1:6*16-8,
         labels = c("T>G", "T>C", "T>A" ,"C>T", "C>G", "C>A"),
         xpd = NA, adj=1, col = basesCol, font = 2)
    segments(y0=1:6*16-16, y1 = 1:6*16, x0=-0.045, lwd = 4, 
             col = basesCol, xpd = NA, lend=1, cex = 1.1)
    title(ylab = "Mutation type", mgp = c(4,1,0), xpd = NA, cex.lab = 1.1)
  }
  abline(h=seq(16,80,length.out = 5), lty = 2, lwd = 1.5)
  abline(h=0, lty = 1, lwd = 2)
  #abline(h=96, lty = 1, lwd = 1.8)
})
title(xlab = "Proportion of variants", outer = T, mgp = c(2,1,0), cex.lab = 1.1)
dev.off()
#####

# Compare GLM and RF #####
#pdf("fig/modelEvaluation/compareRFandGLMperformance.pdf",width=16, height=20, pointsize = 34)
png("fig/modelEvaluation/compareRFandGLMperformance.png", 
    height=1600, width=1300, pointsize=30)#, res = 200)
par(mfrow = c(length(tissues),3), mar = c(0.5,4,0.1,0.1), 
    mgp = c(2.5,1,0), oma = c(3,1,0.05,0.05))
#predictorCols = setNames(rainbow(length(predictorOrder)),predictorOrder )
plotDump = sapply(tissues, function(tissue){
  print(tissue)
  # roc and auroc
  plot(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
       ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]],lwd = 2, #
       las = 1,xlab = "",  ylab = "", type = "l", xaxt = "n",yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7, cex.axis = 1.4)
  abline(0,1, col = "grey", lty = 2)
  plot(ROC_PR_glm_concat_sig[[tissue]]$roc, add = T,
       lwd = 2, col = "red", lty = "4414")
  mtext(side=2, text=t2T[tissue], line = 3.7)
  if(tissue == tail(tissues,1)){   
    axis(1, mgp = c(2,0.7,0), las = 1,
         at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis=1.4)
    axis(1, mgp = c(2,0.7,0), las = 1,at = 1, cex.axis=1.4)
    title(xlab = "FPR", line = 2, xpd = NA, cex.lab = 1.4)
  }
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    title(ylab = "TPR", xpd = NA, mgp = c(2.5,1,0), cex.lab = 1.4)
  }
  AUCs = data.frame(RF = sapply(ROC_PR_RF_perChr[[tissue]],function(x){x$auc}),
                    GLM = sapply(ROC_PR_glm_perChr_sig[[tissue]],function(x){x$auc}))
  lim = par("plt")
  smallPlot(expr={
    boxplot(AUCs, las = 1, 
            border = c("black", "red"), xaxt = "n", yaxt = "n")
    title(ylab = "AUROC", line=2.5)
    mtext(text=c("RF", "GLM"),side=1,  at=1:2, cex = 0.55, line=0.2,
          col=c("black", "red"), font = 2)
    axis(2, mgp = c(3,0.7,0), las = 1,cex=0.7,
         at = par("yaxp")[1:2], lwd = 0, lwd.ticks=1)
  }, 
  x1 = lim[1]+0.72*(lim[2]-lim[1]), x2 = lim[2]-0.02*(lim[2]-lim[1]),
  y1 = lim[3]+0.2*(lim[4]-lim[3]), y2 = lim[4]-0.45*(lim[4]-lim[3]), xpd = F,
  mar = c(0,0,0,0), border = "transparent")
  
  #pr
  plot(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
       ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], lwd = 2, 
       las = 1,xlab = "", xaxt = "n", type = "l", 
       ylab = "",  ylim = c(0,1), mgp =  c(2,0.7,0),  xaxt = "n",yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7, cex.axis = 1.4)
  plot(ROC_PR_glm_concat_sig[[tissue]]$pr, lwd = 2, 
       col = "red", add = T, lty = "4414")
  rect(xleft = -0.02, xright=0.1, ybottom=0.5, ytop=1, lwd = 1.5)
  if(tissue == tail(tissues,1)){   
    axis(1, mgp = c(2,0.7,0), las = 1,
         at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
    axis(1, mgp = c(2,0.7,0), las = 1,at = 1, cex.axis = 1.4)
    title(xlab = "Recall", line = 2, xpd = NA, cex.lab = 1.4)
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    title(ylab = "Precision", xpd = NA, 
          mgp =  c(2.5,0.7,0), cex.lab = 1.4)
  }
  lim = par("plt")
  xlims = lim[2]-lim[1]
  ylims = lim[4]-lim[3]
  inset = 0.05
  smallPlot(expr={
    plot(x = ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
         y = ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]],
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
  rfpreds = do.call(rbind,predPerTissueRF[[tissue]])
  glmpreds = do.call(rbind,predPerTissueGLMsig[[tissue]])
  rfpredssplit = split(rfpreds$pred, rfpreds$label)
  glmpredssplit = split(glmpreds$pred, rfpreds$label)
  plotDat = list("TN.RF" = rfpredssplit$`0`,
                 "TN.GLM" = glmpredssplit$`0`,
                 # "gap" = c(0.5,0.6),
                 "TP.RF" = rfpredssplit$`1`,
                 "TP.GLM" = glmpredssplit$`1`)
  sinaplot(plotDat, col = c("black", "red", "black", "red"), #
           xaxt = "n",yaxt = "n",las = 1,  ylab = "", cex = 0.8, ylim= c(0,1))
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7, cex.axis = 1.4)
  abline(h=0.5, col = "grey", lty = 2, lwd = 2)
  abline(v=2.5)
  boxplot(plotDat, border = c("dark grey", "dark red", "dark grey", "dark red"),
          add = T, boxwex = 0.2, outline = F, col = rgb(0,0,0,alpha = 0),
          ann = F, yaxt = "n", xaxt = "n", mgp = c(2,0.7,0))
  if(tissue == tail(tissues,1)){   
    axis(1,at = c(1.5,3.5), labels=c("0", "1"), mgp = c(2,0.7,0), cex.axis=1.4)
    title(xlab = "True labels", line = 2, xpd = NA, 
          mgp = c(2,0.7,0), cex.lab = 1.4)
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    title(ylab = "Prediction", xpd = NA, 
          mgp =  c(2.5,0.7,0), cex.lab = 1.4)
  }
})
dev.off()  
#####

# ROC, PR and nMuts overview #####
#pdf(paste0("fig/modelEvaluation/GLM_ROCsOverChrsAllTissues.pdf"),
#    width=17, height=5, pointsize = 26)
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


#pdf(paste0("fig/modelEvaluation/GLM_PRsOverChrsAllTissues.pdf"),
#    width=17, height=5, pointsize = 26)
png(paste0("fig/modelEvaluation/GLM_PRsOverChrsAllTissues_301122.png"),
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

# collection of plots for every tissue #####
plotDump = sapply(tissues, function(tissue){
  #pdf(paste0("fig/modelEvaluation/perfGLMsummary_", tissue, ".pdf"),
  #    width=12, height=10, pointsize = 20)
  png(paste0("fig/modelEvaluation/perfGLMsummary_", tissue, ".png"),
      width=1650, height=1300, pointsize = 28)
  m = rbind(c(1,1,2,2,5,5,5,5), 
            c(1,1,2,2,5,5,5,5), 
            c(3,3,4,4,5,5,5,5), 
            c(3,3,4,4,5,5,5,5), 
            c(6,6,6,6,5,5,5,5),
            c(6,6,6,6,5,5,5,5))
  par(oma = c(0.1,0.1,1.3,0.1), mar = c(4,4,0.5,0))
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
  #percMut = dataInfos[[tissue]]$percTPperChr
  #plot(percMut, AUCs, las = 1, 
  #     ylab = "AUROC", xlab = "percent TP", mgp = c(2.7,1,0))
  #mod = lm(AUCs~percMut)
  #abline(mod, col = "grey", lty = 2, lwd = 2)
  #p = summary(mod)$coefficients[2,4]
  #r2 = summary(mod)$adj.r.squared
  #rp = vector('expression',2)
  #rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
  #                   list(MYVALUE = format(r2,dig=3)))[2]
  #rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
  #                   list(MYOTHERVALUE = format(p, digits = 2)))[2]
  
  #legend('bottomright', legend = rp, bty = 'n')
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
  par(mar = c(4,17.5,0,0.7))
  imps = glm_imps[[tissue]]
  imps = imps[apply(imps,1,function(x){!all(is.na(x))}),]
  imps = imps[!rownames(imps) %in% rarePredictors,]
  
  # color-code for sig. features 
  pvals = rowMeans(glm_pvals[[tissue]])
  pvals = pvals[!is.na(pvals)]
  pvals[!names(pvals) %in% rarePredictors]
  pvals = pvals[pvals < 0.05]
  #pval_colors = rep(NA, nrow(imps))
  #names(pval_colors)=p2P[rownames(imps)]
  #pval_colors[p2P[names(pvals)]] = "tomato"
  
  bxplt=boxplot(t(imps), las = 1, drop =T,
                horizontal = T, cex.axis = 0.5, 
                xaxt = "n", yaxt = "n")#, col = pval_colors)
  axis(1)
  
  # Add * to sig. features
  y_lab = sapply(p2P[rownames(imps)], function(name){
    if(name %in% p2P[names(pvals)]){
      return(paste0(name, " *"))
    } else {
      return(name)
    }
  })
  axis(2, at=seq(1, length(rownames(imps))), labels= y_lab, las = 2)
  mtext("GLM coefficients", side = 1, cex= 0.65, line = 2.5)
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
  par(mar = c(4,8,0.5,0))
  barplot(rev(mutDat[,1]), horiz=T, col = basesCol,
          las = 1, xlim = c(-max(mutDat),max(mutDat)),names.arg=rev(labels),
          cex.lab=1,cex.names=1, xlab = "count", xaxt = "n")
  barplot(rev(mutDat[,2]), add = T, horiz = T, col = basesCol,xaxt = "n", yaxt = "n")
  temp = axTicks(1)
  axis(1, at = temp, labels=abs(temp))
  title(t2T[tissue], outer = T)
  dev.off()
})
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

ggplot(subMeanPvals, aes(x = tissue, y = predictor, colour="")) + 
  geom_raster(aes(fill=direction_cutoffs)) +
  #scale_fill_viridis_c(option = "cividis") +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(na.omit(meanPvals$direction_cutoffs)),0,max(na.omit(meanPvals$direction_cutoffs))),
    labels = c(round(min(na.omit(meanPvals$direction_cutoffs), digits = -1)), 0, round(max(na.omit(meanPvals$direction_cutoffs)), digits = -1)))+
  geom_text(aes(label = ifelse(pvalue < 0.05, yes = "*", no = " ")), nudge_y = -0.4, col = "gray7")+
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Effect score")+ # labs(fill = "Log p-value &\neffect direction")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(subMeanPvals$predictor)])+
  scale_colour_manual(name = element_blank(), values=c("grey50", "grey50"))+
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))  + theme_minimal()
#ggsave("fig/GLMpvals_heatmap_noTFs_120822.png", height=8, width=7,, bg = "white", dpi = 200)
ggsave("fig/GLMpvals_heatmap_noTFs.pdf", height=8, width=7, bg = "white", dpi = 200) # Add grey box to legend by hand 
######

# cross-tissue application #####
# collect ROCs
ROCs = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = GLMperformanceCrossTissue[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs$AUROC = as.numeric(ROCs$AUROC)
ROCs$rank <- NA

ROCs_subTrain = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subTrain[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subTrain$AUROC = as.numeric(ROCs_subTrain$AUROC)
ROCs_subTrain$rank = NA

ROCs_subTest = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subTest[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subTest$AUROC = as.numeric(ROCs_subTest$AUROC)
ROCs_subTest$rank = NA

ROCS_subBoth = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subBoth[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCS_subBoth$AUROC = as.numeric(ROCS_subBoth$AUROC)
ROCS_subBoth$rank = NA

RocsPreds = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
  t(sapply(tissues, function(predtissue){
    temp = CrossTissuePreds[[muttissue]][[predtissue]]$auc@y.values[[1]]
    return(c("predictorSource" = predtissue, 
             "MutationAndModelSource" = muttissue, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
RocsPreds$AUROC = as.numeric(RocsPreds$AUROC)
RocsPreds$rank <- NA

RocsMuts = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
  t(sapply(tissues, function(predtissue){
    temp = CrossTissueMuts[[muttissue]][[predtissue]]$auc@y.values[[1]]
    return(c("MutationSource" = predtissue, 
             "PredictorAndModelSource" = muttissue, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
RocsMuts$AUROC = as.numeric(RocsMuts$AUROC)
RocsMuts$rank <- NA

load("data/rdata/GLMCrossTissue_subTrainNsamples.RData")
load("data/rdata/GLMCrossTissue_subTestNsamples.RData")
load("data/rdata/GLMCrossTissue_subBothNsamples.RData")
tissues = c("brain","breast", "colon","kidney","liver", "lung","ovary", #"esophagus",
            "prostate", "skin")
ROCs_subNTrain = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subTrain[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subNTrain$AUROC = as.numeric(ROCs_subNTrain$AUROC)
ROCs_subNTrain$rank = NA

ROCs_subNTest = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subTest[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subNTest$AUROC = as.numeric(ROCs_subNTest$AUROC)
ROCs_subNTest$rank = NA

ROCS_subNBoth = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = subBoth[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUROC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCS_subNBoth$AUROC = as.numeric(ROCS_subNBoth$AUROC)
ROCS_subNBoth$rank = NA

# Add rowrank to each dataset
tissues = c("brain","breast", "colon","kidney","liver", "lung","ovary", "esophagus",
            "prostate", "skin")
for(tissue in tissues){
  ROCs[ROCs$testedOn == tissue,]$rank[order(-ROCs[ROCs$testedOn == tissue,]$AUROC)] = 1:length(tissues)
  ROCs_subTrain[ROCs_subTrain$testedOn == tissue,]$rank[order(-ROCs_subTrain[ROCs_subTrain$testedOn == tissue,]$AUROC)] = 1:length(tissues)
  ROCs_subTest[ROCs_subTest$testedOn == tissue,]$rank[order(-ROCs_subTest[ROCs_subTest$testedOn == tissue,]$AUROC)] = 1:length(tissues)
  ROCS_subBoth[ROCS_subBoth$testedOn == tissue,]$rank[order(-ROCS_subBoth[ROCS_subBoth$testedOn == tissue,]$AUROC)] = 1:length(tissues)
  RocsPreds[RocsPreds$MutationAndModelSource == tissue,]$rank[order(-RocsPreds[RocsPreds$MutationAndModelSource == tissue,]$AUROC)] = 1:length(tissues)
  RocsMuts[RocsMuts$MutationSource == tissue,]$rank[order(-RocsMuts[RocsMuts$MutationSource == tissue,]$AUROC)] = 1:length(tissues)
  if(tissue != "esophagus"){
    ROCs_subNTrain[ROCs_subNTrain$testedOn == tissue,]$rank[order(-ROCs_subNTrain[ROCs_subNTrain$testedOn == tissue,]$AUROC)] = 1:length(tissues)
    ROCs_subNTest[ROCs_subNTest$testedOn == tissue,]$rank[order(-ROCs_subNTest[ROCs_subNTest$testedOn == tissue,]$AUROC)] = 1:length(tissues)
   ROCS_subNBoth[ROCS_subNBoth$testedOn == tissue,]$rank[order(-ROCS_subNBoth[ROCS_subNBoth$testedOn == tissue,]$AUROC)] = 1:length(tissues)
  }
}

# Plotting
upperL = 0.69
lowerL = 0.48
ggplot(ROCs, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)
#ggsave(filename="fig/CrossTissue_191022.png", 
#       width = 12, height = 12, units = "cm", bg = "white", dpi = 200)
ggsave(filename="fig/CrossTissue.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(RocsPreds, aes(x = predictorSource, y = MutationAndModelSource, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  labs(x = "Predictors tested on", y = "Model & mutations trained and tested on") + 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)
#ggsave(filename="fig/CrossTissuePreds_glm_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissuePreds.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(RocsMuts, aes(y = MutationSource, x = PredictorAndModelSource, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  labs(y = "Mutations tested on", x = "Model & predictors trained and tested on") + 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)
#ggsave(filename="fig/CrossTissueMuts_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissueMuts.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCs_subTrain, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)
#ggsave(filename="fig/CrossTissue_subTrain_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subTrain.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCs_subTest, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)
#ggsave(filename="fig/CrossTissue_subTest_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subTest.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCS_subBoth, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)
#ggsave(filename="fig/CrossTissue_subBoth_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subBoth.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCs_subNTrain, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)
#ggsave(filename="fig/CrossTissue_subTrainNsamples_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subTrainNsamples.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCs_subNTest, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)
#ggsave(filename="fig/CrossTissue_subTestNsamples_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subTestNsamples.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)

ggplot(ROCS_subNBoth, aes(trainedOn, testedOn, fill= AUROC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=c(lowerL,upperL)) +
  theme_minimal() +
  geom_text(aes(label = rank), color = "white", size = 2.5)+
  labs(y = "Mutations and predictors tested on", x = "Model trained on") +  
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = t2T)
#ggsave(filename="fig/CrossTissue_subBothNsamples_191022.png", 
#       width = 12, height = 12, units = "cm")
ggsave(filename="fig/CrossTissue_subBothNsamples.pdf", 
       width = 13, height = 10, units = "cm", bg = "white", dpi = 200)
