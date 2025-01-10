library(ggplot2)
library(data.table) 
library(ROCR)
library(ranger) 
library(berryFunctions) # for smallPlot
library(corrplot)
library(sinaplot)
library(plotrix) # for point labels

dir.create("fig/performance", showWarnings=F)

tissues = c("luad", "skin", "colon", "ovary",
            "kidney", "prostate", "breast")
tissueCols = setNames(rainbow(length(tissues)), tissues)
chrCols = setNames(rainbow(22), sort(paste0("chr", 1:22)))
# manually sort the importances by function of the predictors 
predictorOrder = {c("ref", "strand", "precedingBase","followingBase" ,
                    "aPhasedRepeats_100bp","directRepeats_100bp",
                    "shortTandemRepeats_100bp","gQuadruplex_100bp",
                    "mirrorRepeats_100bp","zDNAmotifs_100bp",
                    "repeatMasker", "Trf", "mappability_100mer",
                    "mappability_24mer", "mappability_40mer",
                    "conservation_100bp","GTEx_eqtl",
                    "meth" , "methylation_bins10kb",  
                    "GCcontent",  "GCcontent_bins10kb",
                    "replDirection", "replSlope", "replTiming",  
                    "HiC_inTAD", "HiCints", "HiC_ints","HiC_compLabels", 
                    "HiC_compPCA", "HiC_FIRE", "HiC_TADbound", 
                    "DNAaccessibility_tissue_100bp",   
                    "DNAaccessibility_tissue_bins10kb",
                    "DNAaccessibility_UCSC_bins10kb",
                    "NsomeGm12878_100bp", "NsomeGm12878_bins10kb", 
                    "NsomeK562_100bp", "NsomeK562_bins10kb", 
                    "H3K27ac_100bp", "H3K27ac_bins10kb", "UCSC_H3k27ac_bins10kb",              
                    "H3K27me3_100bp","H3K27me3_bins10kb",
                    "H3K36me3_100bp", "H3K36me3_bins10kb",
                    "H3K4me1_100bp", "H3K4me1_bins10kb","UCSC_H3k4me1_bins10kb",
                    "H3K4me3_100bp", "H3K4me3_bins10kb","UCSC_H3k4me3_bins10kb",
                    "H3K9me3_100bp", "H3K9me3_bins10kb", 
                    "H3K9ac_100bp", "H3K9ac_bins10kb",
                    "EP300_100bp", "EP300_bins10kb",
                    "cancerExpr","healthyExpr",
                    "POLR2A_100bp", "POLR2A_bins10kb", 
                    "POLR2AphosphoS5_100bp", "POLR2AphosphoS5_bins10kb", 
                    "TF_BS_100bp", "TF_BS_bins10kb", 
                    "ETS_BS_100bp", "ETS_BS_bins10kb",
                    "CTCF_100bp", "CTCF_bins10kb")} 
predictorCols = setNames(rainbow(length(predictorOrder)), predictorOrder)


# load RF predictions for each tissue #####
predPerTissue = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/RFmodel/predictions_withBinWise.RData"))
   return(predictions)
}, simplify=F)
save(predPerTissue, file = "data/rdata/predPerTissue.RData")
#####
# get mutation rate and nPositions #####
dataInfos = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
   chrCount = as.integer(table(removed$chr))
   chrPerc = sapply(split(data$mutated, removed$chr),function(x){mean(x==1)})
   data$inexon = NULL
   cors = cor(data[sapply(data, is.numeric)], use = "pair")
   return(list(nMuts = nrow(data), 
               percTP = mean(data$mutated == 1), 
               nMutsPerChr = chrCount,
               percTPperChr = chrPerc, 
               cors = cors))
}, simplify = F)
save(dataInfos, file = "data/rdata/dataInfos.RData")
#####
# compute ROC, PR, and AUROC for each chromosome #####
ROC_PR_RF_perChr = sapply(names(predPerTissue), function(tissue){ #iterate through tissues
   pred = predPerTissue[[tissue]]
   # get performances
   res = lapply(pred, function(x){ #iterate through chromosomes
      perf = prediction(x$pred_1, x$labels)
      roc = performance(perf, "tpr", "fpr")
      pr = performance(perf,"prec", "rec")
      auc = performance(perf,"auc")@y.values[[1]]
      return(list(roc = roc, pr = pr, auc = auc))
   })
   return(res)
}, simplify=F)
save(ROC_PR_RF_perChr, file = "data/rdata/ROC_PR_RF_perChr.RData")
#####
# compute ROC, PR, and AUROC for all chromosomes concatenated #####
ROC_PR_RF_concat = sapply(names(predPerTissue), function(tissue){
   predConcat = do.call(rbind,predPerTissue[[tissue]])
   ROC_rf = performance(prediction(predConcat$pred_1, 
                                   predConcat$labels), 
                        "tpr", "fpr")
   AUC_rf = performance(prediction(predConcat$pred_1, 
                                   predConcat$labels), 
                        "auc")
   PR_rf = performance(prediction(predConcat$pred_1, 
                                  predConcat$labels), 
                       "prec", "rec")
   return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_RF_concat, file = "data/rdata/ROC_PR_RF_concat.RData")
#####
# get ROC, PR, and AUROC for each chr from glm #####
ROC_PR_glm_perChr = sapply(names(predPerTissue), function(tissue){
   aucList = get(load(paste0("data/nina/rdata_withBinwise/", tissue,
                         "/data_glm/performROC_AUC_glm_Loco-CV.RData")))
   rocList = get(load(paste0("data/nina/rdata_withBinwise/", tissue,
               "/data_glm/performROC_glm_Loco-CV.RData")))
   prList = get(load(paste0("data/nina/rdata_withBinwise/", tissue,
                         "/data_glm/performPR_glm_Loco-CV.RData")))
   perf = sapply(names(aucList), function(x){ #iterate through chromosomes
      return(list(roc = rocList[[x]], 
                  pr = prList[[x]], 
                  auc = aucList[x]))
   }, simplify=F)
   return(perf)
}, simplify=F)
save(ROC_PR_glm_perChr, file = "data/rdata/ROC_PR_glm_perChr.RData")
#####
# get ROC, PR, and AUROC concatenated from glm #####
ROC_PR_glm_concat = sapply(tissues, function(tissue){
   ROC = get(load(paste0("data/nina/rdata_withBinwise/", tissue, 
                     "/data_glm/performROC_concat_glm.RData")))
   PR = get(load(paste0("data/nina/rdata_withBinwise/", tissue, 
                    "/data_glm/performPR_concat_glm.RData")))
   AUC = get(load(paste0("data/nina/rdata_withBinwise/", tissue, 
                     "/data_glm/performAUC_concat_glm.RData")))
   return(list(roc = ROC, pr = PR, auc = AUC))
}, simplify=F)
save(ROC_PR_glm_concat, file = "data/rdata/ROC_PR_glm_concat.RData")
#####
# get predictor importances and pvals from rf#####
rf_imps = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, 
               "/RFmodel/imp_withBinWise.RData"))
   t(sapply(predictorOrder, function(i){
      if(i %in% rownames(imp)){(imp[i,])} else{rep(NA,ncol(imp))}
   }))
}, simplify=F)
rf_pvals = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, 
               "/RFmodel/imp_pvals_withBinWise.RData"))
   sapply(predictorOrder, function(i){
      if(i %in% rownames(pvals)){(pvals[i,2])} else{NA}
   })
})
save(rf_imps, rf_pvals, file = "data/rdata/RF_impsAndPvals.RData")
#####
# get predictor importances and pvals from glm #####
glm_imps = sapply(tissues, function(tissue){
   load(paste0("data/nina/rdata_withBinwise/", tissue, 
               "/data_glm/coefs_allChr_glm.RData"))
   coefs_allChr = as.matrix(coefs_allChr)
   coefsSorted = t(sapply(predictorOrder, function(i){
      if(i %in% rownames(coefs_allChr)){
         (coefs_allChr[i,])
      } else if(i %in% c("ref", "precedingBase", "folloringBase")){ 
         temp = coefs_allChr[grep(i,row.names(coefs_allChr)),]
         colMeans(temp)
      } else{
         setNames(rep(NA,ncol(coefs_allChr)),
                           colnames(coefs_allChr)) 
      }
   }))
   return(coefsSorted)
}, simplify=F)
glm_pvals = sapply(tissues, function(tissue){
   load(paste0("data/nina/rdata_withBinwise/", tissue,
               "/data_glm/pValues_features_allChr_glm.RData"))
   pValue_features_allChr = as.matrix(pValue_features_allChr)
   pvalsSorted = t(sapply(predictorOrder, function(i){
      if(i %in% rownames(pValue_features_allChr)){
         (pValue_features_allChr[i,])
      } else{
         rep(NA,ncol(pValue_features_allChr))
      }
   }))
   pvalsSorted[pvalsSorted == 0] = 2.2e-16
   return(pvalsSorted)
}, simplify=F)
save(glm_imps, glm_pvals, file = "data/rdata/GLM_impsAndPvals.RData")
#####

# load #####
load("data/rdata/predPerTissue.RData")
load("data/rdata/dataInfos.RData")
load("data/rdata/ROC_PR_RF_perChr.RData")
load("data/rdata/ROC_PR_RF_concat.RData")
load("data/rdata/ROC_PR_glm_perChr.RData")
load("data/rdata/ROC_PR_glm_concat.RData")
load("data/rdata/RF_impsAndPvals.RData")
load("data/rdata/GLM_impsAndPvals.RData")
#####


# collection of plots for every tissue #####
plotDump = sapply(tissue, function(tissue){
   png(paste0("fig/performance/perfRFsummary_", tissue, ".png"),
       width=1200, height=1000, pointsize = 25)
   m = rbind(c(1,2,7,7), c(3,4,7,7), c(5,6,7,7))
   par(oma = c(0.1,0.1,1.3,0.1), mar = c(4,4,0.5,0))
   layout(m)
   # roc over Chrs 
   plot(NA, xlim = c(0,1), ylim = c(0,1),
        mgp = c(2.7,1,0), las = 1,
        xlab = "False positive rate", ylab = "True positive rate")
   abline(0,1, col = "grey", lty = 2)
   dump = sapply(seq_len(length(chrCols)), function(i){
      plot(ROC_PR_RF_perChr[[tissue]][[i]]$roc, 
           col = chrCols[i], add = T)
   })
   plot(ROC_PR_RF_concat[[tissue]]$roc, 
        col = "black", add = T, lwd = 2)
   # pr over Chrs
   plot(NA, xlim = c(0,1), ylim = c(0,1),
        mgp = c(2.7,1,0), las = 1,
        xlab = "Recall", ylab = "Precision")
   dump = sapply(seq_len(length(chrCols)), function(i){
      plot(ROC_PR_RF_perChr[[tissue]][[i]]$pr, 
              col = chrCols[i], add = T)
   })
   plot(ROC_PR_RF_concat[[tissue]]$pr, 
        col = "black", add = T, lwd = 2)
   # TP vs TN preds violin plot
   predsConcat = do.call(rbind,predPerTissue[[tissue]])
   sinaplot(pred_1 ~ labels, predsConcat, las = 1, 
            xlab = "True value", ylab = "Predicted value", 
            col = c("grey", "black"), maxwidth = 0.9,mgp = c(2.7,1,0))
   # legend
   plot.new()
   legend("center", legend=c(names(chrCols), "concat"), cex = 0.75,
          col = c(chrCols, "black"), lty=1, ncol=2, bty = "n")
   # AUC vs mutrate
   AUCs = sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc})
   percMut = dataInfos[[tissue]]$percTPperChr
   plot(percMut, AUCs, las = 1, col = tissueCols[tissue],
        ylab = "AUC", xlab = "percent TP", mgp = c(2.7,1,0))
   mod = lm(AUCs~percMut)
   abline(mod)
   p = summary(mod)$coefficients[2,4]
   legend("topright", bty = "n",
          legend=paste0("p=", format(p,digits=2)))
   # AUC vs nPos 
   nMuts = dataInfos[[tissue]]$nMutsPerChr
   plot(nMuts, AUCs, las = 1, col = tissueCols[tissue],mgp = c(2.7,1,0),
        ylab = "AUC", xlab = "*1000 positions per chr")
   mod = lm(AUCs~nMuts)
   abline(mod)
   p = summary(mod)$coefficients[2,4]
   legend("topright", bty = "n",
          legend=paste0("p=", format(p, digits=2)))
   # plot of coefficients
   par(mar = c(4,11,0,0))
   coeffs = rf_imps[[tissue]]
   pvals = rf_pvals[,tissue]
   coeffs = coeffs[!is.na(pvals),]
   pvals = pvals[!is.na(pvals)]
   bxplt=boxplot(t(coeffs), las = 1, drop =T,
                col = c("white", "dark grey")[(pvals<0.05)+1],
                horizontal = T, cex.axis = 0.75, xaxt = "n",
                xlim = c(2,nrow(coeffs)-1))
   axis(1)
   mtext("permutation importance", side = 1, cex= 0.8, line = 2.5)
   segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
             x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
   legend("topright", fill = c("white", "dark grey"),
          legend = c("not significant", "significant"), bty = "n")
   title(tissue, outer = T)
   dev.off()
})
#####


# collection of plots comparison with glm #####
# tissues x measures
png("fig/performance/compareRFandGLMperformance.png", 
    height=1200, width=1200, pointsize=30)
par(mfrow = c(7,4), mar = c(0.1,4,0.1,0.1), 
    mgp = c(2.5,1,0), oma = c(3,1,0.05,0.05))
plotDump = sapply(tissues, function(tissue){
   print(tissue)
   # roc
   plot(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
        ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]],lwd = 2, #
        las = 1,xlab = "",  ylab = "TPR", type = "l")
   abline(0,1, col = "grey")
   plot(ROC_PR_glm_concat[[tissue]]$roc, add = T,
        lwd = 2, col = "red", lty = "4414")
   mtext(side=2, text=tissue, line = 3.7)
   if(tissue == tail(tissues,1)){   
      axis(1)
      title(xlab = "FPR", line = 2, xpd = NA)
   }
   #pr
   plot(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
        ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], lwd = 2, 
        las = 1,xlab = "", xaxt = "n", type = "l", 
        ylab = "Precision",  ylim = c(0,1))
   plot(ROC_PR_glm_concat[[tissue]]$pr, lwd = 2, 
        col = "red", add = T, lty = "4414")
   rect(xleft = -0.02, xright=0.1, ybottom=0.5, ytop=0.9, lwd = 1.5)
   if(tissue == tail(tissues,1)){   
      axis(1)
      title(xlab = "Recall", line = 2, xpd = NA)
   }
   lim = par("plt")
   inset = 0.05
   smallPlot(expr={
      plot(x = ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
           y = ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]],
           xlab = "", ylab = "", xaxt = "n", yaxt = "n",
           xlim = c(0,0.1),ylim = c(0.5,0.9), 
           lwd = 2, type = "l")
      lines(x = ROC_PR_glm_concat[[tissue]]$pr@x.values[[1]],
            y = ROC_PR_glm_concat[[tissue]]$pr@y.values[[1]],
            col = "red", lty = 2, lwd = 2)
      box(lwd = 2)
   }, 
   x1 = lim[1]+0.3*(lim[2]-lim[1]), x2 = lim[2]-inset,
   y1 = lim[3]+0.3*(lim[2]-lim[1]), y2 = lim[4]-inset, xpd = F,
   mar = c(0,0,0,0), border = "transparent")
   
   # boxplot of AUCs
   AUCs = data.frame(RF = sapply(ROC_PR_RF_perChr[[tissue]],function(x){x$auc}),
                GLM = sapply(ROC_PR_glm_perChr[[tissue]],function(x){x$auc}))
   boxplot(AUCs, las = 1, ylab = "AUROC", 
            border = c("black", "red"), xaxt = "n")
   sinaplot(AUCs, las = 1, ylab = "AUROC", add = T,
           col = c("black", "red"), xaxt = "n")
   if(tissue == tail(tissues,1)){   
      axis(1, at=1:2, labels=c("RF", "GLM"))
   }
   # violin of predictions
   rfpreds = do.call(rbind,predPerTissue[[tissue]])[,2:3]
   glmpreds = do.call(c,sapply(names(chrCols), function(chrom){
      get(load(paste0("data/nina/rdata_withBinwise/", tissue, "/",
                     chrom,"/yhat_",chrom, "_glm.RData"))) }))
   rfpredssplit = split(rfpreds$pred_1, rfpreds$labels)
   glmpredssplit = split(glmpreds, rfpreds$labels)
   plotDat = list("TN.RF" = rfpredssplit$`0`,
               "TN.GLM" = glmpredssplit$`0`,
               # "gap" = c(0.5,0.6),
               "TP.RF" = rfpredssplit$`1`,
               "TP.GLM" = glmpredssplit$`1`)
   sinaplot(plotDat, col = c("black", "red", "black", "red"), #
            xaxt = "n",las = 1,  ylab = "Prediction")
   abline(h=0.5, col = "grey", lty = 2, lwd = 2)
   abline(v=2.5)
   boxplot(plotDat, border = c("dark grey", "dark red", "dark grey", "dark red"),
           add = T, boxwex = 0.2, outline = F, col = rgb(0,0,0,alpha = 0),
           ann = F, yaxt = "n", xaxt = "n")
   if(tissue == tail(tissues,1)){   
      axis(1,at = c(1.5,3.5), labels=c("0", "1"))
      title(xlab = "True labels", line = 2, xpd = NA)
   }   
})
dev.off()

png("fig/performance/compareRFandGLMcoefficients.png", 
    height=1200, width=400, pointsize=20)
par(mfrow = c(7,3), mar = c(1.5,3.5,0.1,0.1), oma = c(3,1,0.05,0.05))
plotDump = sapply(tissues, function(tissue){
   #coefficients
   rf_imp = rf_imps[[tissue]]
   glm_imp = abs(glm_imps[[tissue]])
   plot(rf_imp, glm_imp, las = 1, pch = 1, cex = 0.5,
        col = NA, # predictorCols[rownames(rf_imp)],
        ylab = "", mgp = c(2,0.4,0), tcl=-0.3) #, xaxt = "n"
   rfMeans = rowMeans(rf_imp, na.rm = T)
   glmMeans= rowMeans(glm_imp, na.rm = T)
   rfSD = apply(rf_imp,1,sd, na.rm = T)
   glmSD = apply(glm_imp,1,sd, na.rm = T)
   arrows(x0 = rfMeans-rfSD, x1 = rfMeans+rfSD, 
          y0 = glmMeans, y1 = glmMeans, code=0)
   arrows(x0 = rfMeans, x1 = rfMeans, 
          y0 = glmMeans-glmSD, y1 = glmMeans+glmSD, code=0)
   points(rfMeans, glmMeans,
          bg = predictorCols[rownames(rf_imp)],
          pch = 21, col = "black")
   outliers = glmMeans >= 0.1
   text(names(rfMeans[outliers]), pos=4, cex = 0.7,
        x=rfMeans[outliers], y = glmMeans[outliers]+0.005)
   abline(lm(glmMeans~rfMeans, na.action="na.exclude"))
   if(tissue == tissues[length(tissues)]){
      title(xlab = "RF permutation importance", xpd = NA, 
            mgp = c(2,1,0), cex.lab = 1.5)
   }
   if(tissue == tissues[ceiling(length(tissues)/2)]){
      title(ylab = "abs(glm coefficient)", xpd = NA, 
            mgp = c(2,1,0), cex.lab = 1.5)
   }
   title(ylab = tissue, mgp = c(3.2,1,0), xpd = NA, cex.lab = 1.5)
   # same plot, but outliers removed
   toExclude = c("repeatMasker", "GTEx_eqtl", "Trf")
   rf_imp = rf_imp[!rownames(rf_imp) %in% toExclude,]
   glm_imp = glm_imp[!rownames(glm_imp) %in% toExclude,]
   plot(rf_imp, glm_imp, las = 1, pch = 1, cex = 0.5,
        col = NA,
        xlab = "", ylab = "",
         mgp = c(4,0.4,0), tcl=-0.3)
   rfMeans = rowMeans(rf_imp, na.rm = T)
   glmMeans= rowMeans(glm_imp, na.rm = T)
   rfSD = apply(rf_imp,1,sd, na.rm = T)
   glmSD = apply(glm_imp,1,sd, na.rm = T)
   arrows(x0 = rfMeans-rfSD, x1 = rfMeans+rfSD, 
          y0 = glmMeans, y1 = glmMeans, code=0)
   arrows(x0 = rfMeans, x1 = rfMeans, 
          y0 = glmMeans-glmSD, y1 = glmMeans+glmSD, code=0)
   points(rfMeans, glmMeans,
          bg = predictorCols[rownames(rf_imp)],
          pch = 21, col = "black")
   # outliers = glmMeans >= (mean(glmMeans, na.rm = T) +
   #                            sd(glmMeans, na.rm=T)) |
   #    rfMeans >= (mean(rfMeans, na.rm = T) +
   #                   sd(rfMeans, na.rm=T))
   # text(names(rfMeans[outliers]), pos=4, cex = 0.7,
   #      x=rfMeans[outliers], y = glmMeans[outliers]+0.001,
   #      col = predictorCols[outliers])
   abline(lm(glmMeans~rfMeans, na.action="na.exclude"))
   if(tissue == tissues[length(tissues)]){
      title(xlab = "RF permutation importance", xpd = NA, 
            mgp = c(2,1,0), cex.lab = 1.5)
   }
   if(tissue == tissues[ceiling(length(tissues)/2)]){
      title(ylab = "abs(glm coefficient)", xpd = NA, 
            mgp = c(2,1,0), cex.lab = 1.5)
   }
   # p-values
   rf_pval = -log10(rf_pvals[,tissue])
   glm_pval = -log10(glm_pvals[[tissue]])
   glm_pvalSD = apply(glm_pval,1,sd)
   plot(rowMeans(glm_pval), rf_pval, las = 1, 
        ylab = "-log10(RF p-value)",
        bg = predictorCols, pch = 21, col = "black",
        mgp = c(4,0.4,0), tcl=-0.3)
   segments(y0=rf_pval, x0 = rowMeans(glm_pval)-glm_pvalSD, 
            y1 = rf_pval, x1 = rowMeans(glm_pval)+glm_pvalSD)
   abline(h=-log10(0.05), v = -log10(0.05))
   # disagree = (rf_pval >= -log10(0.05) | 
   #                rowMeans(glm_pval) >= -log10(0.05))  & 
   #                (rf_pval) < 3.9 
   # disagree[is.na(disagree)] = F
   # text(names(rf_pval[disagree]), pos=4, cex = 0.7,
   #      y=rf_pval[disagree]-0.05, x = rowMeans(glm_pval)[disagree],
   #      col = predictorCols[disagree])
   if(tissue == tissues[length(tissues)]){
      title(xlab = "-log10(GLM p-value)", xpd = NA, 
            mgp = c(2,1,0), cex.lab = 1.5)
   }
   if(tissue == tissues[ceiling(length(tissues)/2)]){
      title(ylab = "-log10(RF p-value)", xpd = NA, 
            mgp = c(2,1,0), cex.lab = 1.5)
   }
})


# collection of plot, comparison of tissues #####
png("fig/performance/compareRFperformanceTissues.png", 
    width=1200, height=400, pointsize=20)
par(mfrow = c(1,3), mar = c(4,4,0.1,0.1), mgp = c(2.5,1,0))
plot(NA, xlim = c(0,1), ylim = c(0,1), 
         xlab = "False positive rate", 
         ylab = "True positive rate", las = 1)
abline(0,1, col = "dark grey", lty = 2)
plotDump = sapply(tissues, function(tissue){
   lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]], 
         ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]], 
         col = tissueCols[tissue], lwd = 2)
})
legend("bottomright", col=tissueCols, legend=names(tissueCols), lwd = 1)
plot(NA, xlim = c(0,1), ylim = c(0,1), 
     xlab = "Recall", 
     ylab = "Precision", las = 1)
plotDump = sapply(tissues, function(tissue){
   lines(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]], 
         ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], 
         col = tissueCols[tissue], lwd = 2)
})
AUCs = sapply(tissues, function(tissue){
   ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
})
AUCsPerChr = sapply(tissues, function(tissue){
   sapply(ROC_PR_RF_perChr[[tissue]], function(x){
      x$auc
   })
})
AUCsSD = apply(AUCsPerChr,2,sd)
nMuts = sapply(dataInfos, function(x){x$nMuts})

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
text(x=nMuts, y = colMeans(AUCsPerChr), labels=tissues, pos=4)
dev.off()


# coeff barplot
impsMelted = do.call(rbind,lapply(names(rf_imps), function(tissue){
   x = as.data.frame(rf_imps[[tissue]])
   long <- reshape(x, direction = "long",
                   ids = row.names(x), idvar = "predictors", 
                   times = names(x), timevar = "chr",
                   varying = names(x), v.names = "importance")
   long$tissue = tissue
   return(long)
}))
ggplot(impsMelted, aes(x = predictors, y = importance, col = tissue)) +
   geom_boxplot(show.legend=F)  + 
   coord_flip() +
   facet_wrap(~ tissue, ncol=7,scales="free_x") +
   theme_minimal() #+
ggsave("fig/RFimportances_barplot.png")
# coeff heatmap
meanImps = do.call(rbind,lapply(tissues, function(tissue){
   imp = rowMeans(rf_imps[[tissue]])
   p = rf_pvals[,tissue]
   d = data.frame(importance = imp, pval = p, tissue)
   d$predictor = rownames(d)
   return(d)
}))
meanImps$predictor = factor(meanImps$predictor,levels=rev(predictorOrder))
ggplot(meanImps, aes(x = tissue, y = predictor)) + 
   geom_raster(aes(fill=importance)) +
   scale_fill_gradient(low="grey90", high="red",na.value="grey") +
   geom_text(aes(label = ifelse(pval<=0.01, 
                                yes = "*", no = NA)),
             color = "black",  nudge_y = -0.45)
ggsave("fig/RFimportances_heatmap.png")
meanImpsScaled = do.call(rbind,lapply(tissues, function(tissue){
   imp = rowMeans(rf_imps[[tissue]])
   imp = (imp-min(imp, na.rm = T))/(max(imp, na.rm = T)-min(imp, na.rm = T))
   p = rf_pvals[,tissue]
   d = data.frame(importance = imp, pval = p, tissue)
   d$predictor = rownames(d)
   return(d)
}))
ggplot(meanImpsScaled, aes(x = tissue, y = predictor)) + 
   geom_raster(aes(fill=importance)) +
   scale_fill_gradient(low="grey90", high="red",na.value="grey") +
   geom_text(aes(label = ifelse(pval<=0.01, 
                                yes = "*", no = NA)),
             color = "black",  nudge_y = -0.45)
ggsave("fig/RFimportances_heatmapScaled.png")
#####



