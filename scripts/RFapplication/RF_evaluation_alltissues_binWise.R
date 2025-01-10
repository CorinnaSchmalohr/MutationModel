library(ggplot2)
library(data.table) 
library(ROCR)
library(ranger) 
library(berryFunctions) # for smallPlot
library(corrplot)
tissues = c("luad", "skin", "colon", "ovary",
            "kidney", "prostate", "breast")
tissueCols = setNames(rainbow(length(tissues)), tissues)
chrCols = setNames(rainbow(22), paste0("chr", 1:22))
dir.create("fig/performance", showWarnings=F)
predPerTissue = lapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/RFmodel/predictions_withBinWise.RData"))
   return(predictions)
})
names(predPerTissue) = tissues



# compute ROC, PR, and AUROC for each chromosome #####
png(paste0("fig/performance/ROCandPRoverChrs_allTissues.png"), 
    height=2000, width=2000,pointsize=40)
par(mfrow = c(4,4), mar = c(3.2,4,1.8,0))
ROC_PR_RF_perChr = lapply(names(predPerTissue), function(tissue){ #iterate through tissues
   print(tissue)
   pred = predPerTissue[[tissue]]
   # get performances
   temp = lapply(pred, function(x){ #iterate through chromosomes
      temp = prediction(x$pred_1, x$labels)
      roc = performance(temp, "tpr", "fpr")
      pr = performance(temp,"prec", "rec")
      auc = performance(temp,"auc")
      return(list(roc = roc, pr = pr, auc = auc))
   })
   # plot ROCs
   colors = rainbow(length(temp))
   lwd = 2
   plotExpression = expression({
      sapply(seq(length(temp)), function(i){
         if(i==1){
            plot(temp[[i]]$roc,  main = tissue,
                 col = colors[i], yaxis.las = 1, lwd = lwd,
                 mgp = c(2.3,1,0))
         } else{
            plot(temp[[i]]$roc,  col = colors[i], add = T, lwd = lwd,
                 mgp = c(2.3,1,0))
         }
      })
      abline(0,1, col = "grey", lty = 2, lwd = lwd)
   })
   png(paste0("fig/performance/ROCoverChrs_", tissue, ".png"),
       height=1200, width=1200,pointsize=30)
   eval(plotExpression)
   legend("bottomright", legend=names(temp), 
          col = colors, lty=1, ncol=3, cex=0.8, lwd = lwd)
   dev.off()
   lwd = 1
   eval(plotExpression)
   # plot PRs
   lwd = 2
   plotExpression = expression({
      sapply(seq(length(temp)), function(i){
         if(i==1){
            plot(temp[[i]]$pr,  main = tissue, col = colors[i],
                 ylim = c(0,1), yaxis.las = 1, lwd = lwd, mgp = c(2.3,1,0))
         } else{
            plot(temp[[i]]$pr, col = colors[i], add = T, lwd = lwd,
                 mgp = c(2.3,1,0))
         }
      })
   })
   png(paste0("fig/performance/PRoverChrs_", tissue, ".png"), 
       height=1200, width=1200,pointsize=30)
   eval(plotExpression)
   legend("bottomright", legend=names(temp), 
          col = colors, lty=1, ncol=3, lwd = lwd)
   dev.off()
   lwd = 1
   eval(plotExpression)
   abline(v=1.1,h = -0.39,xpd = NA)
   return(temp)
})
plot.new()
lwd = 1
legend(x = 0.5, y = 0.9, legend=names(ROC_PR_RF_perChr[[1]]), 
       col = rainbow(length(ROC_PR_RF_perChr[[1]])), 
       lty=1, ncol=3, lwd = lwd, xpd = NA)
dev.off()
names(ROC_PR_RF_perChr) = names(predPerTissue)
save(ROC_PR_RF_perChr, file = "data/rdata/ROC_PR_RF_perChr.RData")
#####



# compute ROC, PR, and AUROC for all chromosomes concatenated #####
png(paste0("fig/performance/ROCandPR_rf_vs_glm_allTissues.png"),
    height=2000, width=2000,pointsize=30)
par(mfrow = c(4,4), mar = c(3.2,4,1.8,0))
ROC_PR_RF_concat = lapply(names(predPerTissue), function(tissue){
   # concatenate predictions
   predConcat = do.call(rbind,predPerTissue[[tissue]])
   # calculate perfromance measures
   ROC_rf = performance(prediction(predConcat$pred_1, 
                                   predConcat$labels), 
                        "tpr", "fpr")
   AUC_rf = performance(prediction(predConcat$pred_1, 
                                   predConcat$labels), 
                        "auc")
   PR_rf = performance(prediction(predConcat$pred_1, 
                                  predConcat$labels), 
                       "prec", "rec")
   # load results for glm (created by nina)
   ROC_glm = load(paste0("data/nina/rdata_withBinwise/", tissue, 
                         "/data_glm/performROC_concat_glm.RData"))
   PR_glm = load(paste0("data/nina/rdata_withBinwise/", tissue, 
                        "/data_glm/performPR_concat_glm.RData"))
   # compare results from RF with results from glm
   # ROC
   plotExpression = expression({
      plot(get(ROC_glm), lwd = 2, main = tissue, mgp = c(2.3,1,0))
      plot(ROC_rf, add = T, col = "red", lty = 2, lwd = 2, mgp = c(2.3,1,0))
      abline(0,1, col = "grey")})
   png(paste0("fig/performance/ROC_rf_vs_glm_", tissue, ".png"),
       height=1200, width=1200,pointsize=30)
   eval(plotExpression)
   legend("bottomright", legend=c("glm", "RF"), col = c("black", "red"),
          lty = 1:2, lwd = 2, xpd=NA, inset=0.1)
   dev.off()
   eval(plotExpression)
   
   #PR
   # plotExpression = expression({
   #    plot(get(PR_glm), lwd = 2, main = tissue, ylim = c(0,1), mgp = c(2.3,1,0))
   #    plot(PR_rf, add = T, col = "red", lty = 2, lwd = 2, mgp = c(2.3,1,0))})
   plotExpression = expression({
      plot(get( PR_glm), lwd = 2, main = tissue, ylim = c(0,1),
           mgp = c(2.3,1,0),yaxis.las = 1)
      plot(PR_rf, add = T, col = "red", lty = 2, lwd = 2,
           mgp = c(2.3,1,0))
      lim = par("plt")
      inset = 0.1
      smallPlot(expr={
         plot(x = get(PR_glm)@x.values[[1]],
              y = get(PR_glm)@y.values[[1]],
              xlab = "", ylab = "",mgp = c(0,0.3,0), xlim = c(0,0.1),
              ylim = c(0.5,0.9), lwd = 2,
              type = "l", lab = c(1,2,7), las = 1, tcl = -0.2)
         lines(x = PR_rf@x.values[[1]], 
               y = PR_rf@y.values[[1]],col = "red", lty = 2, lwd = 2)
      }, 
      x1 = lim[1]+inset, x2 = lim[2]*0.7,
      y1 = lim[3]+inset, y2 = lim[4]*0.6, xpd = F,
      mar = c(0,0,0,0), border = "transparent")
   })
   png(paste0("fig/performance/PR_rf_vs_glm_", tissue, ".png"),
       height=1200, width=1200,pointsize=30)
   eval(plotExpression)
   legend("bottomright", legend=c("glm", "RF"), col = c("black", "red"),
          lty = 1:2, lwd = 2)
   dev.off()
   eval(plotExpression)
   return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
})
plot.new()
lwd = 1
legend(x = 0.5, y = 0.9, legend=c("glm", "RF"), col = c("black", "red"),
       lty = 1:2, lwd = 2)
dev.off()
names(ROC_PR_RF_concat) = names(predPerTissue)
save(ROC_PR_RF_concat, file = "data/rdata/ROC_PR_RF_concat.RData")
# plot ROCs of all tissues
png("fig/performance/ROC_rf_alltissues.png",
    height=1200, width=1200,pointsize=30)
plot(NA, xlim = c(0,1), ylim = c(0,1), las = 1, 
     xlab = "False positive rate", ylab = "True positive rate")
abline(0,1, lwd =2, lty = 2, col = "grey")
dumpVar = sapply(names(ROC_PR_RF_concat), function(tissue){
   plot(ROC_PR_RF_concat[[tissue]]$roc, add = T, 
        col = tissueCols[tissue], lwd = 2)
})
legend("bottomright", col = tissueCols, legend=names(tissueCols), 
       lwd = 2)
dev.off()
# plot PRs of all tissues
png("fig/performance/PR_rf_alltissues.png",
    height=1200, width=1200,pointsize=30)
plot(NA, xlim = c(0,1), ylim = c(0,1), las = 1, 
     xlab = "Recall", ylab = "Precision")
colors = rainbow(length(ROC_PR_RF_concat))
dumpVar = sapply(names(ROC_PR_RF_concat), function(tissue){
   plot(ROC_PR_RF_concat[[tissue]]$pr, add = T, 
        col = tissueCols[tissue], lwd = 2)
   return(NULL)
})
legend("bottomright", col = tissueCols, legend=names(tissueCols), 
       lwd = 2)
dev.off()
#####


# correlate AUC with mutation rate and nPositions ########
load("data/rdata/ROC_PR_RF_perChr.RData")
png("fig/performance/AUCvsNpos_alltissues.png",height = 500, width = 1000)
par(mfrow = c(2,4))
dumpVar = sapply(tissues, function(tissue){
   AUCs =sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc@y.values[[1]]})
   load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
   nMuts = table(removed$chr)/1000
   plot(as.integer(nMuts), AUCs, las = 1, 
        ylab = "AUC", xlab = "k Positions", main = tissue)
   abline(lm(AUCs~nMuts))
   legend("topright", bty = "n",
          legend=paste0("p=", 
                        format(summary(lm(AUCs~nMuts))$coefficients[2,4],
                               digits=2)))
})
dev.off()

png("fig/performance/AUCvsTPpercentage_alltissues.png")
par(mfrow = c(2,4))
dumpVar = sapply(tissues, function(tissue){
   AUCs =sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc@y.values[[1]]})
   load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
   nMuts = sapply(split(as.numeric(data$mutated)-1, removed$chr), mean)
   plot(nMuts, AUCs, las = 1, 
        ylab = "AUC", xlab = "% TP Positions", main = tissue)
   abline(lm(AUCs~nMuts))
   legend("topright", bty = "n",
          legend=paste0("p=", 
                        format(summary(lm(AUCs~nMuts))$coefficients[2,4],
                               digits=2)))
})
dev.off()

load("data/rdata/ROC_PR_RF_concat.RData")
AUCs = sapply(ROC_PR_RF_concat, function(x)x$auc@y.values[[1]])
nMuts = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
   return(nrow(data))
})
png("fig/performance/AUCsVsNpos_concat.png", width=800, height=800,pointsize=25)
plot(nMuts, AUCs, las = 1, pch = 19)
abline(lm(AUCs ~ nMuts))
dev.off()
######

# compare AUCs to glm #####
png("fig/performance/AUC_rf_vs_glm_alltissues.png",
    height=1200, width=1200,pointsize=30)
AUCsBoth = sapply(tissues, function(tissue){
   rf = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
   glm = load(paste0("data/nina/rdata_withBinwise/", tissue, 
                     "/data_glm/performAUC_concat_glm.RData"))
   return(c(rf = rf, glm = get(glm)@y.values[[1]]))
})
barplot(AUCsBoth, beside=T, col = rep(tissueCols, each = 2),
        density=c(NA,20), las = 1, ylab = "AUC")
legend("topright", legend=c("RF", "GLM"), 
       fill = "grey30", density=c(NA,20), inset=0)
dev.off()

AUCsPerChrBoth = lapply(tissues, function(tissue){
   rf = sapply(ROC_PR_RF_perChr[[tissue]], function(x){
      x$auc@y.values[[1]]
   })
   glm = get(load(paste0("data/nina/rdata_withBinwise/", tissue,
                         "/data_glm/performROC_AUC_glm_Loco-CV.RData")))
   return(cbind(rf = rf, glm = glm))
})

names(AUCsPerChrBoth) = tissues
temp = melt(AUCsPerChrBoth)
colnames(temp) = c("chr", "method", "auc", "tissue")
g = ggplot(temp, aes(tissue, auc, fill = method), ylim = c(0,1)) +
   geom_boxplot() +
   labs(x = "", y = "AUC") +
   theme_minimal() +
   theme(text=element_text(size=20),
         plot.title = element_text(hjust=0.5),
         legend.title=element_blank())
ggsave(g,filename="fig/performance/AUC_rf_vs_glm_Boxplot_alltissues.png")
######


# TPs higher prediction than TNs? #####
dumpVar = lapply(names(predPerTissue), function(tissue){
   temp = do.call(rbind,predPerTissue[[tissue]])
   g = ggplot(temp, aes(labels, pred_1), ylim = c(0,1)) +
      geom_violin(stat="ydensity") +
      labs(x = "True value", y = "Predicted value") +
      ggtitle(tissue) +
      scale_y_continuous(limits = c(0, 1))+
      theme_minimal() +
      theme(text=element_text(size=20),
            plot.title = element_text(hjust=0.5))
   ggsave(g,file = paste0("fig/performance/predsTPvsTN_", tissue, ".png"))
   
   temp = melt(predPerTissue[[tissue]])
   temp = temp[temp$variable == "pred_1",]
   g = ggplot(temp, aes(x = L1, y = value, fill = labels), 
              ylim = c(0,1),) +
      geom_violin(position=position_dodge(0.8)) +
      labs(x = "Chromosome", y = "Predicted value") +
      ggtitle(tissue) +
      scale_y_continuous(limits = c(0, 1)) +
      theme_minimal() +
      theme(text=element_text(size=20),
            axis.text.x=element_text(size=rel(1), angle=35,hjust=1),
            plot.title = element_text(hjust=0.5))
   ggsave(g,file = paste0("fig/performance/predsTPvsTNperChr_", tissue, ".png"),
          width=15, height=7)
})

#####

# corrplots #####
dumpVar = sapply(tissues, function(tissue){
   print(tissue)
   
   load(paste0("data/rdata/", tissue, 
               "/completeData_withBinwise.RData"))
   data$context = NULL
   data$pentamer = NULL
   data$trimer = NULL
   data$septamer = NULL
   data$inexon = NULL
   png(paste0("fig/performance/corrplot_binWise_", tissue, ".png"),
       height = 1200, width = 1200, pointsize = 30)
   corrplot(cor(data[,sapply(data, is.numeric)]),
            tl.col = "black", tl.cex = 0.6)
   dev.off()
   return(NULL)
})
#####

# get predictor importances from rf and glm#####
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
rf_imps = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, 
               "/RFmodel/imp_withBinWise.RData"))
   t(sapply(predictorOrder, function(i){
      if(i %in% rownames(imp)){(imp[i,])} else{rep(NA,ncol(imp))}
   }))
}, simplify=F)
glm_imps = sapply(tissues, function(tissue){
   load(paste0("data/nina/rdata_withBinwise/", tissue, 
               "/data_glm/coefs_allChr_glm.RData"))
   coefs_allChr = as.matrix(coefs_allChr)
   temp = t(sapply(predictorOrder, function(i){
      if(i %in% rownames(coefs_allChr)){(coefs_allChr[i,])} else{rep(NA,ncol(coefs_allChr))}
   }))
   return(temp)
}, simplify=F)
#####


# get predictor p-values #####
rf_pvals = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, 
               "/RFmodel/imp_pvals_withBinWise.RData"))
   sapply(predictorOrder, function(i){
      if(i %in% rownames(pvals)){(pvals[i,2])} else{NA}
   })
})
glm_pvals = sapply(tissues, function(tissue){
   load(paste0("data/nina/rdata_withBinwise/", tissue,
               "/data_glm/pValues_features_allChr_glm.RData"))
   pValue_features_allChr = as.matrix(pValue_features_allChr)
   temp = t(sapply(predictorOrder, function(i){
      if(i %in% rownames(pValue_features_allChr)){
         (pValue_features_allChr[i,])
      } else{
         rep(NA,ncol(pValue_features_allChr))
      }
   }))
   temp[temp == 0] = 2.2e-16
   return(temp)
}, simplify=F)
#####


# plotting of importances and pvalues #####
# barplot of importances, each tissue one column
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
   theme_minimal() +
   theme(text=element_text(size=20))   
ggsave("fig/RFpermutationImportances_binWise_alltissues.png", 
       height=12, width=12)
# image plot of mean (over chrs) importance values
meanImps = as.data.frame(sapply(tissues, function(tissue){
   rowMeans(rf_imps[[tissue]])
}))
meanImpsMelted = reshape(meanImps,  direction = "long",
                         ids = row.names(meanImps), idvar = "predictors", 
                         times = names(meanImps), timevar = "tissue",
                         varying = names(meanImps), v.names = "importance")
ggplot(meanImpsMelted, aes(x = tissue, y = predictors)) + 
   geom_raster(aes(fill=importance)) +
   scale_fill_gradient(low="grey90", high="red",na.value="grey")+
   theme(text=element_text(size=20))   
ggsave("fig/permutationImportancesMean_binWise_alltissues.png",
       height=14, width=11)
# image plot of importance values (like above), but scaled across tissues
meanImpsScaled = as.data.frame(apply(meanImps,2,function(x){
   (x-min(x,na.rm = T)) / (max(x,na.rm = T)-min(x,na.rm = T))
}))
meanImpsScaledMelted = reshape(meanImpsScaled,  direction = "long",
                               ids = row.names(meanImpsScaled), idvar = "predictors", 
                               times = names(meanImpsScaled), timevar = "tissue",
                               varying = names(meanImpsScaled), v.names = "importance")
ggplot(meanImpsScaledMelted, aes(x = tissue, y = predictors)) + 
   geom_raster(aes(fill=importance)) +
   scale_fill_gradient(low="grey90", high="red",na.value="grey")+
   theme(text=element_text(size=20))   
ggsave("fig/permutationImportancesMeanScaled_binWise_alltissues.png",
       height=14, width=11)
# scatter plots comparing tissues in  panels. 
# upper triangle chr-wise values, lower triangle mean with sd
impsReformatted = data.frame(sapply(names(rf_imps), function(tissue){
   x = as.data.frame(rf_imps[[tissue]])
   temp = reshape(x,  direction = "long",
                  ids = row.names(x), idvar = "predictors",
                  times = names(x), timevar="chr", 
                  varying = list(names(x)), v.names = "importance")
   if(tissue == names(rf_imps)[1]){
      temp = temp[,c("predictors", "chr", "importance")]
      colnames(temp) = c("predictors", "chr", tissue)
      return(temp)
   } else
      return(temp[,"importance"])
}, simplify=F))
colnames(impsReformatted)[1:3] = do.call(rbind,
                                         strsplit(colnames(impsReformatted)[1:3], 
                                                  split=".",fixed = T))[,2]
predictorCols = setNames(rainbow(length(predictorOrder)), predictorOrder)
png("fig/compareRFpredictorValuesBetweenTissue.png",
    height=1000, width=1000, pointsize=25)
pairs(impsReformatted[,-(1:2)],gap=0.5, oma = c(4,4,2,2),
      panel = function(x,y){
         points(x,y,pch = ".", cex = 2,
                col = predictorCols[impsReformatted$predictor])
         abline(lm(y~x))},
      lower.panel=function(x,y){
         tempX = split(x,impsReformatted$predictor)
         tempY = split(y,impsReformatted$predictor)
         Xmeans = sapply(tempX, mean, na.rm = T)
         Ymeans = sapply(tempY, mean, na.rm = T)
         Xsd = sapply(tempX, sd, na.rm = T)
         Ysd = sapply(tempY, sd, na.rm = T)
         arrows(x0=Xmeans, x1 = Xmeans,y0 = Ymeans-Ysd,
                y1 = Ymeans+Ysd, code = 0, lwd = 1.5)
         arrows(x0=Xmeans-Xsd, x1 = Xmeans+Xsd,y0 = Ymeans, 
                y1 = Ymeans, code = 0, lwd = 1.5)
         points(Xmeans, Ymeans, cex = 1.5,
                col = predictorCols[names(tempX)], pch = 1)
      })
mtext(side=1, "permutation importance", line = 3.4)
mtext(side=2, "permutation importance", line = 3)
dev.off()
# compare rf and glm coefficients
png("fig/compareGLMandRFpredictorValues.png", width=800, height=800, pointsize=25)
par(mfrow = c(2,4), mar = c(4,4,2,0))
dumpVar = sapply(tissues, function(tissue){
   rf_imp = rf_imps[[tissue]]
   glm_imp = abs(glm_imps[[tissue]])
   plot(rf_imp, glm_imp, las = 1, pch = 1, cex = 0.5,
        col = predictorCols[rownames(rf_imp)],
        main = tissue, xlab = "RF permutation importance", 
        ylab = "glm coefficient", mgp = c(2.5,1,0))
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
   text(names(rfMeans[glmMeans>=0.1]), pos=4, cex = 0.7,
        x=rfMeans[glmMeans>=0.1], y = glmMeans[glmMeans>=0.1])
   abline(lm(glmMeans~rfMeans, na.action="na.exclude"))
})
plot.new()
legend("topleft", col=predictorCols, cex = 0.75, inset=-0.1,
       legend=names(predictorCols), pch = 19, ncol = 3, xpd = NA)
dev.off()
# compare rf and glm coefficients, without repeatMasker, GTEx_eqtl, and Trf
png("fig/compareGLMandRFpredictorValuesWoutOutliers.png", 
    width=800, height=800, pointsize=25)
toExclude = c("repeatMasker", "GTEx_eqtl", "Trf")
par(mfrow = c(2,4), mar = c(4,4,2,0))
dumpVar = sapply(tissues, function(tissue){
   rf_imp = rf_imps[[tissue]]
   glm_imp = abs(glm_imps[[tissue]])
   rf_imp = rf_imp[!rownames(rf_imp) %in% toExclude,]
   glm_imp = glm_imp[!rownames(glm_imp) %in% toExclude,]
   plot(rf_imp, glm_imp, las = 1, pch = 1, cex = 0.5,
        col = predictorCols[rownames(rf_imp)],
        main = tissue, xlab = "RF permutation importance", 
        ylab = "glm coefficient", mgp = c(2.5,1,0))
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
   text(names(rfMeans[glmMeans>=0.1]), pos=4, cex = 0.7,
        x=rfMeans[glmMeans>=0.1], y = glmMeans[glmMeans>=0.1])
   abline(lm(glmMeans~rfMeans, na.action="na.exclude"))
})
plot.new()
legend("topleft", col=predictorCols, cex = 0.75, inset=-0.1,
       legend=names(predictorCols), pch = 19, ncol = 3, xpd = NA)
dev.off()
#####




