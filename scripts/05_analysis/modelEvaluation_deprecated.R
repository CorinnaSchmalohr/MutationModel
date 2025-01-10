library(ggplot2)
library(data.table) 
library(ROCR)
library(ranger) 
library(berryFunctions) # for smallPlot
library(corrplot)
library(sinaplot)
library(RColorBrewer)
library(plotrix) # for point labels
library(readxl)
source("lib/general_function.R")
source("scripts/05_analysis/00_NamesAndColors.R")

dir.create("fig/modelEvaluation", showWarnings=F)

# get mutation rate and nPositions, percTP, correlation etc. (dataInfos) #####
dataInfos = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue, "_Muts_mapped_processed.RData"))
  chrCount = sapply(table(datchroms), as.integer)
  chrPerc = sapply(split(dat$mutated, datchroms),function(x){mean(x==1)})
  cors = cor(dat[sapply(dat, is.numeric)], use = "pair")
  return(list(nMuts = nrow(dat), 
              percTP = mean(dat$mutated == 1), 
              nMutsPerChr = chrCount,
              percTPperChr = chrPerc, 
              cors = cors))
}, simplify = F)
save(dataInfos, file = "data/processedData/dataInfos.RData")
#####

# load predictions for each model and each tissue #####
# random forest
predPerTissueRF = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue,
              "_predictions.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueRF, file = "data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")
# load sig. GLM predictions for each tissue
predPerTissueGLMsig = sapply(tissues, function(tissue){
  load(file = paste0("data/Modeling/exomeTrainData/GLM/", tissue, 
                     "_predictions_sig.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueGLMsig, file = "data/Modeling/exomeTrainData/GLM/predPerTissueGLM_sig.RData")
predPerTissueLasso = sapply(tissues, function(tissue){
  load(file = paste0("data/Modeling/exomeTrainData/Lasso/", 
                     tissue, "_predictions_sig.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueGLMsig, file = "data/Modeling/exomeTrainData/Lasso/predPerTissueLasso.RData")
#####


# compute ROC, PR, and AUROC for each chromosome for RF #####
# RF
ROC_PR_RF_perChr = sapply(names(predPerTissueRF), function(tissue){ #iterate through tissues
  pred = predPerTissueRF[[tissue]]
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
save(ROC_PR_RF_perChr, file = "data/Modeling/exomeTrainData/RF/ROC_PR_RF_perChr.RData")
# glm
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
save(ROC_PR_glm_perChr_sig, file = "data/Modeling/exomeTrainData/GLM/ROC_PR_glm_perChr_sig.RData")
# lasso
ROC_PR_lasso_perChr = sapply(names(predPerTissueLasso), function(tissue){
  pred = predPerTissueLasso[[tissue]]
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
save(ROC_PR_lasso_perChr, file = "data/Modeling/exomeTrainData/Lasso/ROC_PR_glm_perChr_sig.RData")
#####



# compute ROC, PR, and AUROC for all chromosomes concatenated #####
ROC_PR_RF_concat = sapply(names(predPerTissueRF), function(tissue){
  predConcat = do.call(rbind,predPerTissueRF[[tissue]])
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
save(ROC_PR_RF_concat, file = "data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
#  glm
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
save(ROC_PR_glm_concat_sig, file = "data/Modeling/exomeTrainData/GLM/ROC_PR_glm_concat_sig.RData")
# lasso
ROC_PR_lasso_concat = sapply(names(predPerTissueLasso), function(tissue){
  predConcat = do.call(rbind,predPerTissueLasso[[tissue]])
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
save(ROC_PR_lasso_concat, file = "data/Modeling/exomeTrainData/Lasso/ROC_PR_lasso_concat.RData")
#####


# get predictor order #####
tab = read_xlsx("data/rawdata/dataMappingAlltissues.xlsx", 
                sheet="allTissues", col_names=T)
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)
tab$NA. = NULL
tab[tab == "NA"] = NA
# for predictors where we want multiple ranges, expand table
tab = apply(tab,1,function(x){
  if(is.na(x["range"])){
    return(x)
  } else{
    rangeWindow = strsplit(x["range"],";")[[1]]
    rangeInds = which(names(ranges) == rangeWindow[2]):which(names(ranges) == rangeWindow[1])
    subRanges = ranges[rangeInds]
    t(sapply(names(subRanges), function(r,y){
      y["range"] = subRanges[r]
      if(length(subRanges)>1){
        y["Name"] = paste0(y["Name"]," ",r)
        y["abbreviation"] = paste0(y["abbreviation"],"_",r)
      }
      return(y)
    },y=x))
  }
})
tab = do.call(rbind, tab)
tab = as.data.frame(tab)
predictorOrder = tab[,1:3] # Group, Name, abbreviation
save(predictorOrder, file = "data/processedData/predictorOrder.RData")
#####


# get predictor importances #####
# rf
rf_gini = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue,
              "_importances_gini.RData"))
  names(imp) = names(chrCols)
  imp = as.data.frame(imp)
  temp = imp[predictorOrder$abbreviation,]
  rownames(temp) = predictorOrder$abbreviation
  return(temp)
}, simplify=F)
save(rf_gini, file = "data/Modeling/exomeTrainData/RF/RF_imps.RData")
# glm predictor importances
glm_imps = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/",
                tissue, "_", cr, ".RData"))
    logR$coefficients[predictorOrder$abbreviation]
  })
  return(temp)
}, simplify=F)
# glm pvals
glm_pvals = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/",
                tissue, "_", cr, ".RData"))
    pvals = coef(summary(logR))[predictorOrder$abbreviation,4]
  })
  temp[temp == 0] = 2e-16
  return(temp)
}, simplify=F)
save(glm_imps, glm_pvals, 
     file = "data/Modeling/exomeTrainData/GLM/GLM_impsAndPvals.RData")

# glm significant coefficients
glm_imps_sig = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/", 
                tissue, "_", cr, "_sig.RData"))
    logR$coefficients[predictorOrder$abbreviation]
  })
  return(temp)
}, simplify=F)
# glm significant pvals
glm_pvals_sig = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/GLM/",
                tissue, "_", cr, "_sig.RData"))
    pvals = coef(summary(logR))[predictorOrder$abbreviation,4]
  })
  temp[temp == 0] = 2e-16
  return(temp)
}, simplify=F)
save(glm_imps_sig, glm_pvals_sig, 
     file = "data/Modeling/exomeTrainData/GLM/GLM_impsAndPvals_sig.RData")
# lasso
lasso_stability = apply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/Lasso/",
                tissue, "_", cr, ".RData")) # sp and stab
    sp$x[,stab$lpos]
  })
  temp = as.data.frame(temp)
  temp = temp[predictorOrder$abbreviation,]
  rownames(temp) = predictorOrder$abbreviation
  return(temp)
})
lasso_imp = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/exomeTrainData/Lasso/",
                tissue, "_", cr, "_sig.RData")) # logR, sigFeatures
    logR$coefficients[,4]
  })
  temp = temp[predictorOrder$abbreviation,]
  rownames(temp) = predictorOrder$abbreviation
  return(temp)
})
save(lasso_stability, lasso_imp, 
     file = "data/Modeling/exomeTrainData/Lasso/Lasso_impAndStab.RData")
#####

# load data #####
#load("data/rdata_old/RFmodel/ROC_PR_RF_concat.RData")
#ROC_PR_RF_concat_old = ROC_PR_RF_concat
#load("data/rdata_old//dataInfos.RData")
#dataInfos_old = dataInfos
#load("data/rdata_old/RFmodel/RF_impsAndPvals.RData")
#rf_gini_old = rf_gini
#rm(rf_imps)

load("data/rdata/dataInfos.RData")

load("data/rdata/RFmodel/predPerTissueRF.RData")
load("data/rdata/GLMmodel/predPerTissueGLM_sig.RData")
load("data/rdata/RFmodel/ROC_PR_RF_perChr.RData")
load("data/rdata/RFmodel/ROC_PR_RF_concat.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_perChr_sig.RData")
load("data/rdata/GLMmodel/ROC_PR_glm_concat_sig.RData")
load("data/rdata/RFmodel/RF_imps.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals.RData")
load("data/rdata/GLMmodel/GLM_impsAndPvals_sig.RData")
load("data/rdata/nPerTissue.RData")
predictorOrder = rownames(rf_gini[[1]])
#####


#plot n samples per tissue #####
library(sinaplot)
png("fig/nMutsPerSample.png", width=1500, height=800, pointsize=25)
sinaplot(nPerTissue, xaxt = "n", log = T)
sinaplot(nPerTissue, log=T, las = 1, scale=T, maxwidth=1,
         col = tissueCols, xaxt = "n", ylab = "n mutations per sample")
boxplot(sapply(nPerTissue, as.numeric), add = T, 
        xaxt = "n", yaxt = "n", outline=F, boxwex = 0.2, col = NULL)
axis(1, at=1:10, labels=t2T)
dev.off()
#####



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
    title(ylab = "Mutation type", mgp = c(4,1,0), xpd = NA, cex.lab = 1.1)
  }
  abline(h=seq(16,80,length.out = 5), lty = 2, lwd = 1.5)
  abline(h=0, lty = 1, lwd = 2)
  #abline(h=96, lty = 1, lwd = 1.8)
})
title(xlab = "Proportion of variants", outer = T, mgp = c(2,1,0), cex.lab = 1.1)
dev.off()
#####



# collection of plots for every tissue #####
plotDump = sapply(tissues, function(tissue){
  png(paste0("fig/modelEvaluation/perfRFsummary_", tissue, ".png"),
      width=1300, height=1000, pointsize = 25)
  m = rbind(c(1,1,2,2,5,5,5,5), 
            c(1,1,2,2,5,5,5,5), 
            c(3,3,4,4,5,5,5,5), 
            c(3,3,4,4,5,5,5,5), 
            c(6,6,7,7,5,5,5,5),
            c(6,6,7,7,5,5,5,5))
  par(oma = c(0.1,0.1,1.3,0.1), mar = c(4,4,0.5,0))
  layout(m)
  # 1. roc over Chrs 
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2.7,1,0), las = 1,
       xlab = "False positive rate", ylab = "True positive rate")
  abline(0,1, col = "darkgrey", lty = 2)
  dump = sapply(seq_len(length(chrCols)), function(i){
    plot(ROC_PR_RF_perChr[[tissue]][[i]]$roc, 
         col = chrCols[i], add = T)
  })
  plot(ROC_PR_RF_concat[[tissue]]$roc, 
       col = "black", add = T, lwd = 2, main = "ROC")
  # 2. pr over Chrs
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2.7,1,0), las = 1,
       xlab = "Recall", ylab = "Precision")
  dump = sapply(seq_len(length(chrCols)), function(i){
    plot(ROC_PR_RF_perChr[[tissue]][[i]]$pr, 
         col = chrCols[i], add = T)
  })
  plot(ROC_PR_RF_concat[[tissue]]$pr, 
       col = "black", add = T, lwd = 2, main = "PR")
  # 3.  AUC vs mutrate
  AUCs = sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc})
  percMut = dataInfos[[tissue]]$percTPperChr
  plot(percMut, AUCs, las = 1, 
       ylab = "AUROC", xlab = "percent TP", mgp = c(2.7,1,0))
  mod = lm(AUCs~percMut)
  abline(mod)
  p = summary(mod)$coefficients[2,4]
  legend("topright", bty = "n",
         legend=paste0("p=", format(p,digits=2)))
  # 4. AUC vs nPos 
  nMuts = dataInfos[[tissue]]$nMutsPerChr
  plot(nMuts, AUCs, las = 1, mgp = c(2.7,1,0),
       ylab = "AUROC", xlab = "*1000 positions per chr")
  mod = lm(AUCs~nMuts)
  abline(mod)
  p = summary(mod)$coefficients[2,4]
  legend("topright", bty = "n",
         legend=paste0("p=", format(p, digits=2)))
  # 5. legend
  #plot.new()
  #plot.new()
  #legend("center", legend=c(names(chrCols), "concat"), cex = 0.75,
  #        col = c(chrCols, "black"), lty=1, ncol=2, bty = "n")
  # 6. volcano plot/compare permutation imp and gini
  # pvals = rf_pvals[,tissue]
  # pvals = runif(nrow(coeffs))
  #coeffs = rf_imps[[tissue]]
  #gini = rf_gini[[tissue]]
  #coeffMeans = rowMeans(coeffs, na.rm = T)
  #giniMeans = rowMeans(gini, na.rm = T)
  #coeffSD = apply(coeffs,1,sd, na.rm = T)
  # giniSD = apply(gini, 1,sd, na.rm = T)
  #plot(coeffs,gini, col = rainbow(nrow(gini)),
  #     xlab = "permutation importance", 
  #      ylab = "gini corrected")
  #arrows(x0 = coeffMeans-coeffSD, x1 = coeffMeans+coeffSD, 
  #       y0 = giniMeans, y1 = giniMeans, code=0)
  #arrows(x0 = coeffMeans, x1 = coeffMeans, 
  #        y0 = giniMeans-giniSD, y1 = giniMeans+giniSD, code=0)
  # points(coeffMeans, giniMeans,
  #       bg = rainbow(nrow(gini)),
  #      pch = 21, col = "black", cex = 1.2)
  # plot(rowMeans(coeffs), gini,mgp = c(2.7,1,0),
  #      xlab = "importance", ylab = "-log10(p-value)")
  
  
  # 7. plot of importances
  par(mar = c(4,17.5,0,0.7))
  imps = rf_gini[[tissue]]
  imps = imps[apply(imps,1,function(x){!all(is.na(x))}),]
  bxplt=boxplot(t(imps), las = 1, drop =T,
                horizontal = T, cex.axis = 0.75, xaxt = "n",
                xlim = c(2,nrow(imps)-1), yaxt = "n")
  axis(1)
  axis(2, at=seq(1, length(rownames(imps))), labels= p2P[rownames(imps)], las = 2)
  mtext("Gini importance", side = 1, cex= 0.8, line = 2.5)
  segments(y0= 1:nrow(imps), y1 = 1:nrow(imps),
           x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
  
  
  # 8. TP vs TN preds violin plot
  par(mar = c(4,4,0.5,0))
  predsConcat = do.call(rbind,predPerTissueRF[[tissue]])
  sinaplot(pred ~ label, predsConcat, las = 1, 
           xlab = "True value", ylab = "Predicted value", 
           col = c("light grey", "dark grey"), maxwidth = 0.9,mgp = c(2.7,1,0))
  boxplot(pred ~ label, predsConcat, add = T, yaxt = "n",xaxt = "n",
          col = "white", boxwex = c(0.5,0.5), outline = F)
  # 9. mutation type
  load(paste0("data/procData/traindata/traindata_", tissue, ".RData"))
  muts = data$muts
  toExclude = data$pred$ConsensusExcludable == 1 | 
    data$pred$repeatMasker ==1 | 
    data$pred$tandemRepeatFinder == 1
  muts = muts[!toExclude,]
  muts = data$muts[ data$muts$mutated == 1,]
  muttypes = paste0(muts$ref, ">", muts$alt)
  muttypes = table(muttypes)
  mutDat = cbind(muttypes[1:6], -muttypes[12:7])
  labels = paste(names(muttypes)[1:6], 
                 names(muttypes)[12:7], sep=" / ")
  par(mar = c(4,8,0.5,0))
  barplot(mutDat[,1], horiz=T, col = rainbow(6),
          las = 1, xlim = c(-max(mutDat),max(mutDat)),names.arg=labels,
          cex.lab=1.3,cex.names=1.3, xlab = "count", xaxt = "n")
  barplot(mutDat[,2], add = T, horiz = T, col = rainbow(6),
          xaxt = "n", yaxt = "n")
  temp = axTicks(1)
  axis(1, at = temp, labels=abs(temp))
  title(t2T[tissue], outer = T)
  dev.off()
  # separately: corrplot
  png(paste0("fig/modelEvaluation/corrplot_", tissue, ".png"),
      width=1550, height=1300, pointsize = 25)
  cors = dataInfos[[tissue]]$cors
  rownames(cors) = p2P[rownames(cors)]
  colnames(cors) = p2P[colnames(cors)]
  corrplot(corr = cors, tl.cex=0.8, tl.col="black")
  dev.off()
})
#####

# collection of plots for every tissue improved #####
png(paste0("fig/modelEvaluation/RF_ROCsOverChrsAllTissues.png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.7,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  abline(0,1, col = "grey", lty = 2)
  dump = sapply(names(chrCols), function(cr){
    plot(ROC_PR_RF_perChr[[tissue]][[cr]]$roc, 
         col = rgb(0,0,0,0.5), add = T)
  })
  lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]], 
        ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2.5)
  text(x=0.05,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0,0))
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0))
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0))
  }
})
title(ylab = "TPR", mgp = c(2,0.7,0), xpd = NA, outer=T)
title(xlab = "FPR", mgp = c(1.8,0.7,0), xpd = NA, outer=T)
legend("bottomright", lty = 1, lwd = c(1,2.5), bty = "n",
       col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"))
dev.off()


png(paste0("fig/modelEvaluation/RF_PRsOverChrsAllTissues.png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.8,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  dump = sapply(names(chrCols), function(cr){
    plot(ROC_PR_RF_perChr[[tissue]][[cr]]$pr, 
         col = rgb(0,0,0,0.5), add = T)
  })
  plot(ROC_PR_RF_concat[[tissue]]$pr, 
       col = tissueCols[tissue], add = T, lwd = 2.5, main = "")
  text(x=0.5,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0.5,0))
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0))
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0))
  }
})
title(xlab = "Recall", mgp = c(1.8,0.7,0), xpd = NA, outer = T)
title(ylab = "Precision", mgp = c(2,0.7,0), xpd = NA, outer = T)
legend("bottomright", lty = 1, lwd = c(1,2.5), bty = "n",
       col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"))
dev.off()

# 7. plot of coefficients
png(paste0("fig/modelEvaluation/GiniOverChrsAllTissues.png"),
    width=1200, height=900, pointsize=20)
par(mar = c(3,0,1.5,0.1), oma = c(0,11,0.1,0.1),mfrow = c(1,10))
plotDump = sapply(tissues, function(tissue){
  coeffs = rf_gini[[tissue]]
  bxplt=boxplot(t(coeffs), las = 1, drop =T,yaxt = "n", xaxt = "n",
                horizontal = T, cex.axis = 0.75, main = t2T[tissue], 
                xlim = c(3,nrow(coeffs)-2), names = NA, xaxs = "i")
  axis(1,labels = F)
  axis(1, mgp = c(2,0.5,0), cex.axis = 0.75, 
       at = par("xaxp")[1], tick=F)
  axis(1, mgp = c(2,0.5,0), cex.axis = 0.75, 
       at = par("xaxp")[2], hadj = 0.8,  tick=F)
  if(tissue == tissues[1]){
    axis(2,at = 1:nrow(coeffs), labels=rownames(coeffs), las = 1,
         cex.axis = 0.75)
  }
  segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
           x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
})
title(xlab = "Corrected gini importance", 
      xpd = NA, outer = T, line = -1.5)
title(ylab = "Predictors", outer = T, line = 10)
dev.off()

#png(paste0("fig/modelEvaluation/ImpOverChrsAllTissues.png"),
#    width=1200, height=900, pointsize=20)
#par(mar = c(3,0,1.5,0.1), oma = c(0,11,0.1,0.1),mfrow = c(1,10))
#plotDump = sapply(tissues, function(tissue){
#   coeffs = rf_imps[[tissue]]
#   bxplt=boxplot(t(coeffs), las = 1, drop =T,yaxt = "n", xaxt = "n",
#                 horizontal = T, cex.axis = 0.75, 
#                 xlim = c(3,nrow(coeffs)-2), names = NA, xaxs = "i")
#   title(main = t2T[tissue], xpd = NA)
#   axis(1,labels = F)
#   axis(1, mgp = c(2,0.5,0), cex.axis = 0.75, 
#        at = par("xaxp")[1], tick=F)
#   axis(1, mgp = c(2,0.5,0), cex.axis = 0.75, 
#        at = par("xaxp")[2], hadj = 0.8,  tick=F)
#   if(tissue == tissues[1]){
#      axis(2,at = 1:nrow(coeffs), labels=rownames(coeffs), las = 1,
#           cex.axis = 0.75)
#   }
#   segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
#            x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
#})
#title(xlab = "Permutation importance", 
#      xpd = NA, outer = T, line = -1.5)
#title(ylab = "Predictors", outer = T, line = 10)
#dev.off()
#####


# collection of plots comparison with glm #####
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
# coefficients against coefficients
#rf_imp = rf_gini[[tissue]]
#rf_imp = rf_imp[apply(rf_imp,1,function(x)!all(is.na(x))),]
#glm_imp = -log10(glm_pvals[[tissue]])
## glm_imp = abs(glm_imps[[tissue]])
#glm_imp = glm_imp[apply(glm_imp,1,function(x)!all(is.na(x))),]
#
#commonPreds = intersect(rownames(rf_imp), rownames(glm_imp))
#rf_imp = rf_imp[commonPreds,]
#glm_imp = glm_imp[commonPreds,]
#plot(rf_imp, glm_imp, las = 1, pch = 1, cex = 0.5,xaxt = "n",yaxt = "n",
#     col =  NA, xlab = "RF gini importance", ylab = "",
#     mgp = c(2,0.7,0), tcl=-0.3) #, xaxt = "n"
#axis(2, mgp = c(2,0.7,0), las = 1,
#     at = par("yaxp")[1], lwd = 0, lwd.ticks=1)
#axis(2, mgp = c(2,0.7,0), las = 1,padj=0.7,
#     at = par("yaxp")[2], lwd = 0, lwd.ticks=1)
#rfMeans = rowMeans(rf_imp, na.rm = T)
#glmMeans = rowMeans(glm_imp, na.rm = T)
#rfSD = apply(rf_imp,1,sd, na.rm = T)
#glmSD = apply(glm_imp,1,sd, na.rm = T)
#arrows(x0 = rfMeans-rfSD, x1 = rfMeans+rfSD, 
#       y0 = glmMeans, y1 = glmMeans, code=0)
#arrows(x0 = rfMeans, x1 = rfMeans, 
#       y0 = glmMeans-glmSD, y1 = glmMeans+glmSD, code=0)
#points(rfMeans, glmMeans,
#       bg = predictorCols[rownames(rf_imp)],
#       pch = 21, col = "black")
#abline(lm(glmMeans~rfMeans, na.action="na.exclude"))
#if(tissue == tissues[length(tissues)]){
#   title(xlab = "RF gini importance", xpd = NA, 
#         mgp = c(2,0.7,0), cex.lab = 1.2)
#}
#if(tissue == tissues[ceiling(length(tissues)/2)]){
#   title(ylab = "-log10(glm p-value)", xpd = NA, 
#         mgp = c(2,1,0), cex.lab = 1.2)
#}
#outliers = which(rfMeans > par("xaxp")[2]*0.66)
#if(length(outliers)>1){
#   text(p2P[names(rfMeans[outliers])], pos=2, cex = 0.7,
#        x=rfMeans[outliers], y = glmMeans[outliers]+0.005)
#}
#outliers = which(glmMeans > par("yaxp")[2]*0.66)
#if(length(outliers)>1){
#   text(p2P[names(rfMeans[outliers])], pos=4, cex = 0.7,
#        x=rfMeans[outliers], y = glmMeans[outliers]+0.005)
#}

#####

# comparison of tissues performance #####
png("fig/modelEvaluation/compareRFperformanceTissues.png", 
    width=2400, height=800, pointsize=35)
par(mfrow = c(1,3), mar = c(4,4,0.1,0.1), mgp = c(2.5,1,0))
# roc
plot(NA, xlim = c(0,1), ylim = c(0,1), 
     xlab = "False positive rate", 
     ylab = "True positive rate", las = 1)
abline(0,1, col = "dark grey", lty = 2)
plotDump = sapply(tissues, function(tissue){
  lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]], 
        ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2)
})
legend("bottomright", col=tissueCols, legend=t2T[names(tissueCols)], lwd = 2)
# pr
plot(NA, xlim = c(0,1), ylim = c(0,1), 
     xlab = "Recall", 
     ylab = "Precision", las = 1)
plotDump = sapply(tissues, function(tissue){
  lines(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]], 
        ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2)
})
# AUCs vs nMuts
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
text(x=nMuts, y = colMeans(AUCsPerChr), labels=t2T[tissues], pos=4)
dev.off()


# what feature correlates most with AUROC?
png("fig/modelEvaluation/AUROCvsNmutsComparison.png",
    width=2400, height=800, pointsize=35)
par(mfrow = c(1,3))
# plot again: AUCs vs nMuts per tissue
plot(nMuts, colMeans(AUCsPerChr), las = 1,
     bg = tissueCols, pch = 21, cex = 1.5, 
     ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
              max(colMeans(AUCsPerChr)+AUCsSD)), 
     xlim = c(min(nMuts), max(nMuts)+65000), 
     xlab = "Number of positions", ylab = "AUROC")
abline(lm(colMeans(AUCsPerChr) ~ nMuts), col = "grey", lty = 2)
p = summary(lm(colMeans(AUCsPerChr) ~ nMuts))$coefficients[2,4]
legend("bottomright", legend=paste0("p-value = ", format(p,digits=3)))
arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nMuts,
       y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nMuts, code = 3,
       angle = 90, length = 0.1)
points(nMuts, colMeans(AUCsPerChr),
       bg = tissueCols, pch = 21, cex = 2)
text(x=nMuts, y = colMeans(AUCsPerChr), labels=t2T[tissues], pos=4, offset=1)
# AUC vs nMuts per sample
medianNmuts =  sapply(nPerTissue, median)
quantNmuts = sapply(nPerTissue, quantile,probs=c(0.25, 0.75))
plot(medianNmuts, colMeans(AUCsPerChr), las = 1,
     bg = tissueCols, pch = 21, cex = 1.5, 
     xlim = c(min(quantNmuts), max(quantNmuts)),
     ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
              max(colMeans(AUCsPerChr)+AUCsSD)), 
     xlab = "Number of  mutations per sample", ylab = "AUROC")
abline(lm(colMeans(AUCsPerChr) ~ medianNmuts), col = "grey", lty = 2)
p = summary(lm(colMeans(AUCsPerChr) ~ medianNmuts))$coefficients[2,4]
legend("bottomright", legend=paste0("p-value = ", format(p,digits=3)))
arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = medianNmuts,
       y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = medianNmuts, code = 3,
       angle = 90, length = 0.1)
arrows(y0 = colMeans(AUCsPerChr), x0 = quantNmuts[2,],
       y1 = colMeans(AUCsPerChr), x1 = quantNmuts[1,],
       code = 3, angle = 90, length = 0.1)
points(medianNmuts, colMeans(AUCsPerChr),
       bg = tissueCols, pch = 21, cex = 2)
text(x=medianNmuts, y = colMeans(AUCsPerChr), labels=t2T[tissues], pos=4, offset=1)
# AUC vs nSamples per tissue
nSamples = sapply(nPerTissue, length)
plot(nSamples, colMeans(AUCsPerChr), las = 1,
     bg = tissueCols, pch = 21, cex = 1.5, 
     xlim = c(min(nSamples), max(nSamples)+150),
     ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
              max(colMeans(AUCsPerChr)+AUCsSD)), 
     xlab = "Number of samples per tissue", ylab = "AUROC")
abline(lm(colMeans(AUCsPerChr) ~ nSamples), col = "grey", lty = 2)
p = summary(lm(colMeans(AUCsPerChr) ~ nSamples))$coefficients[2,4]
legend("bottomright", legend=paste0("p-value = ", format(p,digits=3)))
arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nSamples,
       y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nSamples, code = 3,
       angle = 90, length = 0.1)
points(nSamples, colMeans(AUCsPerChr),
       bg = tissueCols, pch = 21, cex = 2)
text(x=nSamples, y = colMeans(AUCsPerChr), labels=t2T[tissues], pos=4, offset=1)
dev.off()
#####


# permutation importance heatmaps #####
# raw permutation importance values
#meanImps = do.call(rbind,lapply(tissues, function(tissue){
#   imp = rowMeans(rf_imps[[tissue]])
#   d = data.frame(importance = imp,tissue)
#   d$predictor = rownames(d)
#   return(d)
#}))
#meanImps$predictor = factor(meanImps$predictor,levels=predictorOrder)
#ggplot(meanImps, aes(x = tissue, y = predictor)) + 
#   geom_raster(aes(fill=importance)) +
#   scale_fill_gradient(low="grey90", high="red",na.value="grey") #+
#ggsave("fig/RFimportances_heatmap.png", height=8, width=7)
# removing effect and TFs only present in liver
#subMeanImps = meanImps[!meanImps$predictor %in% c(rarePredictors, "effect"),]
#ggplot(subMeanImps, aes(x = tissue, y = predictor)) + 
#   geom_raster(aes(fill=importance)) +
#   scale_fill_gradient(low="grey90", high="red",na.value="grey") #+
#ggsave("fig/RFimportances_heatmap_noEffect_noTFs.png", height=8, width=7)
# scale across tissues
#meanImpsScaled = do.call(rbind,lapply(tissues, function(tissue){
#   imp = rowMeans(rf_imps[[tissue]])
#   imp = (imp-min(imp, na.rm = T))/(max(imp, na.rm = T)-min(imp, na.rm = T))
#   d = data.frame(importance = imp, tissue)
#   d$predictor = rownames(d)
#   return(d)
#}))
#meanImpsScaled$predictor = factor(meanImpsScaled$predictor,
#                                  levels=predictorOrder)
#ggplot(meanImpsScaled, aes(x = tissue, y = predictor)) + 
#   geom_raster(aes(fill=importance)) +
#   scale_fill_gradient(low="grey90", high="red",na.value="grey") 
#ggsave("fig/RFimportances_heatmapScaled.png", height=8, width=7)
# scale + remove effect+TFs
#subMeanImpsScaled = meanImpsScaled[!meanImpsScaled$predictor %in% c(rarePredictors, "effect"),]
#ggplot(subMeanImpsScaled, aes(x = tissue, y = predictor)) + 
#   geom_raster(aes(fill=importance)) +
#   scale_fill_gradient(low="grey90", high="red",na.value="grey") #+
#ggsave("fig/RFimportances_heatmapScaled_noEffect_noTFs.png", height=8, width=7)

# coeff barplot
#impsMelted = do.call(rbind,lapply(names(rf_imps), function(tissue){
#   x = as.data.frame(rf_imps[[tissue]])
#   long <- reshape(x, direction = "long",
#                   ids = row.names(x), idvar = "predictors", 
#                   times = names(x), timevar = "chr",
#                   varying = names(x), v.names = "importance")
#   long$tissue = tissue
#   return(long)
#}))
#impsMelted$predictor = factor(impsMelted$predictor,
#                                  levels=predictorOrder)
#ggplot(impsMelted, aes(x = predictors, y = importance, col = tissue)) +
#   geom_boxplot(show.legend=F)  + 
#   coord_flip() +
#   facet_wrap(~ tissue, ncol=10,scales="free_x") +
#   theme_minimal() #+
#ggsave("fig/RFimportances_barplot.png", height=8, width=13)


#####

# gini importance heatmaps #####
meanImps = do.call(rbind,lapply(tissues, function(tissue){
  imp = rowMeans(rf_gini[[tissue]], na.rm = T)
  d = data.frame(gini_corrected = imp, tissue)
  d$predictor = rownames(d)
  return(d)
}))
meanImps$predictor = factor(meanImps$predictor,
                            levels=predictorOrder)
meanImps$gini_corrected = log(meanImps$gini_corrected+1)
ggplot(meanImps, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=gini_corrected)) +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Gini\nimportance")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(meanImps$predictor)])+
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))
ggsave("fig/RFgini_heatmap.png", height=8, width=7)
# remove effect + TFs
#meanImpsSub = meanImps[!meanImps$predictor %in% c("effect", rarePredictors),]
meanImpsSub = meanImps[!meanImps$predictor %in% rarePredictors,]
ggplot(meanImpsSub, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=gini_corrected)) +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Gini\nimportance")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(meanImpsSub$predictor)])+
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))
ggsave("fig/RFgini_heatmap_noTFs.png", height=8, width=7)
# scale 
meanImpsScaled = do.call(rbind,lapply(tissues, function(tissue){
  imp = rowMeans(rf_gini[[tissue]])
  imp = (imp-min(imp, na.rm = T))/(max(imp, na.rm = T)-min(imp, na.rm = T))
  d = data.frame(gini_corrected = imp, tissue)
  d$predictor = rownames(d)
  return(d)
}))
meanImpsScaled$predictor = factor(meanImpsScaled$predictor,
                                  levels=predictorOrder)
meanImpsScaled$gini_corrected = log(meanImpsScaled$gini_corrected+1)
ggplot(meanImpsScaled, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=gini_corrected)) +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Gini\nimportance")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(meanImpsScaled$predictor)])+
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))
ggsave("fig/RFgini_heatmapScaled.png", height=8, width=7)
# scale after removing effect+TFs
subMeanImpsScaled = do.call(rbind,lapply(tissues, function(tissue){
  imp = rowMeans(rf_gini[[tissue]], na.rm = T)
  #imp = imp[!names(imp) %in% c(rarePredictors, "effect")]
  imp = imp[!names(imp) %in% rarePredictors]
  imp = (imp-min(imp, na.rm = T))/(max(imp, na.rm = T)-min(imp, na.rm = T))
  d = data.frame(gini_corrected = imp, tissue)
  d$predictor = rownames(d)
  return(d)
}))
subMeanImpsScaled$predictor = factor(subMeanImpsScaled$predictor,
                                     levels=predictorOrder)
subMeanImpsScaled$gini_corrected = log(subMeanImpsScaled$gini_corrected+1)

ggplot(subMeanImpsScaled, aes(x = tissue, y = predictor)) + 
  geom_raster(aes(fill=gini_corrected)) +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Tissue") + labs(fill = "Gini\nimportance")+ 
  scale_x_discrete(labels=t2T, guide = guide_axis(n.dodge = 2)) + 
  scale_y_discrete(labels = p2P[as.character(subMeanImpsScaled$predictor)])+
  theme(plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0))
ggsave("fig/RFgini_heatmapScaled_noTFs.png", height=8, width=7)
#####



# predictor versus mutation rate - effect direction #####
preds = predictorOrder[!predictorOrder %in% rarePredictors]
for(tissue in tissues){
  print(tissue)
  load(paste0("data/procData/traindata/traindata_processed_",
              tissue, ".RData"))
  png(paste0("fig/predictor_effectDirection_", tissue, ".png"), 
      width=2150, height=1350, pointsize=20)
  par(mfrow = c(6,10), mar = c(2.5,3,1.5,0.5))
  sapply(preds, function(x){
    if(! x %in% colnames(dat)){
      plot.new()
    }else{
      nLevels = length(unique(dat[,x]))
      if(nLevels == 2){
        pd = cbind(table(dat[,x][dat$mutated == 0])/sum(dat$mutated == 0), NA, NA)
        barplot(pd, names.arg=c("TN", "TP", ""),
                las = 1, main = p2P[x], col = c("firebrick 2", "dark red"))
        pd = cbind(NA, table(dat[,x][dat$mutated == 1])/sum(dat$mutated == 1), NA)
        barplot(pd, col = c( "steelblue3", "dark blue"), add= T, yaxt = "n", width=c(1,1,0.2))
        legend("topright", fill = c("grey90", "grey50"), legend=c(0,1), 
               bty = "n")
      } else if (nLevels <5){
        barplot(cbind(table(dat[,x][dat$mutated == 0]),
                      table(dat[,x][dat$mutated == 1]), NA), 
                names.arg=c("TN", "TP", ""),
                las = 1, main = p2P[x], col = rainbow(nLevels), legend.text=T,
                args.legend = list(bty = "n"))
      }
      else   if(is.numeric(dat[,x])){
        temp = density(log(dat[,x][dat$mutated == 1]))
        temp2 = density(log(dat[,x][dat$mutated == 0]))
        plot(temp, col = "firebrick 2", lwd = 2, main = p2P[x], xlab = "", las = 1,
             ylim = c(0,max(c(temp$y, temp2$y))))
        lines(temp2, col = "dark blue", lwd = 2, lty = 2)
      }  else
        plot(1:10, main = p2P[x]) 
    }
  })
  
  
  dev.off()
}
hist(dat$effect[dat$mutated == 1], ylim = c(-15000, 15000), col = "dark red",
     main = "", xlab = "effect", yaxt = "n")
temp = hist(dat$effect[dat$mutated == 0], plot=F)
temp$counts = -temp$counts
plot(temp, add = T, col = "dark blue")
test = axTicks(2)
axis(2, at=test, labels=abs(test), las = 2)
legend("topleft", legend=c("TP", "TN"), fill = c("dark red", "dark blue"), cex = 1.4, bty="n")
######



# Compare new, data w/o UTR with the old one #####
png(paste0("fig/modelEvaluation/compare_OldVSnoUTR_", tissue, ".png"),
    width=1200, height=1000, pointsize = 25)
m = rbind(c(1,1,1,2,2,2,3,3,3), 
          c(1,1,1,2,2,2,3,3,3), 
          c(1,1,1,2,2,2,3,3,3), 
          c(4,4,4,2,2,2,3,3,3), 
          c(4,4,4,2,2,2,3,3,3),
          c(4,4,4,2,2,2,3,3,3))
par(oma = c(0.1,0.1,1.3,0.5), mar = c(4,4,0.5,1))
layout(m)

# ROC
plot(ROC_PR_RF_concat[[tissue]]$roc, 
     col = "red",lwd = 2)
plot(ROC_PR_RF_concat_old[[tissue]]$roc, 
     col = "black", add = T, lwd = 2)
abline(coef = c(0,1), lty = 2)
legend("bottomright", lty = 1, bty = "n",
       col = c("red", "black"), legend=c("AUC_noUTR", "AUC_old"))

# No UTR gini imp
par(mar = c(4,9.5,0,0))
coeffs = rf_gini[[tissue]]
coeffs = coeffs[apply(coeffs,1,function(x){!all(is.na(x))}),]
bxplt=boxplot(t(coeffs), las = 1, drop =T,
              horizontal = T, cex.axis = 0.75, xaxt = "n",
              xlim = c(2,nrow(coeffs)-1))
axis(1)
mtext("gini importance", side = 1, cex= 0.8, line = 2.5)
segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
         x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
text(x = 1400, y = 45, "No UTR")

# Old RF gini imp
par(mar = c(4,9.5,0,0))
coeffs_old = rf_gini_old[[tissue]]
coeffs_old = coeffs_old[rownames(coeffs_old) %in% rownames(coeffs),]
bxplt=boxplot(t(coeffs_old), las = 1, drop =T,
              horizontal = T, cex.axis = 0.75, xaxt = "n",
              xlim = c(2,nrow(coeffs)-1))
axis(1)
mtext("gini importance", side = 1, cex= 0.8, line = 2.5)
segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
         x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
text(x = 10000, y = 45, "Old data")

# Data sizes
par(mar = c(1,1,1,1))
plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
text(x = 0.34, y = 0.9, paste0("nMuts_noUTR: ", dataInfos[[tissue]]$nMuts, "\n nMuts_old: ", dataInfos_old[[tissue]]$nMuts), 
     cex = 1.1, col = "black", family="serif", font=1, adj=0.5)
dev.off()
#####

# collection of plots for every tissue old vs UTR#####
# ROC
png(paste0("fig/modelEvaluation/compare_OldVSnoUTR_ROC.png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.7,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  abline(0,1, col = "grey", lty = 2)
  plot(ROC_PR_RF_concat[[tissue]]$roc, 
       col = "red", add = T)
  plot(ROC_PR_RF_concat_old[[tissue]]$roc, 
       col = "black", add = T)
  text(x=0.05,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0,0))
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0))
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0))
  }
})
title(ylab = "TPR", mgp = c(2,0.7,0), xpd = NA, outer=T)
title(xlab = "FPR", mgp = c(1.8,0.7,0), xpd = NA, outer=T)
legend("bottomright", lty = 1, bty = "n",
       col = c("red", "black"), legend=c("No UTR", "Old data"))
dev.off()
# PR
png(paste0("fig/modelEvaluation/compare_OldVSnoUTR_PR.png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.7,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  plot(ROC_PR_RF_concat[[tissue]]$pr, 
       col = "red", add = T)
  plot(ROC_PR_RF_concat_old[[tissue]]$pr, 
       col = "black", add = T)
  text(x=0.5,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0.5,0))
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0))
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0))
  }
})
title(xlab = "Recall", mgp = c(1.8,0.7,0), xpd = NA, outer = T)
title(ylab = "Precision", mgp = c(2,0.7,0), xpd = NA, outer = T)
legend("bottomright", lty = 1, bty = "n",
       col = c("red", "black"), legend=c("No UTR", "Old data"))
dev.off()
# Coefs
png(paste0("fig/modelEvaluation/compare_OldVSnoUTR_coeffs_", tissue, ".png"),
    width=1200, height=1000, pointsize = 25)
par(mfrow = c(1,2), par(mar = c(3,9.5,1,1)))
coeffs = rf_gini[[tissue]]
coeffs = coeffs[apply(coeffs,1,function(x){!all(is.na(x))}),]
bxplt=boxplot(t(coeffs), las = 1, drop =T,
              horizontal = T, cex.axis = 0.6, xaxt = "n",
              xlim = c(2,nrow(coeffs)-1), main = "no UTR")
axis(1)
segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
         x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")

coeffs_old = rf_gini_old[[tissue]]
coeffs_old = coeffs_old[rownames(coeffs_old) %in% rownames(coeffs),]
bxplt=boxplot(t(coeffs_old), las = 1, drop =T,
              horizontal = T, cex.axis = 0.6, xaxt = "n",
              xlim = c(2,nrow(coeffs)-1), main = "old data")
axis(1)
segments(y0= 1:nrow(coeffs), y1 = 1:nrow(coeffs),
         x0=-5,x1 = bxplt$stats[1,], lty = 2,col = "grey")
title(xlab = "Gini importance", xpd = NA, outer = T, line = -1)
dev.off()
#####