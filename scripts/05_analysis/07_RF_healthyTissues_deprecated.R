# preparation #####
library(ROCR)
library(ranger) 
nThreads = 6
source("./scripts/05_analysis/00_NamesAndColors.R")
tissues = c("brain","breast", "colon","esophagus", 
            "kidney", "liver", "lung","ovary", 
            "prostate", "skin") 
######


# various tissues #####
print("various Tissues")
tempTissues = c("brain", "colon","esophagus",
                "liver","liver_Blokzijl", "lung",
                "prostate", "skin")
perfsVarious = sapply(tempTissues, function(tissue){
   print(tissue)
   # load healthy tissue trainingData
   load(paste0("data/procData/validationData/validationData_processed_",
               tissue, ".RData")) #dat,datchroms
   #colnames(dat)[colnames(dat) == "aPhased_repeats"] = "aPhased_repeates"
   # iterate through chromosomes
   predPerChr = sapply(unique(datchroms), function(cr){
      cat(cr, ' ')
      # load RF trained on cancer data
      if(tissue == "liver_Blokzijl"){
         load(paste0("data/rdata/RFmodel/liver_", cr, ".RData"))
      } else{
         load(paste0("data/rdata/RFmodel/", tissue, "_", cr, ".RData"))
      }
      # subset test data to chromosome
      testData = dat[datchroms == cr,]
      # get predictions
      p = predict(rf, data = testData,verbose=F)
      temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
      return(temp)
   }, simplify = F); cat('\n')
   # compute ROC per chromosome
   perfPerChr = sapply(predPerChr, function(x){
      if(nrow(x)<3){return(NULL)}
      temp = prediction(x$pred,x$label)
      roc = performance(temp,  "tpr", "fpr")
      auc = performance(temp, "auc")
      pr = performance(temp, "prec", "rec")
      return(list(roc = roc, pr = pr, auc = auc))
   }, simplify = F)
   # compute concatenated ROC
   predConcat = do.call(rbind,predPerChr)  
   temp = prediction(predConcat$pred, predConcat$label)
   roc = performance(temp,  "tpr", "fpr")
   auc = performance(temp, "auc")
   pr = performance(temp, "prec", "rec")
   perfConcat = list(roc = roc, pr = pr, auc = auc)

   return(list(perfPerChr = perfPerChr, perfConcat = perfConcat))
}, simplify = F)
save(perfsVarious, file="data/rdata/perf_variousTissues.RData")
# #####


# GTEx #####
print("GTEx")
tempTissues = c("brain", "breast", "colon", "esophagus", "kidney", "liver",
                "lung", "ovary", "prostate", "skin")
perfsGTEx = sapply(tempTissues,  function(tissue){
   print(tissue)
   # load healthy tissue trainingData
   load(paste0("data/procData/validationData/validationData_GTEx_processed_",
               tissue, ".RData")) #dat,datchroms
   #colnames(dat)[colnames(dat) == "aPhased_repeats"] = "aPhased_repeates"
   # iterate through chromosomes
   predPerChr = sapply(unique(datchroms), function(cr){
      cat(cr, ' ')
      # load RF trained on cancer data
      load(paste0("data/rdata/RFmodel/", tissue, "_", cr, ".RData"))
      # subset test data to chromosome
      testData = dat[datchroms == cr,]
      # get predictions
      p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
      temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
      return(temp)
   }, simplify = F); cat('\n')
   # compute ROC per chromosome
   perfPerChr = sapply(predPerChr, function(x){
      if(nrow(x)<3){return(NULL)}
      temp = prediction(x$pred,x$label)
      roc = performance(temp,  "tpr", "fpr")
      auc = performance(temp, "auc")
      pr = performance(temp, "prec", "rec")
      return(list(roc = roc, pr = pr, auc = auc))
   }, simplify = F)
   # compute concatenated ROC
   predConcat = do.call(rbind,predPerChr)
   temp = prediction(predConcat$pred, predConcat$label)
   roc = performance(temp,  "tpr", "fpr")
   auc = performance(temp, "auc")
   pr = performance(temp, "prec", "rec")
   perfConcat = list(roc = roc, pr = pr, auc = auc)

   return(list(perfPerChr = perfPerChr, perfConcat = perfConcat))
}, simplify = F)
save(perfsGTEx, file="data/rdata/perf_GTEx.RData")
#####


# Moore #####
print("Moore")
tempTissues = c("colon","esophagus", "kidney","liver", "lung", 
                "prostate", "skin")
perfsMoore = sapply(tempTissues, function(tissue){
   print(tissue)
   # load healthy tissue trainingData
   load(paste0("data/procData/validationData/validationData_Moore_processed_",
               tissue, ".RData")) #dat,datchroms
   #colnames(dat)[colnames(dat) == "aPhased_repeats"] = "aPhased_repeates"
   # iterate through chromosomes
   predPerChr = sapply(unique(datchroms), function(cr){
      cat(cr, ' ')
      # load RF trained on cancer data
      load(paste0("data/rdata/RFmodel/", tissue, "_", cr, ".RData")) 
      # subset test data to chromosome
      testData = dat[datchroms == cr,]
      # get predictions
      p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
      temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
      return(temp)
   }, simplify = F); cat('\n')
   # compute ROC per chromosome
   perfPerChr = sapply(predPerChr, function(x){
      if(nrow(x)<3){return(NULL)}
      temp = prediction(x$pred,x$label)
      roc = performance(temp,  "tpr", "fpr")
      auc = performance(temp, "auc")
      pr = performance(temp, "prec", "rec")
      return(list(roc = roc, pr = pr, auc = auc))
   }, simplify = F)
   # compute concatenated ROC
   predConcat = do.call(rbind,predPerChr)
   temp = prediction(predConcat$pred, predConcat$label)
   roc = performance(temp,  "tpr", "fpr")
   auc = performance(temp, "auc")
   pr = performance(temp, "prec", "rec")
   perfConcat = list(roc = roc, pr = pr, auc = auc)
   
   return(list(perfPerChr = perfPerChr, perfConcat = perfConcat))
}, simplify = F)
save(perfsMoore, file="data/rdata/perf_Moore.RData")
######



# Li #####
print("Li")
tempTissues = c("colon","esophagus", "liver", "lung")
perfsLi = sapply(tempTissues, function(tissue){
   print(tissue)
   # load healthy tissue trainingData
   load(paste0("data/procData/validationData/validationData_Li_processed_",
               tissue, ".RData")) #dat,datchroms
   #colnames(dat)[colnames(dat) == "aPhased_repeats"] = "aPhased_repeates"
   # iterate through chromosomes
   predPerChr = sapply(unique(datchroms), function(cr){
      cat(cr, ' ')
      # load RF trained on cancer data
      load(paste0("data/rdata/RFmodel/", tissue, "_", cr, ".RData")) 
      # subset test data to chromosome
      testData = dat[datchroms == cr,]
      # get predictions
      p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
      temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
      return(temp)
   }, simplify = F); cat('\n')
   # compute ROC per chromosome
   perfPerChr = sapply(predPerChr, function(x){
      if(nrow(x)<3){return(NULL)}
      temp = prediction(x$pred,x$label)
      roc = performance(temp,  "tpr", "fpr")
      auc = performance(temp, "auc")
      pr = performance(temp, "prec", "rec")
      return(list(roc = roc, pr = pr, auc = auc))
   }, simplify = F)
   # compute concatenated ROC
   predConcat = do.call(rbind,predPerChr)
   temp = prediction(predConcat$pred, predConcat$label)
   roc = performance(temp,  "tpr", "fpr")
   auc = performance(temp, "auc")
   pr = performance(temp, "prec", "rec")
   perfConcat = list(roc = roc, pr = pr, auc = auc)
   
   return(list(perfPerChr = perfPerChr, perfConcat = perfConcat))
}, simplify = F)
save(perfsLi, file="data/rdata/perf_Li.RData")
#####



# plotting #####
load("data/rdata/perf_variousTissues.RData") # perfsVarious
load("data/rdata/perf_GTEx.RData") #perfsGTEx
load("data/rdata/perf_Moore.RData") #perfsMoore
load("data/rdata/perf_Li.RData") #perfsLi
load("data/rdata/RFmodel/ROC_PR_RF_concat.RData")
load("data/rdata/RFmodel/ROC_PR_RF_perChr.RData")
colGradients = sapply(tissueCols, function(x){
   cl = colorRampPalette(c(x, "black"))(4)[1:3]
})


# plot 1: one panel per tissue, ROC over chrs, thick line concat, 
# dotted line training data performance
dumpVar = sapply(c("perfsVarious", "perfsGTEx", "perfsMoore", "perfsLi"), function(ds){
   print(ds)
   temp = get(ds)
   t = names(temp)
   names(t) = t
   t[t == "liver_Blokzijl"] = "liver"
   
   png(paste0("fig/healthyTissues_",ds, "_perf.png"),
       width=2000, height=1800, res = 200)
   par(mfrow = c(ceiling(length(temp)/2),4), mar = c(3,3,1.5,0.5), oma = c(0.1,0.1,0.1,0.1))
   dumpVar2 = sapply(names(temp), function(tissue){
      cat(tissue, ' ')
      # ROC
      x = temp[[tissue]]$perfPerChr
      plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "FPR",
           ylab = "TPR", mgp = c(2.2,1,0), las = 1)
      abline(0,1, col = "darkgrey", lwd = 3)
      dumpVar3 = sapply(x, function(y){
         if(is.null(y)){return(NA)}
         plot(y$roc, add = T, col = colGradients[1,t[tissue]], lwd = 1)
      })
      x = temp[[tissue]]$perfConcat
      plot(x$roc, add = T, col = colGradients[2,t[tissue]], lwd = 3)
      plot(ROC_PR_RF_concat[[t[tissue]]]$roc, col = colGradients[3,t[tissue]],
           add = T, lwd = 3, lty = 2)

      # PR
      x = temp[[tissue]]$perfPerChr
      plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Recall",
           ylab = "Precision", mgp = c(2.2,1,0), las = 1)
      dumpVar3 = sapply(x, function(y){
         if(is.null(y)){return(NA)}
         plot(y$pr, add = T, col = colGradients[1,t[tissue]], lwd = 1)
      })
      x = temp[[tissue]]$perfConcat
      plot(x$pr, add = T, col = colGradients[2,t[tissue]], lwd = 3)
      plot(ROC_PR_RF_concat[[t[tissue]]]$pr, col =  colGradients[3,t[tissue]],
           add = T, lwd = 3, lty = 2)
      
      mtext(side=3,text= if(tissue %in% names(t2T)){t2T[tissue]} else{"Liver Blokzijl"}, at=-0.2, adj=0.5, line=0.2)
   }); cat('\n')
   dev.off()
})


# plot2:  3 panels.
dumpVar = sapply(c("perfsVarious", "perfsGTEx", "perfsMoore", "perfsLi"), function(ds){
   print(ds)
   temp = get(ds)
   t = names(temp)
   names(t) = t
   t[t == "liver_Blokzijl"] = "liver"
   
   png(paste0("fig/healthyTissues_",ds, "_overview.png"),
       width=2000, height=2000, res = 200)
   par(mfrow = c(1,3), mar = c(3,3,0.5,0.5), oma = c(0.1,0.1,0.1,0.1))
   layout(matrix(c(1,2,3,3), byrow = T, ncol = 2))
   # first panel: ROC concat, colors=tissues, thick line healthy, dotted line training data
   plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "FPR",
        ylab = "TPR", mgp = c(2.2,1,0), las = 1)
   abline(0,1, col = "darkgrey")
   dumpVar2 = sapply(names(t), function(tissue){
      x = temp[[tissue]]$perfConcat
      plot(x$roc, add = T, col = colGradients[1,t[tissue]], lwd = 2)
      plot(ROC_PR_RF_concat[[t[tissue]]]$roc, col = colGradients[2,t[tissue]],
           add = T, lwd = 2, lty = 2)
   })
   legend("bottomright", col = tissueCols[t], legend = t2T[names(t)], lwd = 2)
   # second panel: same as first panel, only PR
   plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Recall",
        ylab = "Precision", mgp = c(2.2,1,0), las = 1)
   dumpVar2 = sapply(names(t), function(tissue){
      x = temp[[tissue]]$perfConcat
      plot(x$pr, add = T, col = colGradients[1,t[tissue]], lwd = 3)
      plot(ROC_PR_RF_concat[[t[tissue]]]$pr, col =  colGradients[2,t[tissue]],
           add = T, lwd = 3, lty = 2)
   })
   legend("bottomright", col = c("grey70", "grey35"), lty = 1:2,
          legend = c("test", "training"), lwd = 2)
   # third panel: mean CV AUCs as barplot. tissues next to each other in different colors.
   plotDat2 = do.call(rbind,sapply(names(t), function(tissue){
      print(tissue)
      hRes = sapply(temp[[tissue]]$perfPerChr,function(y){
         if(!is.null(y)){y$auc@y.values[[1]]}else{NA}
      })
      cRes = sapply(ROC_PR_RF_perChr[[t[tissue]]], function(x){x$auc})
      rbind(data.frame(tissue, source = "healthy", auc = hRes),
            data.frame(tissue, source = "cancer", auc = cRes))
   }, simplify = F))
   #par(mar = c(4,3,0.5,8))
   bp = boxplot(auc ~ source:tissue,plotDat2, plot = F)
   bxp(bp, at = rep(1:length(t)*2, each = 2)+c(-0.3,0.3), las = 2,
       boxfill = colGradients[c(3,1),t],  ylim = c(0,1),
       outline = F, xlim = c(1.5,length(t)*2), show.names=F) #
   title(ylab = "AUC", line = 2.2)
   text(y=-0.08, x = 1:length(t)*2,labels = t2T[names(t)], xpd = NA,  srt = 35, adj = 1)
   abline(v = 0:length(t)*2 +1, lty = 2, col = "grey")
   concatAUCs = sapply(temp, function(x){
      x$perfConcat$auc@y.values[[1]]
   })
   concatAUCs2 = sapply(t, function(tissue){
      ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
   })
   points(x = 1:length(t)*2+0.3, y = concatAUCs, bg = colGradients[1,t], pch = 21)
   points(x = 1:length(t)*2-0.3, y = concatAUCs2, bg = colGradients[3,t], pch = 21)
   legend("bottomright", legend=c( "training", "test", "concat"), 
          fill = c("grey35", "grey70", NA),border=c(1,1,NA), xpd = NA,
          pch = c(NA, NA, 21), pt.bg=c(NA, NA, "white"))
   dev.off()
})


# plot 4: perf vs. data size
load("data/rdata/GTEx_nPerSample.RData")
load("data/rdata/Moore_nPerSample.RData")
load("data/rdata/Li_nPerSample.RData")
load("data/rdata/healthyTissues_nPerSample.RData")
#####

