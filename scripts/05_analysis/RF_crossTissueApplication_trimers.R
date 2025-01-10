
# preparation ##### 
source("scripts/05_analysis/00_NamesAndColors.R")
library(ranger)
library(ROCR)
library(Biostrings)
library(ggplot2)
chroms = paste0("chr", 1:22)

# prepare trimers and count per tissue
bases = c("A", "C", "G", "T")
trimers = sort(apply(expand.grid(bases, bases, bases),1,paste0, collapse = ""))
compl = cbind(trimers,
              as.character(reverseComplement(DNAStringSet(trimers))))
compl = t(apply(compl,1,sort))
compl = unique(compl)
rownames(compl) = paste(compl[,1], compl[,2], sep = "/")
#####


# only test data for specific trimer  #####
print("subset testData to trimer")
testTrimer = sapply(tissues, function(tissue2predict){
   # load test Data
   cat("testData: ", tissue2predict, "\n")
   load(paste0("data/procData/traindata/traindata_processed_",
               tissue2predict, ".RData"))
   testDat = dat
   testChroms = datchroms
   testTrimers = paste0(testDat$context_precedingBase, 
                        testDat$context_ref,
                        testDat$context_followingBase)
   rm(dat, datchroms)
   cat("trainData: ")
   res = sapply(tissues, function(predictingTissue){
      cat(predictingTissue,' ')
      
      # load train Data
      load(paste0("data/procData/traindata/traindata_processed_",
                  predictingTissue, ".RData"))
      trainDat = dat
      trainChroms = datchroms
      rm(dat, datchroms)
      
      # make sure we have same set of predictors
      predictors = colnames(trainDat)
      testDat = testDat[,colnames(testDat) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testDat)]
      for(x in toAdd){
         testDat[x] = rep(mean(trainDat[,x]), nrow(testDat))
      }
      
      # iterate through chromosomes and get predictions
      predsPerChr = sapply(chroms, function(cr){
         # load  RF
         load(paste0("data/rdata/RFmodel/", predictingTissue,
                     "_", cr, ".RData"))
         # subset test Data to chr
         subDat = testDat[testChroms == cr,]
         subTrimers = testTrimers[testChroms == cr]
         # get predictions on chromosome
         yHat = predict(rf,data=subDat,type="response", num.threads=10)
         # return 3 columns: prediction, true value, trimer
         return(data.frame(prediction = yHat$predictions[,2], 
                           trueVal = subDat$mutated, 
                           trimer = subTrimers))
      }, simplify=F)
      predsPerChr = do.call(rbind,predsPerChr)
      
      # get performance per trimer
      ROCperTrimer = apply(compl,1,function(trimer){
         subPreds = predsPerChr[predsPerChr$trimer %in% trimer,]
         temp = prediction(predictions=subPreds$prediction, 
                           labels=subPreds$trueVal)
         auc = performance(temp,"auc")
         return(auc@y.values[[1]])
      })
      return(ROCperTrimer)
   })
   cat('\n')
   return(res)
}, simplify = F)
save(testTrimer, file = "data/rdata/crossTissue_testTrimer.RData")
#####


# only train data for specific trimer  #####
print("subset trainData to trimer")
trainTrimer = sapply(tissues, function(predictingTissue){
   cat("trainData: ", predictingTissue, "\n")
   load(paste0("data/procData/traindata/traindata_processed_",
               predictingTissue, ".RData"))
   trainDat = dat
   trainChroms = datchroms
   rm(dat, datchroms)
   datTrimer = paste0(trainDat$context_precedingBase, 
                      trainDat$context_ref,
                      trainDat$context_followingBase)
   # train a model for only one trimer, excluding chr1 and chr2
   dumpTrimers = apply(compl,1,function(trimer){
      # train a model only on one trimer
      subtrainDat = trainDat[datTrimer %in% trimer & 
                                !trainChroms %in% c("chr1", "chr2"),]
      rf = ranger(mutated ~ ., data = subtrainDat, 
                  write.forest = T, seed = 1234, num.threads =  10,
                  respect.unordered.factors = 'partition',
                  probability = T, verbose=F, importance = "impurity_corrected")
      save(rf, file = paste0("data/rdata/RFmodel/", predictingTissue,
                             "_", trimer[1], "_",trimer[2],
                             "_chr1and2.RData"))
      return(NA)
   })
   
   # iterate through test tissues and test trimer-specific models on whole data
   cat("testData: ")
   ROCs = sapply(tissues, function(tissue2predict){
      # load test Data
      cat(tissue2predict,' ')
      load(paste0("data/procData/traindata/traindata_processed_",
                  tissue2predict, ".RData"))
      testDat = dat[datchroms %in% c("chr1", "chr2"),]
      rm(dat, datchroms)
      # make sure we have same set of predictors
      predictors = colnames(trainDat)
      testDat = testDat[,colnames(testDat) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testDat)]
      for(x in toAdd){
         testDat[x] = rep(mean(trainDat[,x]), nrow(testDat))
      }
      # iterate through trimers and get predictions
      ROCsperTrimer = apply(compl,1,function(trimer){
         # train a model only on one trimer
         load(paste0("data/rdata/RFmodel/", predictingTissue,
                     "_", trimer[1], "_",trimer[2],
                     "_chr1and2.RData"))
         yHat = predict(rf,data=testDat,type="response", num.threads=10)
         temp = prediction(predictions=yHat$predictions[,2], 
                           labels=testDat$mutated)
         auc = performance(temp,"auc")
         return(auc@y.values[[1]])
      })
      return(ROCsperTrimer)
   })
   cat('\n')
   return(ROCs)
}, simplify = F)
save(trainTrimer, file = "data/rdata/crossTissue_trainTrimer.RData")
######


# both train and test data for specific trimer  #####
bothTrimer = sapply(tissues, function(predictingTissue){
   cat("trainData: ", predictingTissue, "\n")
   load(paste0("data/procData/traindata/traindata_processed_",
               predictingTissue, ".RData"))
   trainDat = dat
   trainChroms = datchroms
   rm(dat, datchroms)
   datTrimer = paste0(trainDat$context_precedingBase, 
                      trainDat$context_ref,
                      trainDat$context_followingBase)
   
   # iterate through test tissues and test trimer-specific models on whole data
   cat("testData: ")
   ROCs = sapply(tissues, function(tissue2predict){
      # load test Data
      cat(tissue2predict, " ")
      load(paste0("data/procData/traindata/traindata_processed_",
                  tissue2predict, ".RData"))
      testDat = dat[datchroms %in% c("chr1", "chr2"),]
      testTrimers = paste0(testDat$context_precedingBase,
                           testDat$context_ref,
                           testDat$context_followingBase)
      rm(dat, datchroms)
      # make sure we have same set of predictors
      predictors = colnames(trainDat)
      testDat = testDat[,colnames(testDat) %in% predictors]
      toAdd = predictors[!predictors %in% colnames(testDat)]
      for(x in toAdd){
         testDat[x] = rep(mean(trainDat[,x]), nrow(testDat))
      }
      # iterate through trimers and get predictions
      ROCsperTrimer = apply(compl,1,function(trimer){
         subtestDat = testDat[testTrimers %in% trimer,]
         # train a model only on one trimer
         load(paste0("data/rdata/RFmodel/", predictingTissue,
                     "_", trimer[1], "_",trimer[2],
                     "_chr1and2.RData"))
         yHat = predict(rf,data=subtestDat,type="response", num.threads=10)
         temp = prediction(predictions=yHat$predictions[,2], 
                           labels=subtestDat$mutated)
         auc = performance(temp,"auc")
         return(auc@y.values[[1]])
      })
      return(ROCsperTrimer)
   })
   cat('\n')
   return(ROCs)
}, simplify = F)
save(bothTrimer, file = "data/rdata/crossTissue_bothTrimer.RData")

######


# collect ROCs ######
# order trimers by # 1. reference base (C or T), 
# 2. following base, 3. preceding base
head(compl)
trims = t(apply(compl, 1,function(x){
   temp = which(substr(x,2,2) %in% c("C", "T"))
   refTrim = x[temp]
   otherTrim = x[-temp]
   toDisplay = paste(x[temp], x[-temp], sep= "/")
   prec = substr(refTrim,1,1) 
   ref = substr(refTrim,2,2)
   foll = substr(refTrim,3,3)
   return(c(toDisplay, refTrim, otherTrim, prec, ref, foll))
}))
colnames(trims) = c("toDisplay", "trimer1", "trimer2","prec", "ref", "foll")
trimOrder = trims[order(trims[,"ref"], trims[,"prec"], trims[,"foll"]),]
# test data specific trimer
load("data/rdata/crossTissue_testTrimer.RData")
testROCs = do.call(rbind,sapply(tissues, function(tissue2predict){
   temp = testTrimer[[tissue2predict]]
   temp2 = do.call(rbind,sapply(colnames(temp), function(tis) {
      data.frame(testData = tissue2predict, trainData = tis, 
                 trimer = rownames(temp), AUC = temp[,tis])
   }, simplify=F))
   return(temp2)
}, simplify = F))
testPlot = cbind(testROCs, trimer2 = trimOrder[testROCs$trimer,"toDisplay"],
                 preceding = trimOrder[(testROCs$trimer),"prec"],
                 reference = trimOrder[testROCs$trimer,"ref"],
                 following = trimOrder[testROCs$trimer,"foll"])
# train data specific trimer
load("data/rdata/crossTissue_trainTrimer.RData")
trainROCs = do.call(rbind,sapply(tissues, function(predictingTissue){
   temp = trainTrimer[[predictingTissue]]
   temp2 = do.call(rbind,sapply(colnames(temp), function(tis) {
      data.frame(testData = tis, trainData = predictingTissue,
                 trimer = rownames(temp), AUC = temp[,tis])
   }, simplify=F))
   return(temp2)
}, simplify = F))
trainPlot = cbind(trainROCs, trimer2 = trimOrder[trainROCs$trimer,"toDisplay"],
                 preceding = trimOrder[(trainROCs$trimer),"prec"],
                 reference = trimOrder[trainROCs$trimer,"ref"],
                 following = trimOrder[trainROCs$trimer,"foll"])
# both data specific trimer
load("data/rdata/crossTissue_bothTrimer.RData")
bothROCs = do.call(rbind,sapply(tissues, function(predictingTissue){
   temp = bothTrimer[[predictingTissue]]
   temp2 = do.call(rbind,sapply(colnames(temp), function(tis) {
      data.frame(testData = tis, trainData = predictingTissue,
                 trimer = rownames(temp), AUC = temp[,tis])
   }, simplify=F))
   return(temp2)
}, simplify = F))
bothPlot = cbind(bothROCs, trimer2 = trimOrder[bothROCs$trimer,"toDisplay"],
                 preceding = trimOrder[(bothROCs$trimer),"prec"],
                 reference = trimOrder[bothROCs$trimer,"ref"],
                 following = trimOrder[bothROCs$trimer,"foll"])
######


# visualize ROCs as heatmaps ######
ggplot(testPlot, aes(trainData, testData, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c() +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1)) +
   facet_wrap(~ trimer2, ncol = 8)
ggsave(filename="fig/crossTissue_testTrimerAUCs.png", 
       width = 27, height = 15, units = "cm")
ggplot(trainPlot, aes(trainData, testData, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c() +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1)) +
   facet_wrap(~ trimer2, ncol = 8)
ggsave(filename="fig/crossTissue_trainTrimerAUCs.png", 
       width = 27, height = 15, units = "cm")
ggplot(bothPlot, aes(trainData, testData, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c() +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1)) +
   facet_wrap(~ trimer2, ncol = 8)
ggsave(filename="fig/crossTissue_bothTrimerAUCs.png", 
       width = 27, height = 15, units = "cm")
#####


# visualize trimer content per tissue #####
bases = c("A", "C", "G", "T")
trimers = sort(apply(expand.grid(bases, bases, bases),1,paste0, collapse = ""))
nTrims = sapply(tissues, function(tissue){
   cat(tissue, ' ')
   load(paste0("data/procData/traindata/traindata_processed_",
               tissue, ".RData"))
   tri = paste0(dat$context_precedingBase, dat$context_ref,
                dat$context_followingBase)
   table(factor(tri, levels=trimers))
})
png("fig/nTrimerPerTissue.png", width=1200, height=1200, pointsize=30)
par(mar = c(4,5,0.2,2))
sub = t(nTrims[compl[,1],])
sub2= t(nTrims[compl[,2],])
trimNames = paste(colnames(sub), colnames(sub2), sep = "/")
barplot(-sub2, col = rainbow(10), horiz=T, las = 1, ylab = "trimer",
        xlab = "n per tissue", legend.text=T, cex.names=0.6, names.arg=trimNames,
        xlim = c(-max(rowSums(nTrims)), max(rowSums(nTrims))), mgp = c(4,1,0),
        args.legend=list(x = "bottomleft", ncol = 2, cex = 0.8, bty = "n"))
title(xlab = "n per tissue", mgp = c(2.4,1,0))
barplot(sub, col = rainbow(10), horiz=T, add = T, xaxt = "n", yaxt = "n")
abline(v = 0, lwd = 2)
dev.off()

png("fig/percTrimerPerTissue.png", width=1200, height=1200, pointsize=30)
par(mar = c(4,5,0.2,1))
percTrims = apply(nTrims,2,function(x){x/sum(x)})
sub = t(percTrims[compl[,1],])
sub2= t(percTrims[compl[,2],])
barplot(-sub2, col = rainbow(10), horiz=T, las = 1, ylab = "trimer",
        xlab = "n per tissue", legend.text=T, cex.names=0.6, names.arg=trimNames,
        xlim = c(-max(rowSums(percTrims)), max(rowSums(percTrims))), mgp = c(4,1,0),
        args.legend=list(x = "bottomleft", ncol = 2, cex = 0.8, bty = "n"))
title(xlab = "n per tissue", mgp = c(2.4,1,0))
barplot(sub, col = rainbow(10), horiz=T, add = T, xaxt = "n", yaxt = "n")
abline(v = 0, lwd = 2)
dev.off()

png("fig/propTrimerPerTissue.png", width=1200, height=1200, pointsize=30)
par(mar = c(6.5,4,1,0.5))
sub = t(percTrims[compl[,1],])
sub2= t(percTrims[compl[,2],])
res = t(sub+sub2)
rownames(res) = paste(colnames(sub), colnames(sub2), sep = "/")
barplot(cbind(res, NA, NA), col = rainbow(nrow(res)), las = 2,
        ylab = "proportion per tissue", mgp = c(2.5,0.7,0),
        legend.text=T, args.legend=list(x = "topright", inset = 0,cex = 0.8,
                                        bty = "n", inset = 0))
dev.off()

nTrims2 = sapply(tissues, function(tissue){
   cat(tissue, ' ')
   load(paste0("data/procData/traindata/traindata_processed_",
               tissue, ".RData"))
   tri = paste0(dat$context_precedingBase, dat$context_ref, dat$context_followingBase)
   tri2 = sapply(tri, function(x){
      trimOrder[trimOrder[,"trimer1"] == x | trimOrder[,"trimer2"] == x,"toDisplay"]
   })
   table(factor(tri2, levels=trimOrder[,"toDisplay"]))
})
png("fig/nTrimerPerTissue_sorted.png", width=1200, height=1200, pointsize=30)
par(mar = c(4,6.5,1,0.5))
temp = barplot(t(nTrims2[nrow(nTrims2):1,]), col = rainbow(10), horiz=T,
               las = 1, ylab = "trimer",
               xlab = "", legend.text=T, cex.names=0.6, 
               mgp = c(5,1,0), axisnames = F, 
               args.legend=list(x = "bottomright", ncol = 2, bty = "n"))
axis(2,at = temp, labels = rownames(nTrims)[nrow(nTrims2):1], las = 1, 
     cex.axis= 0.8, family = "mono", font = 2)
title(xlab = "n per tissue", line=2.5)
dev.off()
######

# plot ROC vs data size (trimer content) #####
nTrims2 = sapply(tissues, function(tissue){
   cat(tissue, ' ')
   load(paste0("data/procData/traindata/traindata_processed_",
               tissue, ".RData"))
   tri = paste0(dat$context_precedingBase, dat$context_ref, dat$context_followingBase)
   tri2 = sapply(tri, function(x){
      trimOrder[trimOrder[,"trimer1"] == x | trimOrder[,"trimer2"] == x,"toDisplay"]
   })
   table(factor(tri2, levels=trimOrder[,"toDisplay"]))
})
tissueCols = setNames(rainbow(length(tissues)), nm=tissues)
trimerCols = setNames(rainbow(nrow(trimOrder)), trimOrder[,"toDisplay"])
# reformat plotData, to include data sizes
for(dat in c("testPlot", "trainPlot", "bothPlot")){
   temp = get(dat)
   temp$nTrainData = apply(temp, 1, function(x){
      nTrims2[x["trimer2"],x["trainData"]]
   })
   temp$nTestData = apply(temp, 1, function(x){
      nTrims2[x["trimer2"],x["testData"]]
   })
   assign(paste0(dat,"2"),temp)
}

png("fig/crossTissue_trimer_AUCvstrainDatasize.png",
    width = 1200, height=1000, pointsize=20)
par(mfrow = c(3,4), mar = c(4,4,1,0), oma = c(0,2,2,0))
for(colBy in c("testData", "trainData", "trimer2")){
   for(dat in c("testPlot2", "trainPlot2", "bothPlot2")){
      temp = get(dat)
      if(colBy == "trimer2"){colPal = trimerCols}else{colPal = tissueCols}
      colPal = rgb(col2rgb(colPal)[1,], col2rgb(colPal)[2,],col2rgb(colPal)[3,], 
          alpha = 150, maxColorValue=255, names=names(colPal))
      plot(AUC ~ nTrainData, data=temp, las = 1,
           col = colPal[temp[,colBy]], cex = 0.5, log = "x", 
            mgp= c(2.5,1,0))
      if(dat == "testPlot2"){
         nameTranslate = c("testData" = "colored by test data tissue",
                           "trainData" = "colored by training data tissue", 
                           "trimer2" = "colored by trimer")
         mtext(nameTranslate[colBy], side = 2, line = 4)
      }
      if(colBy == "testData"){
         nameTranslate = c("testPlot2" = "test data trimer-specific", 
                           "trainPlot2" = "training data trimer-specific", 
                           "bothPlot2" = "test and training data trimer-specific")
         mtext(nameTranslate[dat], side = 3, line = 1)
      }
   }
   plot.new()
   if(colBy == "trimer2"){ce = 0.75}else{ce = 1.2}
   legend("topleft", legend = names(colPal), pch = 1, col = colPal, ncol = 2,
          inset = 0, cex = ce)
}
dev.off()

png("fig/crossTissue_trimer_AUCvstestDatasize.png",
    width = 1200, height=1000, pointsize=20)
par(mfrow = c(3,4), mar = c(4,4,1,0), oma = c(0,2,2,0))
for(colBy in c("testData", "trainData", "trimer2")){
   for(dat in c("testPlot2", "trainPlot2", "bothPlot2")){
      temp = get(dat)
      if(colBy == "trimer2"){colPal = trimerCols}else{colPal = tissueCols}
      colPal = rgb(col2rgb(colPal)[1,], col2rgb(colPal)[2,],col2rgb(colPal)[3,], 
                   alpha = 150, maxColorValue=255, names=names(colPal))
      plot(AUC ~ nTestData, data=temp, las = 1,
           col = colPal[temp[,colBy]], cex = 0.5, log = "x", 
           mgp= c(2.5,1,0))
      if(dat == "testPlot2"){
         nameTranslate = c("testData" = "colored by test data tissue",
                           "trainData" = "colored by training data tissue", 
                           "trimer2" = "colored by trimer")
         mtext(nameTranslate[colBy], side = 2, line = 4)
      }
      if(colBy == "testData"){
         nameTranslate = c("testPlot2" = "test data trimer-specific", 
                           "trainPlot2" = "training data trimer-specific", 
                           "bothPlot2" = "test and training data trimer-specific")
         mtext(nameTranslate[dat], side = 3, line = 1)
      }
   }
   plot.new()
   if(colBy == "trimer2"){ce = 0.75}else{ce = 1.2}
   legend("topleft", legend = names(colPal), pch = 1, col = colPal, ncol = 2,
          inset = 0, cex = ce)
}
dev.off()
#####
