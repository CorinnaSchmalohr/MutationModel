source("scripts/05_analysis/00_NamesAndColors.R")
library(ranger)
library(ROCR)
library(ggplot2)
library(RColorBrewer)
library(Biostrings)
chroms = paste0("chr", 1:22)

# compare number of reference bases between tissues ####
nRefs = sapply(tissues, function(tissue){
   cat(tissue, ' ')
   load(paste0("data/procData/traindata/traindata_processed_",
              tissue, ".RData"))
   table(dat$context_ref)
})
colnames(nRefs) = t2T[colnames(nRefs)]
baseCols = brewer.pal(4,"Dark2")
png("fig/nBasesPerTissue.png")
par(mar = c(6,5,2,1))
barplot(nRefs, col = baseCols, las = 2, ylab = "n bases",legend.text=T,
        args.legend=list(x="topleft", bty = "n"), mgp = c(4,0.8,0))
dev.off()

png("fig/percBasesPerTissue.png")
percRefs = apply(nRefs,2,function(x){x/sum(x)})
barplot(percRefs, col = baseCols, las = 2, ylab = "porportion of bases", mgp = c(4,0.8,0))
dev.off()
######


# only test data is specific base  #####
print("ref base for test data")
testRef = sapply(c("A", "C", "G", "T"), function(base){
   print(paste0("reference base", base))
   crossAUCs = sapply(tissues, function(tissue2predict){
      cat("data: ", tissue2predict, "\n")
      load(paste0("data/procData/traindata/traindata_processed_",
                  tissue2predict, ".RData"))
      tempChroms = datchroms[dat$context_ref == base]
      tempData = dat[dat$context_ref == base,]
      colnames(tempData)[colnames(tempData) == "aPhased_repeats"] = "aPhased_repeates"
      rm(dat, datchroms)
      cat("model: ")
      AUCs = sapply(tissues, function(predictingTissue){
         cat(predictingTissue, " ")
         load(paste0("data/procData/traindata/traindata_processed_",
                     predictingTissue, ".RData"))
         # because we have so little data, we have to use all chromosomes
         preds = sapply(chroms, function(cr){
            cat(cr, ' ')
            testData = tempData[tempChroms == cr,]
            truePreds = testData$mutated
            load(paste0("data/rdata/RFmodel/", predictingTissue,
                        "_", cr, ".RData"))
            # make sure we have same set of predictors
            predictors = names(rf$variable.importance)
            testData = testData[,colnames(testData) %in% predictors]
            toAdd = predictors[!predictors %in% colnames(testData)]
            for(x in toAdd){
               testData[x] = rep(mean(dat[,x]), nrow(testData))
            }
            # predict
            yHat = predict(rf,data=testData,type="response", num.threads=10)
            return(cbind(yHat$predictions[,2], truePreds))
         }, simplify = F)
         preds = do.call(rbind,preds)
         cat('\n')

         temp = prediction(preds[,1], preds[,2])
         auc = performance(temp,"auc")
         return(auc@y.values[[1]])
      }, simplify=F)
      cat("\n")
      return(AUCs)
   }, simplify = F)
   save(crossAUCs, file = paste0("data/rdata/CrossTissue_testRef", base, ".RData"))
   return(crossAUCs)
})
#####

# only training data is specific base  #####
print("ref base for training data")

trainRef = sapply(c("A", "C", "G", "T"), function(base){
   print(paste0("reference base", base))
   crossAUCs = sapply(tissues, function(tissue2predict){
      cat("data: ", tissue2predict, "\n")
      load(paste0("data/procData/traindata/traindata_processed_",
                  tissue2predict, ".RData"))
      # tempChroms = datchroms[dat$context_ref == base]
      # tempData = dat[dat$context_ref == base,]
      tempChroms = datchroms
      tempData = dat
      colnames(tempData)[colnames(tempData) == "aPhased_repeats"] = "aPhased_repeates"
      rm(dat, datchroms)
      cat("model: ")
      AUCs = sapply(tissues, function(predictingTissue){
         cat(predictingTissue, " ")
         load(paste0("data/procData/traindata/traindata_processed_",
                     predictingTissue, ".RData"))
         # because we have so little data, we have to use all chromosomes
         preds = sapply(chroms, function(cr){
            cat(cr, ' ')
            testData = tempData[tempChroms == cr,]
            truePreds = testData$mutated
            # train model only on this reference base
            # load(paste0("data/rdata/RFmodel/", predictingTissue,
            #             "_", cr, ".RData"))
            trainData = dat[dat$context_ref == base & datchroms != cr,]
            rf = ranger(mutated ~ ., data = trainData, 
                        write.forest = T, seed = 1234, num.threads =  10,
                        respect.unordered.factors = 'partition',
                        probability = T, verbose=F, importance = "impurity_corrected")
            save(rf, file = paste0("data/rdata/RFmodel/", predictingTissue,
                                   "_", cr, "_ref", base,".RData"))
            # make sure we have same set of predictors
            predictors = rf$forest$independent.variable.names
            testData = testData[,colnames(testData) %in% predictors]
            toAdd = predictors[!predictors %in% colnames(testData)]
            for(x in toAdd){
               testData[x] = rep(mean(dat[,x]), nrow(testData))
            }
            # predict
            yHat = predict(rf,data=testData,type="response", num.threads=10)
            return(cbind(yHat$predictions[,2], truePreds))
         }, simplify = F)
         preds = do.call(rbind,preds)
         cat('\n')
         
         temp = prediction(preds[,1], preds[,2])
         auc = performance(temp,"auc")
         return(auc@y.values[[1]])
      }, simplify=F)
      cat("\n")
      return(AUCs)
   }, simplify = F)
   save(crossAUCs, file = paste0("data/rdata/CrossTissue_trainRef", base, ".RData"))
   return(crossAUCs)
})
#####

#  training AND testing data is specific base  #####
print("ref base for both")
trainBoth = sapply(c("A", "C", "G", "T"), function(base){
   print(paste0("reference base", base))
   crossAUCs = sapply(tissues, function(tissue2predict){
      cat("data: ", tissue2predict, "\n")
      load(paste0("data/procData/traindata/traindata_processed_",
                  tissue2predict, ".RData"))
      tempChroms = datchroms[dat$context_ref == base]
      tempData = dat[dat$context_ref == base,]
      colnames(tempData)[colnames(tempData) == "aPhased_repeats"] = "aPhased_repeates"
      rm(dat, datchroms)
      cat("model: ")
      AUCs = sapply(tissues, function(predictingTissue){
         cat(predictingTissue, " ")
         load(paste0("data/procData/traindata/traindata_processed_",
                     predictingTissue, ".RData"))
         # because we have so little data, we have to use all chromosomes
         preds = sapply(chroms, function(cr){
            cat(cr, ' ')
            testData = tempData[tempChroms == cr,]
            truePreds = testData$mutated
            # train model only on this reference base
            load(paste0("data/rdata/RFmodel/", predictingTissue,
                        "_", cr, "_ref", base,".RData"))
            # make sure we have same set of predictors
            predictors = rf$forest$independent.variable.names
            testData = testData[,colnames(testData) %in% predictors]
            toAdd = predictors[!predictors %in% colnames(testData)]
            for(x in toAdd){
               testData[x] = rep(mean(dat[,x]), nrow(testData))
            }
            # predict
            yHat = predict(rf,data=testData,type="response", num.threads=10)
            return(cbind(yHat$predictions[,2], truePreds))
         }, simplify = F)
         preds = do.call(rbind,preds)
         cat('\n')
         
         temp = prediction(preds[,1], preds[,2])
         auc = performance(temp,"auc")
         return(auc@y.values[[1]])
      }, simplify=F)
      cat("\n")
      return(AUCs)
   }, simplify = F)
   save(crossAUCs, file = paste0("data/rdata/CrossTissue_bothRef", base, ".RData"))
   return(crossAUCs)
})
#####


# load and visualize ######
testRefAUCs = do.call(rbind,lapply(c("A", "C", "G", "T"), function(base){
   load(paste0("data/rdata/CrossTissue_testRef", base, ".RData"))
   temp = do.call(rbind,lapply(names(crossAUCs), function(tissue2predict){
      x = unlist(crossAUCs[[tissue2predict]])
      data.frame(testData = tissue2predict, trainData = names(x), AUC = unname(x), base = base, rank = NA)
   }))
   for(tissue in tissues){
     temp[temp$testData == tissue,]$rank[order(-temp[temp$testData == tissue,]$AUC)] = 1:length(tissues)
   }
   return(temp)
}))
ggplot(testRefAUCs, aes(trainData, testData, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1)) +
   facet_wrap(~ base)+
   labs(y = "Tested on ", x = "Trained on") +  
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)#+
   #geom_text(aes(label = rank), color = "white", size = 2.5)
ggsave(filename="fig/crossTissue_testRefAUCs.png", 
       width = 15, height = 15, units = "cm")

trainRefAUCs = do.call(rbind,lapply(c("A", "C", "G", "T"), function(base){
   load(paste0("data/rdata/CrossTissue_trainRef", base, ".RData"))
   temp = do.call(rbind,lapply(names(crossAUCs), function(tissue2predict){
      x = unlist(crossAUCs[[tissue2predict]])
      data.frame(testData = tissue2predict, trainData = names(x), AUC = unname(x), base = base, rank = NA)
   }))
   for(tissue in tissues){
     temp[temp$testData == tissue,]$rank[order(-temp[temp$testData == tissue,]$AUC)] = 1:length(tissues)
   }
   return(temp)
}))
ggplot(trainRefAUCs, aes(trainData, testData, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1)) +
   facet_wrap(~ base)+
   labs(y = "Tested on ", x = "Trained on") +  
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)
ggsave(filename="fig/crossTissue_trainRefAUCs.png", 
       width = 15, height = 15, units = "cm")

bothRefAUCs = do.call(rbind,lapply(c("A", "C", "G", "T"), function(base){
   load(paste0("data/rdata/CrossTissue_bothRef", base, ".RData"))
   temp = do.call(rbind,lapply(names(crossAUCs), function(tissue2predict){
      x = unlist(crossAUCs[[tissue2predict]])
      data.frame(testData = tissue2predict, trainData = names(x), AUC = unname(x), base = base, rank = NA)
   }))
   for(tissue in tissues){
     temp[temp$testData == tissue,]$rank[order(-temp[temp$testData == tissue,]$AUC)] = 1:length(tissues)
   }
   return(temp)
}))
ggplot(bothRefAUCs, aes(trainData, testData, fill= AUC)) +
   geom_tile() +
   scale_fill_viridis_c(limits=c(0.47,0.67)) +
   theme_minimal() +
   theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1)) +
   facet_wrap(~ base)+
   labs(y = "Tested on ", x = "Trained on") +  
   scale_x_discrete(labels=t2T) + 
   scale_y_discrete(labels = t2T)
ggsave(filename="fig/crossTissue_bothRefAUCs.png", 
       width = 15, height = 15, units = "cm")
#####

# compare heatmap rows directly #####
png("fig/crossTissue_compareRefAUCs.png", width=900, height=1400, pointsize=30)
m = cbind(c(1,2,2,2,2),
          c(3,4,4,4,4),
          c(5,6,6,6,6))
par(oma = c(4,7,1.3,0.7), mar = c(0.1,0.1,0.5,0))
layout(m) 
boxplot(AUC ~ testData, data = testRefAUCs, las = 1, ylim= c(0.46, 0.645),
        horizontal=T, cex.axis = 0.9,  ylab = "", xaxt = "n", yaxt = "n")
axis(2, at = 1:10, labels=t2T[tissues], las = 1)
mtext("Test on reference",side = 3, cex = 0.7)

mtext("Tested on tissue data",side = 2, line = 5.5)
boxplot(AUC ~ testData:base, data = testRefAUCs, las = 1, ylim= c(0.46, 0.645),
        horizontal=T, cex.axis = 0.9, col =rep(baseCols, each=10), ylab = "", yaxt = "n") 
axis(2, at=1:40, labels=rep(t2T[tissues],4), las = 1)
abline(h=10*0:4+0.5, lty = 2)
#legend( "bottomright", legend=rev(unique(testRefAUCs$base)),
#        fill=baseCols,cex=0.8)

boxplot(AUC ~ testData, data = trainRefAUCs, las = 1, ylim= c(0.46, 0.645),
        horizontal=T,  ylab = "", xaxt = "n", yaxt ="n")
mtext("Training on reference",side = 3, cex = 0.7)

boxplot(AUC ~ testData:base, data = trainRefAUCs, las = 1, ylim= c(0.46, 0.645),
        horizontal=T,  col =rep(baseCols, each=10), ylab = "", yaxt ="n")
abline(h=10*0:4+0.5, lty = 2)

boxplot(AUC ~ testData, data = bothRefAUCs, las = 1, ylim= c(0.46, 0.645),
        horizontal=T,   ylab = "", xaxt = "n", yaxt ="n")
mtext("Training & test on reference",side = 3, cex = 0.7)

boxplot(AUC ~ testData:base, data = bothRefAUCs, las = 1, ylim= c(0.46, 0.645),
        horizontal=T,  col =rep(baseCols, each=10), ylab = "", yaxt ="n")
abline(h=10*0:4+0.5, lty = 2)
mtext("Tested on base-specific tissue data",side = 2, line = 5.5,outer=T)
mtext("AUC",side = 1, line = 2.3,outer=T)
legend("bottomright", legend=c("all", "T", "G", "C", "A"), 
       fill = c("grey", baseCols[c(4:1)]), cex=1.1)
dev.off()
#####

