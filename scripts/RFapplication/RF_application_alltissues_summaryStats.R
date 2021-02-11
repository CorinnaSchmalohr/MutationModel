library(corrplot)
library(ggplot2)
library(data.table)
tissues = c("luad", "skin", "colon", "ovary",
            "kidney", "prostate", "breast")
date = "2020-01-14"


# compare predictor importances between tissues #####
imps = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/RFmodel/imp.RData"))
   return(rowMeans(imp))
})
preds = unique(unlist(sapply(imps, names)))
imps = sapply(imps, function(x){
   res = sapply(preds, function(i){x[i]})
   names(res) = preds
   return(res)
})
dat = melt(imps)
colnames(dat) = c("predictor", "tissue", "importance")
p = ggplot(dat, aes(x = tissue, y = predictor)) + 
   geom_raster(aes(fill=importance)) + 
   scale_fill_gradient(low="grey90", high="red")
ggsave(paste0("fig/predictorImportance_tissueComparison_",
              date,".png"),
       p, width=8, height=10)
imps = scale(imps)
dat = melt(imps)
colnames(dat) = c("predictor", "tissue", "importance")
p = ggplot(dat, aes(x = tissue, y = predictor)) + 
   geom_raster(aes(fill=importance)) + 
   scale_fill_gradient(low="grey90", high="red")
ggsave(paste0("fig/predictorImportance_tissueComparison_scaled_",
              date,".png"),
       p, width=8, height=10)
#####

# mutation type per tissue #####
png(paste0("fig/mutType_per_tissues_", date,".png"),
    width=2000, height=800, pointsize=25)
layout(matrix(c(1,1:7),1,8))
temp = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   mutType = paste0(data$ref[data$mutated > 0], ">",
                    removed$alt[data$mutated > 0])
   mutType = table(mutType)
   dat = cbind(mutType[1:6], -mutType[12:7])
   if(tissue == "luad"){
      labels = paste(names(mutType)[1:6], 
                     names(mutType)[12:7], sep=" / ")
      par(mar = c(2,8,3,1))
   }else{
      labels = NA
      par(mar = c(2,1.5,3,1.5))
   }
   barplot(dat[,1], horiz=T, col = rainbow(6),  main = tissue, 
           las = 1, xlim = c(-max(dat),max(dat)),names.arg=labels,
           cex.axis=1.3,cex.names=1.3 )
   barplot(dat[,2], add = T, horiz = T, col = rainbow(6),
           xaxt = "n", yaxt = "n")
})
dev.off()
#####

# corrplot, variable importance per tissue, data size,
# and pred. performance #####
summaryStats = lapply(tissues, function(tissue){
   print(tissue)
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   
   # corrplot
   png(paste0("fig/", tissue, "/corrplot_allpredictors_", tissue, "_",date,".png"),
       height = 1200, width = 1200, pointsize = 30)
   corrplot(cor(data[,sapply(data, is.numeric) & colnames(data) != "inexon"]),
            tl.col = "black", tl.cex = 0.6)
   dev.off()
   
   # variable importance per tissue
   load(paste0("data/rdata/", tissue, "/RFmodel/imp.RData"))
   png(paste0("fig/", tissue, "/RF_variable_importance_", tissue, "_",date,".png"),
       height = 800, width = 500)
   par(mar = c(5,11,1,1))
   boxplot(t(imp), horizontal = T, las = 1, cex.names = 0.8,
           xlab = "variable importance (*1000)", 
           border = "white")
   abline(h = 1:49, lty = 2, col = "grey")
   abline(v = 0)
   boxplot(t(imp), horizontal = T, las = 1, 
           add = T, cex.names = 0.8, names = NA,
           xlab = "variable importance (*1000)")
   dev.off()
   
   
   # performance
   load(paste0("data/rdata/", tissue, "/RFmodel/perf.RData"))
   OOBpredError = unlist(perf)

   return(list(dataSize = nrow(data), OOBpredError = OOBpredError))
})
png(paste0("fig/RFmodel_alltissues_summaryStats_",date,".png"), 
    width=800, height=500, pointsize=25)
layout(matrix(c(1,1,1,2,2),1,5))
par(mar = c(5,5,1,1))
datasize = sapply(summaryStats, function(x){x$dataSize})
names(datasize) = tissues
barplot(datasize, las = 1, xlab = "n positions", horiz=T)
par(mar = c(5,2,1,1))
accuracy = sapply(summaryStats, function(x){1-x$OOBpredError})
colnames(accuracy) = tissues
boxplot(accuracy, las = 1,names.arg=NA, horizontal = T, 
        xlab = "Prediction accuracy\n(1 - OOB prediction Error)")
dev.off()
#####


# scaled vs. unscaled data TODO #####
library(reshape2)

for(tissue in tissues){
   load(paste0("data/rdata/", tissue, "/completeData_unscaled.RData"))
   toScale = sapply(data,is.numeric) & !sapply(data,is.integer)
   unscaled = data[,toScale]
   load(paste0("data/rdata/", tissue, "/completeData.RData"))
   scaled = data[,toScale]
   rm(data, removed)
   temp = melt(scaled)
   temp2 = melt(unscaled)
   temp$scale = "scaled"
   temp2$scale = "unscaled"
   dat = rbind(temp, temp2)
   boxplot(value~variable:scale, data = dat, las = 2,
           boxwex = 0.5, col = c("orange", "yellow"),
           sep = ":", lex.order = TRUE,  yaxs = "i")
}
#####


# are predictions better than random? #####
# permutations

# ROC and PR

# TPs higher prediction than TNs?

#####