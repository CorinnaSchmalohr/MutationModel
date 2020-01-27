# random Forest #####
load("data/rdata/luad/completeData.RData")
data$precedingDimer = NULL
data$followingDimer = NULL
data$trimer = NULL
data$pentamer = NULL
data$septamer = NULL
removed = removed[data$mutated < 2 & 
                     data$repeatMasker == 1,]
data = data[data$mutated < 2 & 
                  data$repeatMasker == 1,]
chroms = unique(removed$chr)
dir.create("fig/luad/RF_variable_importance_reducedData/")
library(ranger)
perf = lapply(chroms, function(cr){
   print(cr)
   rf = ranger(mutated ~ ., data = data, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  10,
               holdout = T, case.weights = as.integer(removed$chr != cr),
               respect.unordered.factors = 'partition')
   
   png(paste0("fig/luad/RF_variable_importance_reducedData/",cr,".png"),
       height = 1500, width = 400)
   par(mar = c(5,10,1,1))
   barplot((rf$variable.importance*1000), las = 1, cex.names = 0.7,
           xlab = "variable importance (*1000)", horiz = T, main = cr)
   dev.off()
   trimsplt = sapply(1:500, function(i){
      used = which(rf$forest$split.varIDs[[i]] == 4)
      split = rf$forest$split.values[[i]][used]
   })
   
   save(rf, file = paste0("data/rdata/luad/RFmodel_class", cr, ".RData"))
   return(list(rf$prediction.error, rf$r.squared, rf$variable.importance))
})
save(perf, file = "data/rdata/luad_RFmodel_class.RData")
#####


load("data/rdata/luad_RFmodel_class.RData")


rss = sapply(perf, function(x){x[[1]]})
mean(rss)

rsqu = sapply(perf, function(x){x[[2]]})
mean(rsqu)

png("fig/luad/RF_variable_importance_only0or1mut.png", height = 800, width = 500)
imp = sapply(perf, function(x){x[[3]]})
par(mar = c(5,11,1,1))
boxplot(t(imp), horizontal = T, las = 1, cex.names = 0.8,
        xlab = "variable importance", names = "", border = "white")
abline(h = 1:49, lty = 2, col = "grey")
abline(v = 0)
boxplot(t(imp), horizontal = T, las = 1, add = T, cex.names = 0.8,
        xlab = "variable importance")
dev.off()

