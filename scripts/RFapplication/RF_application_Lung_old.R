load("data/rdata/luad/completeData.RData")


library(corrplot)
png("fig/luad/corrplot_allpredictors.png", height = 1200, width = 1200, pointsize = 30)
corrplot(cor(data[,sapply(data, is.numeric)]),tl.col = "black", tl.cex = 0.6)
dev.off()


mutType = paste0(data$ref[data$mutated > 0], ">",removed$alt[data$mutated > 0])
png("fig/luad/mutationTypes.png", width = 700, height = 400)
barplot(table(mutType), ylab = "n mutations",
        xlab = "mutation type", main = "LUAD mutations",
        col = rep(c("red", "blue", "dark green", "yellow"), each = 3))
dev.off()

data$precedingDimer = NULL
data$followingDimer = NULL
data$trimer = NULL
data$pentamer = NULL
data$septamer = NULL

library(ranger)
# random Forest #####
chroms = unique(removed$chr)
dir.create("fig/luad/RF_variable_importance/")
perf = lapply(chroms, function(cr){
   print(cr)
   rf = ranger(mutated ~ ., data = data, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  10,
               holdout = T, case.weights = as.integer(removed$chr != cr),
               respect.unordered.factors = 'partition')
   save(rf, file = paste0("data/rdata/luad/RFmodel_", cr, ".RData"))
   return(c(rf$prediction.error, rf$r.squared))
})
save(perf, file = "data/rdata/luad_RFmodel.RData")
#####



# load("data/rdata/luad_RFmodel.RData")
rss = do.call(rbind, perf)[,1]
rsqu = do.call(rbind, perf)[,2]
mean(rss)
mean(rsqu)


imp = sapply(chroms, function(cr){
   print(cr)
   load(paste0("data/rdata/luad/RFmodel_", cr, ".RData"))
   return(rf$variable.importance)
})

png("fig/luad/RF_variable_importance_all.png", height = 800, width = 500)
par(mar = c(5,11,1,1))
boxplot(t(imp), horizontal = T, las = 1, cex.names = 0.8,
        xlab = "variable importance (*1000)", names = "", border = "white")
abline(h = 1:49, lty = 2, col = "grey")
abline(v = 0)
boxplot(t(imp), horizontal = T, las = 1, add = T, cex.names = 0.8,
        xlab = "variable importance (*1000)")
dev.off()



png("fig/luad/Nmutations_by_repeatMasker.png")
boxplot(data$mutated ~ !data$repeatMasker,
        xlab = "in repeatMasker region", ylab = "number of mutations")
text(x = 1, y = 10, labels = paste0("n = ",table(!data$repeatMasker)[1]))
text(x = 2, y = 10, labels = paste0("n = ",table(!data$repeatMasker)[2]))
dev.off()
