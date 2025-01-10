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
library(ranger)
cr = "chr1"  

rf = ranger(mutated ~ ., data = data, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads =  10,
            holdout = T, case.weights = as.integer(removed$chr != cr),
            respect.unordered.factors = 'partition')
splitsRef = sapply(1:rf$num.trees, function(i){
   structure = treeInfo(rf, tree = i)
   table(factor(structure$splitval[structure$splitvarName == "ref"],
                levels = c("1", "1,2", "1,2,3", "1,3", "2", "2,3", "3")))
})
rownames(splitsRef) = sapply(rownames(splitsRef), function(x){
   y = as.integer(unlist(strsplit(x,split = ",")))
   paste(levels(data$ref)[y], collapse = ",")
})


data2 = data
old2new = c("A" = "D", "C" = "B", "G" = "C", "T" = "A")
new2old = names(old2new)
names(new2old) = old2new
levels(data2$ref) = old2new[levels(data2$ref)]
rf2 = ranger(mutated ~ ., data = data2, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads =  10,
            holdout = T, case.weights = as.integer(removed$chr != cr),
            respect.unordered.factors = 'partition')
splitsRef2 = sapply(1:rf2$num.trees, function(i){
   structure = treeInfo(rf2, tree = i)
   table(factor(structure$splitval[structure$splitvarName == "ref"],
                levels = c("1", "1,2", "1,2,3", "1,3", "2", "2,3", "3")))
})
rownames(splitsRef2) = sapply(rownames(splitsRef2), function(x){
   y = as.integer(unlist(strsplit(x,split = ",")))
   z = new2old[levels(data2$ref)][y]
   paste(z[order(z)], collapse = ",")
})

par(mfrow = c(2,1))
barplot(rowSums(splitsRef),
        ylab = "number of splits", main = "refence base")
barplot(rowSums(splitsRef2[order(rownames(splitsRef2)),]),
        ylab = "number of splits", main = "refence base")


