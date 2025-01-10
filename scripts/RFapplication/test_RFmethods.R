tissue = "breast"

load(paste0("data/rdata/", tissue, "/completeData.RData"))

data$precedingDimer = NULL
data$followingDimer = NULL
data$trimer = NULL
data$pentamer = NULL
data$septamer = NULL
data$context = NULL

library(ranger)
cr = "chr3"
# regression
data1 = data
rf1 = ranger(mutated ~ ., data = data1, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads =  10,
            holdout = T, case.weights = as.integer(removed$chr != cr),
            respect.unordered.factors = 'partition', 
            scale.permutation.importance = T)
# classification
data2 = data
data2$mutated = as.factor(data2$mutated)
rf2 = ranger(mutated ~ ., data = data2, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads =  10,
            holdout = T, case.weights = as.integer(removed$chr != cr),
            respect.unordered.factors = 'partition', 
            scale.permutation.importance = T)
#classification with probability
data3 = data2
rf3 = ranger(mutated ~ ., data = data3, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads =  10,
            holdout = T, case.weights = as.integer(removed$chr != cr),
            respect.unordered.factors = 'partition', 
            scale.permutation.importance = T, probability = T)
