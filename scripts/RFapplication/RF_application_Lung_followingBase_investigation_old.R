setwd("/cellnet/MutationModel/")

# without preceding Base #####
load("data/rdata/luad/completeData.RData")
data$precedingDimer = NULL
data$followingDimer = NULL
data$trimer = NULL
data$pentamer = NULL
data$septamer = NULL

data$precedingBase = NULL
removed = removed[data$mutated < 2 & 
                     data$repeatMasker == 1,]
data = data[data$mutated < 2 & 
               data$repeatMasker == 1,]
chroms = unique(removed$chr)
dir.create("fig/luad/RF_variable_importance_woutPreceding/")
library(ranger)
perf = lapply(chroms[1:4], function(cr){
   print(cr)
   rf = ranger(mutated ~ ., data = data, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  10,
               holdout = T, case.weights = as.integer(removed$chr != cr),
               respect.unordered.factors = 'partition')
   
   save(rf, file = paste0("data/rdata/luad/RFmodel_woutPreceding", cr, ".RData"))
   return(list(rf$prediction.error, rf$r.squared, rf$variable.importance))
})
save(perf, file = "data/rdata/luad_RFmodel_woutPreceding.RData")
#####


# withouth following base #####
load("data/rdata/luad/completeData.RData")
data$precedingDimer = NULL
data$followingDimer = NULL
data$trimer = NULL
data$pentamer = NULL
data$septamer = NULL

data$followingBase = NULL
removed = removed[data$mutated < 2 & 
                     data$repeatMasker == 1,]
data = data[data$mutated < 2 & 
               data$repeatMasker == 1,]
chroms = unique(removed$chr)
dir.create("fig/luad/RF_variable_importance_woutFollowing/")
library(ranger)
perf = lapply(chroms[1:4], function(cr){
   print(cr)
   rf = ranger(mutated ~ ., data = data, importance = 'permutation',
               write.forest = T, seed = 1234, num.threads =  10,
               holdout = T, case.weights = as.integer(removed$chr != cr),
               respect.unordered.factors = 'partition')
   save(rf, file = paste0("data/rdata/luad/RFmodel_woutFollowing", cr, ".RData"))
   return(list(rf$prediction.error, rf$r.squared, rf$variable.importance))
})
save(perf, file = "data/rdata/luad_RFmodel_woutFollowing.RData")
#####

# compare performance of RF model when 
# leaving out preceding and following base #####
load("data/rdata/luad_RFmodel_woutPreceding.RData")
prec = perf
load("data/rdata/luad_RFmodel_woutFollowing.RData")
follow = perf
load("data/rdata/luad_RFmodel_class.RData")
normal = perf

png("fig/luad/perf_wout_precedingAndFollowing.png", 
    width = 850, height = 400, pointsize = 15)
par(mfrow = c(1,2))
precRSS = sapply(prec, function(x){x[[1]]})
followRSS = sapply(follow, function(x){x[[1]]})
normalRSS = sapply(normal, function(x){x[[1]]})
c(mean(precRSS), mean(followRSS), mean(normalRSS))
barplot(c(mean(normalRSS), mean(precRSS), mean(followRSS)), 
        ylab = "mean prediction error", 
        names.arg = c("full model", "preceding\nbase omitted", "following\nbase omitted"))

precRsquared = sapply(prec, function(x){x[[2]]})
followRsquared = sapply(follow, function(x){x[[2]]})
normalRsquared = sapply(normal, function(x){x[[2]]})
barplot(c(mean(normalRsquared), mean(precRsquared), mean(followRsquared)), 
        ylab = "R-squared", 
        names.arg = c("full model", "preceding\nbase omitted", "following\nbase omitted"))
dev.off()
#####



# what are the splits done by the RF on following, preceding, and ref base? #####
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
cr = "chr1"
rf = ranger(mutated ~ ., data = data, importance = 'permutation',
            write.forest = T, seed = 1234, num.threads =  10,
            holdout = T, case.weights = as.integer(removed$chr != cr),
            respect.unordered.factors = 'partition')

splitsPrec = sapply(1:rf$num.trees, function(i){
   structure = treeInfo(rf, tree = i)
   table(factor(structure$splitval[which(structure$splitvarName == "precedingBase")],
                levels = c("1", "1,2", "1,2,3", "1,3", "2", "2,3", "3")))
})
splitsRef = sapply(1:rf$num.trees, function(i){
   structure = treeInfo(rf, tree = i)
   table(factor(structure$splitval[which(structure$splitvarName == "ref")],
                levels = c("1", "1,2", "1,2,3", "1,3", "2", "2,3", "3")))
})
splitsFollowing = sapply(1:rf$num.trees, function(i){
   structure = treeInfo(rf, tree = i)
   table(factor(structure$splitval[which(structure$splitvarName == "followingBase")],
                levels = c("1", "1,2", "1,2,3", "1,3", "2", "2,3", "3")))
})

png("fig/luad/splitType_following_ref_preceding.png", 
    height = 800, width = 500, pointsize = 15)
par(mfrow = c(3,1))
barplot(rowSums(splitsPrec),
        names.arg = c("A", "A & C", "A, C & G", "A & G", "C", "C & G", "G"),
        ylab = "number of splits", main = "preceding base")
barplot(rowSums(splitsRef),
        names.arg = c("A", "A & C", "A, C & G", "A & G", "C", "C & G", "G"),
        ylab = "number of splits", main = "refence base")
barplot(rowSums(splitsFollowing),
        names.arg = c("A", "A & C", "A, C & G", "A & G", "C", "C & G", "G"),
        ylab = "number of splits", main = "following base")
dev.off()




#####



# is the dataset imbalanced in sequence context, i.e. do mutated
# spots have a different sequence context than TNs? #####
load("data/rdata/luad/completeData.RData")
png("fig/luad/Nmutations_by_bases.png")
par(mfrow = c(3,1))
dat = rbind(table(data$ref[data$mutated > 0]),
            table(data$ref[data$mutated == 0]))
barplot(dat, ylab = "n positions", beside = T,
        col = c("red", "blue"), border = F, main = "reference base")
legend("topleft", legend = c("mutated", "TNs"),
       fill = c("red", "blue"), border = F, bty = "n")
dat2 = rbind(table(data$followingBase[data$mutated > 0]),
             table(data$followingBase[data$mutated == 0]))
barplot(dat2, ylab = "n positions", beside = T,
        col = c("red", "blue"), border = F, main = "following Base")
legend("topleft", legend = c("mutated", "TNs"),
       fill = c("red", "blue"), border = F, bty = "n")
dat3 = rbind(table(data$precedingBase[data$mutated > 0]),
             table(data$precedingBase[data$mutated == 0]))
barplot(dat3, ylab = "n positions", beside = T,
        col = c("red", "blue"), border = F, main = "preceding Base")
legend("topleft", legend = c("mutated", "TNs"),
       fill = c("red", "blue"), border = F, bty = "n")
dev.off()

png("fig/luad/Nmutations_by_trimer.png", 
    width = 900, height = 900, pointsize = 18)
par(mfrow = c(3,1))
precedingDimer = substr(data$trimer, start = 1,stop = 2)
dat = rbind(table(precedingDimer[data$mutated > 0]),
             table(precedingDimer[data$mutated == 0]))
barplot(dat, ylab = "n positions", beside = T,
        col = c("red", "blue"), border = F, main = "preceding dimer")
legend("topleft", legend = c("mutated", "TNs"),
       fill = c("red", "blue"), border = F, bty = "n")
followingDimer = substr(data$trimer, start = 2,stop = 3)
dat = rbind(table(followingDimer[data$mutated > 0]),
            table(followingDimer[data$mutated == 0]))
barplot(dat, ylab = "n positions", beside = T,
        col = c("red", "blue"), border = F, main = "following dimer")
legend("topleft", legend = c("mutated", "TNs"),
       fill = c("red", "blue"), border = F, bty = "n")

dat = rbind(table(data$trimer[data$mutated > 0]),
             table(data$trimer[data$mutated == 0]))
barplot(dat,las = 2, beside = T, col = c("red", "blue"),
        border = NA, ylab = "n positions", main = "trimer", ylim = c(0,9000))
containsCC = rep(F, ncol(dat))
containsCC[grep("CC", colnames(dat))] = T
containsCC[grep("GG", colnames(dat))] = T
axis(1, labels = colnames(dat)[containsCC],
     at = ((1:ncol(dat))*3 - 1)[containsCC],
     las = 2, tick = F, col.axis = "green")
lastA = rep(F, ncol(dat))
lastA[substr(colnames(dat), 3,3) == "A"]  = T
axis(1, labels = colnames(dat)[lastA],
     at = ((1:ncol(dat))*3 - 1)[lastA],
     las = 2, tick = F, col.axis = "red")
legend("topleft", legend = c("mutated", "TNs"),
       fill = c("red", "blue"), border = F, bty = "n")
legend(x = 10, y = 9000, legend = c("contains CC/GG", "last pos A"),
       text.col = c("green", "red"), bty = "n")
dev.off()
#####


# pyramid plot #####
# load("data/rdata/luad/completeData.RData")

library(plotrix)
dat = cbind(table(data$trimer[data$mutated == 1]),
            table(data$trimer[data$mutated == 0]))
colnames(dat) = c("mutated", "TNs")

f = dat[1:32,]
r = dat[64:33,]
labels = paste(rownames(f), rownames(r), sep = "/")
rownames(f) = labels
rownames(r) = labels
png("fig/luad/Nmutations_by_trimer_FandR.png", 
    width = 700, height = 700, pointsize = 15)
pyramid.plot(lx = f, rx = r,
             labels = labels,
             top.labels = c("forward", "trimer", "reverse"), gap = 2000, 
             unit = "N positions",lxcol = c("red", "blue"), rxcol = c("red", "blue"), 
             raxlab = seq(0,12500,5000), laxlab =  seq(0,12500,5000))
legend("topright", legend = c("mutated", "TNs"), fill = c("red", "blue"))
dev.off()

png("fig/luad/Nmutations_by_trimer_diffFandR.png", 
    width = 700, height = 700, pointsize = 15)
par(mar = c(4,6,2,1))
dif = f - r
barplot(t(dif),  horiz = T, las = 1, beside = T, col = c("red", "blue"))
abline(h = ((1:32)*3 - 1), lty = 2)
axis(side = 1,at = c(-2000, 2000), 
     labels = c("more in reverse", "more in forward"), line = 1.5, tick = F)
legend("topleft", fill = c("red", "blue"), legend = c("mutated", "TNs"))
dev.off()
#####

# is it feasible to train a model on each pentamer separately?
png("fig/luad/nMut_by_pentamer.png")
par(mfrow = c(2,1))
hist(table(data$pentamer[data$mutated > 0]),
     xlab = "n mutations per pentamer",  main = "pentamer")
hist(table(data$trimer[data$mutated > 0]),
     xlab = "n mutations per trimer",  main = "trimer")
dev.off()