load("data/rdata/luad/completeData.RData")
removed = removed[data$mutated < 2 & data$repeatMasker == 1,]
data = data[data$mutated < 2 & data$repeatMasker == 1,]
data$repeatMasker = NULL
library(ranger)
chroms = paste0("chr", 1:22)
names(chroms) = chroms

# one model per reference base #####
print("bases")
bases = c("A", "C", "T", "G")
names(bases) = bases

perfPerRef = lapply(bases, function(b){
   cat(b,' ')
   subdata = data[data$ref == b,]
   subremoved = removed[data$ref == b,]
   subdata$ref = NULL
   subdata$precedingDimer = NULL
   subdata$followingDimer = NULL
   subdata$trimer = NULL
   subdata$pentamer = NULL
   subdata$septamer = NULL
   perf = lapply(chroms, function(cr){
      rf = ranger(mutated ~ ., data = subdata, importance = 'permutation',
                  write.forest = T, seed = 1234, num.threads =  10,
                  holdout = T, case.weights = as.integer(subremoved$chr != cr),
                  respect.unordered.factors = 'partition')
      return(list(rf$variable.importance, rf$prediction.error, rf$r.squared))
   })
   return(perf)
})
save(perfPerRef, file = "data/rdata/luad/RFmodel_perfPerRef.RData")
#####


# one model per trimer #####
print("trimers")
trimers = names(table(data$trimer))
names(trimers) = trimers
chroms = cbind(paste0("chr", 1:11), paste0("chr", 22:12))
perfPerTrimer = lapply(trimers, function(tri){
   cat(tri,' ')
   subdata = data[data$trimer == tri,]
   subremoved = removed[data$trimer == tri,]
   subdata$trimer = NULL
   subdata$precedingDimer = NULL
   subdata$followingDimer = NULL
   subdata$pentamer = NULL
   subdata$septamer = NULL
   subdata$ref = NULL
   subdata$precedingBase = NULL
   subdata$followingBase = NULL
   perf = apply(chroms,1, function(cr){
      rf = ranger(mutated ~ ., data = subdata, importance = 'permutation',
                  write.forest = T, seed = 1234, num.threads =  10,
                  holdout = T, 
                  case.weights = as.integer(!subremoved$chr %in% cr),
                  respect.unordered.factors = 'partition')
      return(list(rf$variable.importance, rf$prediction.error, rf$r.squared))
   })
   return(perf)
})
save(perfPerTrimer, file = "data/rdata/luad/RFmodel_perfPerTrimer.RData")
#####


# one model per pentamer #####
print("pentamer")
pents = table(data$pentamer)
pents = pents[pents >= 500]
pents = names(pents)
names(pents) = pents
perfPerPentamer = lapply(pents, function(pent){
   cat(pent, ' ')
   subdata = data[data$pentamer == pent,]
   subremoved = removed[data$pentamer == pent,]
   subdata$pentamer = NULL
   subdata$precedingDimer = NULL
   subdata$followingDimer = NULL
   subdata$trimer = NULL
   subdata$septamer = NULL
   subdata$ref = NULL
   subdata$precedingBase = NULL
   subdata$followingBase = NULL
   perf = apply(chroms, 1, function(cr){
      rf = ranger(mutated ~ ., data = subdata, importance = 'permutation',
                  write.forest = T, seed = 1234, num.threads =  10,
                  holdout = T, 
                  case.weights = as.integer(!subremoved$chr %in% cr),
                  respect.unordered.factors = 'partition')
      return(list(rf$variable.importance, rf$prediction.error, rf$r.squared))
   })
   return(perf)
})
save(perfPerPentamer, file = "data/rdata/luad/RFmodel_perfPerPentamer.RData")
#####



# analyse ref base #####
load("data/rdata/luad/RFmodel_perfPerRef.RData")
library(ggplot2)
library(reshape2)

dat = melt(table(data$ref))
p1 = ggplot(data = dat, aes(x = Var1, y = value, fill = Var1)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = "reference base", y = "n positions") +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12, color = "black"),
         legend.position = "none")


imp = sapply(perfPerRef, function(x){
   temp = sapply(x, function(y){y[[1]]})
   rowMeans(temp)
})
imp = melt(imp)
imp$Var2 = as.character(imp$Var2)
colnames(imp) = c("predictor", "refBase", "importance")
p2 = ggplot(imp, aes(x = refBase, y = predictor)) + 
   geom_raster(aes(fill=importance)) + 
   scale_fill_gradient(low="grey90", high="red")

pE = sapply(perfPerRef, function(x){
   sapply(x, function(y){
      y[[2]]
   })
})
dat = melt(pE)
dat$Var2 = as.character(dat$Var2)
p3 = ggplot(data = dat, aes(x = Var2, y = value, fill = Var2)) +
   geom_boxplot() +
   theme_classic() +
   labs(x = "reference base", y = "prediction error") +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12, color = "black"),
         legend.position = "none")



library(gridExtra)
gs = list(p1, p2, p3)
lay <- rbind(c(1,2,2),
             c(3,2,2),
             c(NA,2,2))
g = arrangeGrob(grobs = gs, layout_matrix = lay)
ggsave("fig/luad/RFmodel_perRef.png", plot=g)
#####

# split type on following per ref base #####
bases = c("A", "C", "T", "G")
names(bases) = bases
par(mfrow = c(2,2))
rfPerRef = lapply(bases, function(b){
   cat(b,' ')
   subdata = data[data$ref == b,]
   subremoved = removed[data$ref == b,]
   subdata$ref = NULL
   subdata$precedingDimer = NULL
   subdata$followingDimer = NULL
   subdata$trimer = NULL
   subdata$pentamer = NULL
   subdata$septamer = NULL
   rf = ranger(mutated ~ ., data = subdata, importance='none',
               write.forest = T, seed = 1234, num.threads =  10,
               respect.unordered.factors = 'partition')
   splitsFollowing = sapply(1:rf$num.trees, function(i){
      structure = treeInfo(rf, tree = i)
      table(factor(structure$splitval[structure$splitvarName == "followingBase"],
                   levels = c("1", "1,2", "1,2,3", "1,3", "2", "2,3", "3")))
   })
   barplot(rowSums(splitsFollowing), las = 2,
           names.arg = c("A", "A & C", "A, C & G", "A & G", "C", "C & G", "G"),
           ylab = "number of splits", 
           main = paste0("splits on following base\nfor ref base ",b))
   return(splitsFollowing)
})
#####

# analyse trimer #####
load("data/rdata/luad/RFmodel_perfPerTrimer.RData")
library(ggplot2)
library(reshape2)

dat = melt(table(data$trimer))
p1 = ggplot(data = dat, aes(x = Var1, y = value, fill = Var1)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = "trimer", y = "n positions") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 10, color = "black", 
                                    angle = 90, vjust =0.5, hjust = 1),
         axis.text.y = element_text(size = 7, color = "black"),
         legend.position = "none") +
   coord_flip()


imp = sapply(perfPerTrimer, function(x){
   temp = sapply(x, function(y){y[[1]]})
   rowMeans(temp)
})
imp = melt(imp)
colnames(imp) = c("predictor", "trimer", "importance")
p2 = ggplot(imp, aes(x = trimer, y = predictor)) + 
   geom_raster(aes(fill=importance)) + 
   scale_fill_gradient(low="grey90", high="red") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 6, color = "black",
                                  angle = 90, hjust=1)) 
pE = sapply(perfPerTrimer, function(x){
   sapply(x, function(y){
      y[[2]]
   })
})
dat = melt(pE)
p3 = ggplot(data = dat, aes(x = Var2, y = value, 
                            fill = Var2, color = Var2)) +
   geom_boxplot() +
   theme_classic() +
   labs(x = "trimer", y = "prediction error") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 10, color = "black", 
                                    angle = 90, vjust = 0.5),
         axis.text.y = element_text(size = 7, color = "black"),
         legend.position = "none") +
   coord_flip()



library(gridExtra)
gs = list(p1, p2, p3)
lay <- rbind(c(1,3,2,2,2,2))
g = arrangeGrob(grobs = gs, layout_matrix = lay)
plot(g)
ggsave("fig/luad/RFmodel_perTrimer.png", plot=g, width=12, height=5)
#####

# analyse pentamer #####
load("data/rdata/luad/RFmodel_perfPerPentamer.RData")
library(ggplot2)
library(reshape2)
dat = table(data$pentamer)
dat = dat[dat >= 500]
dat = melt(dat)
p1 = ggplot(data = dat, aes(x = Var1, y = value, fill = Var1)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = "pentamer", y = "n positions") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 10, color = "black", 
                                    angle = 90, vjust = 0.5),
         axis.text.y = element_text(size = 4, color = "black"),
         legend.position = "none") +
   coord_flip()

imp = sapply(perfPerPentamer, function(x){
   temp = sapply(x, function(y){y[[1]]})
   rowMeans(temp)
})
imp = melt(imp)
colnames(imp) = c("predictor", "pentamer", "importance")
p2 = ggplot(imp, aes(x = pentamer, y = predictor)) + 
   geom_raster(aes(fill=importance)) + 
   scale_fill_gradient(low="grey90", high="red") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5))

pE = sapply(perfPerPentamer, function(x){
   sapply(x, function(y){
      y[[2]]
   })
})
dat = melt(pE)
p3 = ggplot(data = dat, aes(x = Var2, y = value, 
                       fill = Var2, color = Var2)) +
   geom_boxplot() +
   theme_classic() +
   labs(x = "trimer", y = "prediction error") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 10, color = "black", 
                                    angle = 90, vjust = 0.5),
         axis.text.y = element_text(size = 4, color = "black"),
         legend.position = "none") +
   coord_flip()

library(gridExtra)
gs = list(p1, p2, p3)
lay <- rbind(c(1,3,2,2,2))
g = arrangeGrob(grobs = gs, layout_matrix = lay)
ggsave("fig/luad/RFmodel_perPentamer.png",
       plot=g, width = 14, height = 6)
#####

