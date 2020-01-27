# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad",  "skin", "colon", "ovary",
            "kidney", "prostate", "breast")
tissue = tissues[args]
print(tissue)

load(paste0("data/rdata/", tissue, "/completeData.RData"))
ref = data$ref
trimer = data$trimer
pentamer = data$pentamer
data$ref = NULL
data$precedingDimer = NULL
data$followingDimer = NULL
data$trimer = NULL
data$pentamer = NULL
data$septamer = NULL
data$context = NULL
data$mutated = as.factor(data$mutated)
data$inexon = NULL

library(ranger)
chroms = paste0("chr", 1:22)
names(chroms) = chroms

# one model per reference base #####
# print("bases")
# bases = c("A", "C", "T", "G")
# names(bases) = bases
# 
# perfPerRef = lapply(bases, function(b){
#    cat(b,' ')
#    subdata = data[ref == b,]
#    subremoved = removed[ref == b,]
#    perf = lapply(chroms, function(cr){
#       rf = ranger(mutated ~ ., data=subdata, importance='permutation',
#                   write.forest = T, seed = 1234, num.threads =  8,
#                   case.weights = as.integer(subremoved$chr != cr),
#                   holdout = T, respect.unordered.factors='partition', 
#                   scale.permutation.importance = T, probability = T)
#       return(list(rf$variable.importance, rf$prediction.error))
#    })
#    return(perf)
# })
# save(perfPerRef, file = paste0("data/rdata/", tissue, 
#      "/RFmodel_perfPerRef.RData"))
#####


# one model per trimer #####
# print("trimers")
# trimers = names(table(trimer))
# names(trimers) = trimers
# chroms = cbind(paste0("chr", 1:11), paste0("chr", 22:12))
# perfPerTrimer = lapply(trimers, function(tri){
#    cat(tri,' ')
#    subdata = data[trimer == tri,]
#    subremoved = removed[trimer == tri,]
#    perf = apply(chroms,1, function(cr){
#       rf = ranger(mutated ~ ., data=subdata, importance='permutation',
#                   write.forest = T, seed = 1234, num.threads = 8,
#                   holdout = T, respect.unordered.factors='partition', 
#                   case.weights = as.integer(!subremoved$chr %in% cr),
#                   scale.permutation.importance = T, probability = T)
#       return(list(rf$variable.importance, rf$prediction.error))
#    })
#    return(perf)
# })
# save(perfPerTrimer, file = paste0("data/rdata/", tissue,
#                                   "/RFmodel_perfPerTrimer.RData"))
#####


# one model per pentamer #####
# print("pentamer")
# pents = table(pentamer)
# pents = pents[pents >= 500]
# pents = names(pents)
# names(pents) = pents
# perfPerPentamer = lapply(pents, function(pent){
#    cat(pent, ' ')
#    subdata = data[pentamer == pent,]
#    subremoved = removed[pentamer == pent,]
#    perf = apply(chroms, 1, function(cr){
#       rf = ranger(mutated ~ ., data = subdata, importance = 'permutation',
#                   write.forest = T, seed = 1234, num.threads =  10,
#                   holdout = T, 
#                   case.weights = as.integer(!subremoved$chr %in% cr),
#                   respect.unordered.factors = 'partition')
#       return(list(rf$variable.importance, rf$prediction.error))
#    })
#    return(perf)
# })
# save(perfPerPentamer, file = paste0("data/rdata/", tissue,
#                                     "/RFmodel_perfPerPentamer.RData"))
#####



# analyse ref base #####
load(paste0("data/rdata/", tissue, "/RFmodel_perfPerRef.RData"))
library(ggplot2)
library(reshape2)

dat = melt(table(ref))
p1 = ggplot(data = dat, aes(x = ref, y = value, fill = ref)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = "reference base", y = "n positions") +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12, color = "black"),
         legend.position = "none")

performance = sapply(perfPerRef, function(x){
   sapply(x, function(y){
      y[[2]]
   })
})
dat = melt(performance)
dat$Var2 = as.character(dat$Var2)
p2 = ggplot(data = dat, aes(x = Var2, y = value, fill = Var2)) +
   geom_boxplot() +
   theme_classic() +
   labs(x = "reference base", y = "prediction error") +
   theme(axis.title = element_text(size = 15),
         axis.text = element_text(size = 12, color = "black"),
         legend.position = "none")


imp = sapply(perfPerRef, function(x){
   temp = sapply(x, function(y){y[[1]]})
   rowMeans(temp)
})
imp = melt(imp)
imp$Var1 = as.character(imp$Var1)
colnames(imp) = c("predictor", "refBase", "importance")
imp$predictor = factor(imp$predictor,
                       levels = unique(imp$predictor))
p3 = ggplot(imp, aes(x = refBase, y = predictor)) + 
   geom_raster(aes(fill=importance)) + 
   labs(x = "reference base") +
   scale_fill_gradient(low="grey90", high="red")

imp = lapply(perfPerRef, function(x){
   temp = sapply(x, function(y){y[[1]]})
   return(temp)
})
imp = melt(imp)
colnames(imp) = c("predictor", "chr", "importance", "refBase")
p4 = ggplot(data=imp, aes(x = predictor, y =importance )) + 
   geom_boxplot(aes(color=refBase), outlier.size = 0.5) +
   labs(y = "permutation importance") +
   coord_flip()

library(gridExtra)
gs = list(p1, p2, p3, p4)
lay <- rbind(c(1,3,3,3,4,4,4),
             c(2,3,3,3,4,4,4),
             c(NA,3,3,3,4,4,4))
g = arrangeGrob(grobs = gs, layout_matrix = lay)
ggsave(paste0("fig/", tissue, "/RFmodel_perRef_", tissue, ".png"),
       plot=g, width=12, height=5)
#####


# analyse trimer #####
load(paste0("data/rdata/", tissue, "/RFmodel_perfPerTrimer.RData"))
library(ggplot2)
library(reshape2)

dat = melt(table(trimer))
p1 = ggplot(data = dat, aes(x = trimer, y = value, fill = trimer)) +
   geom_bar(stat = "identity") +
   theme_classic() +
   labs(x = "trimer", y = "n positions") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 10, color = "black", 
                                    angle = 90, vjust =0.5, hjust = 1),
         axis.text.y = element_text(size = 7, color = "black"),
         legend.position = "none") +
   coord_flip()

performance = sapply(perfPerTrimer, function(x){
   sapply(x, function(y){
      y[[2]]
   })
})
dat = melt(performance)
p2 = ggplot(data = dat, aes(x = Var2, y = value, 
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


imp = sapply(perfPerTrimer, function(x){
   temp = sapply(x, function(y){y[[1]]})
   rowMeans(temp)
})
imp = melt(imp)
colnames(imp) = c("predictor", "trimer", "importance")
p3 = ggplot(imp, aes(x = trimer, y = predictor)) + 
   geom_raster(aes(fill=importance)) + 
   scale_fill_gradient(low="grey90", high="red") +
   theme(axis.title = element_text(size = 15),
         axis.text.x = element_text(size = 6, color = "black",
                                  angle = 90, hjust=1)) 


imp = lapply(perfPerTrimer, function(x){
   temp = sapply(x, function(y){y[[1]]})
   return(temp)
})
imp = melt(imp)
colnames(imp) = c("predictor", "chr", "importance", "trimer")
p4 = ggplot(data=imp, aes(x = predictor, y =importance )) + 
   geom_boxplot(aes(color=trimer), outlier.size = 0.5) +
   labs(y = "permutation importance") +
   theme(legend.position='none') +
   coord_flip()

library(gridExtra)
gs = list(p1, p2, p3, p4)
lay <- rbind(c(1,3,3,3,4,4,4),
             c(2,3,3,3,4,4,4),
             c(NA,3,3,3,4,4,4))
g = arrangeGrob(grobs = gs, layout_matrix = lay)
# plot(g)
ggsave(paste0("fig/", tissue, "/RFmodel_perTrimer_", tissue, ".png"), 
       plot=g, width=12, height=5)
#####

# analyse pentamer #####
# load("data/rdata/luad/RFmodel_perfPerPentamer.RData")
# library(ggplot2)
# library(reshape2)
# dat = table(data$pentamer)
# dat = dat[dat >= 500]
# dat = melt(dat)
# p1 = ggplot(data = dat, aes(x = Var1, y = value, fill = Var1)) +
#    geom_bar(stat = "identity") +
#    theme_classic() +
#    labs(x = "pentamer", y = "n positions") +
#    theme(axis.title = element_text(size = 15),
#          axis.text.x = element_text(size = 10, color = "black", 
#                                     angle = 90, vjust = 0.5),
#          axis.text.y = element_text(size = 4, color = "black"),
#          legend.position = "none") +
#    coord_flip()
# 
# imp = sapply(perfPerPentamer, function(x){
#    temp = sapply(x, function(y){y[[1]]})
#    rowMeans(temp)
# })
# imp = melt(imp)
# colnames(imp) = c("predictor", "pentamer", "importance")
# p2 = ggplot(imp, aes(x = pentamer, y = predictor)) + 
#    geom_raster(aes(fill=importance)) + 
#    scale_fill_gradient(low="grey90", high="red") +
#    theme(axis.title = element_text(size = 15),
#          axis.text.x = element_text(size = 4, angle = 90, vjust = 0.5))
# 
# pE = sapply(perfPerPentamer, function(x){
#    sapply(x, function(y){
#       y[[2]]
#    })
# })
# dat = melt(pE)
# p3 = ggplot(data = dat, aes(x = Var2, y = value, 
#                        fill = Var2, color = Var2)) +
#    geom_boxplot() +
#    theme_classic() +
#    labs(x = "trimer", y = "prediction error") +
#    theme(axis.title = element_text(size = 15),
#          axis.text.x = element_text(size = 10, color = "black", 
#                                     angle = 90, vjust = 0.5),
#          axis.text.y = element_text(size = 4, color = "black"),
#          legend.position = "none") +
#    coord_flip()
# 
# library(gridExtra)
# gs = list(p1, p2, p3)
# lay <- rbind(c(1,3,2,2,2))
# g = arrangeGrob(grobs = gs, layout_matrix = lay)
# ggsave("fig/luad/RFmodel_perPentamer.png",
#        plot=g, width = 14, height = 6)
#####

