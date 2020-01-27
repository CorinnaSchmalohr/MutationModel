files = list.files("data/procData/luad/methylation/", pattern = "out")
ids = unique(do.call(rbind,strsplit(files,split = "[.]"))[,1])
meta = read.table("data/rawdata/methylation/luad/metadata.tsv", 
                  sep = "\t", header = T, as.is = T)
info = meta[meta$File.accession %in% ids,
            c("File.accession", "Experiment.accession","Biosample.term.name",
              "Biosample.type", "Biological.replicate.s.")]
meth = sapply(unique(info$Experiment.accession), function(e){
   f = info$File.accession[info$Experiment.accession == e]
   print(f)
   r1 = read.table(paste0("data/procData/luad/methylation/",f[1],".bed.gz_local.out"),
                   sep = "\t", as.is = T, na.strings = ".")
   r2 = read.table(paste0("data/procData/luad/methylation/",f[2],".bed.gz_local.out"),
                   sep = "\t", as.is = T, na.strings = ".")
   r1_1kb = read.table(paste0("data/procData/luad/methylation/",f[1],".bed.gz_1kbp.out"),
                       sep = "\t", as.is = T, na.strings = ".")
   r2_1kb = read.table(paste0("data/procData/luad/methylation/",f[2],".bed.gz_1kbp.out"),
                       sep = "\t", as.is = T, na.strings = ".")
   r1_1Mb = read.table(paste0("data/procData/luad/methylation/",f[1],".bed.gz_1Mbp.out"),
                       sep = "\t", as.is = T, na.strings = ".")
   r2_1Mb = read.table(paste0("data/procData/luad/methylation/",f[2],".bed.gz_1Mbp.out"),
                       sep = "\t", as.is = T, na.strings = ".")
   r1$V6[r1$V4 < 10]  = NA
   r2$V6[r2$V4 < 10]  = NA
   r1_1kb$V6[r1_1kb$V4 < 10]  = NA
   r2_1kb$V6[r2_1kb$V4 < 10]  = NA
   r1_1Mb$V6[r1_1Mb$V4 < 10]  = NA
   r2_1Mb$V6[r2_1Mb$V4 < 10]  = NA
   res = data.frame(perc_local = rowMeans(cbind(r1[,6], r2[,6]), na.rm = T),
                    perc_1kb = rowMeans(cbind(r1_1kb[,6], r2_1kb[,6]), na.rm = T),
                    perc_1Mb = rowMeans(cbind(r1_1Mb[,6], r2_1Mb[,6]), na.rm = T))
   return(res)
})
methylation = do.call(cbind,meth)
colnames(methylation) = paste(rep(info$Experiment.accession[(2*(1:7))],each = 3), 
                              c("local", "1kb", "1Mb"), sep = "_")

dir.create("fig/luad")
dir.create("fig/luad/DNAmethylation")
png("fig/luad/DNAmethylation/methylationCrossplot_local.png",
    height = 800, width = 800)
temp = as.data.frame(meth[1,])
colnames(temp) = info$Biosample.term.name[(2*(1:7))]
plot(temp, pch = 19, cex = 0.2, main = "local")
dev.off()
png("fig/luad/DNAmethylation/methylationCrossplot_1kbp.png",
    height = 800, width = 800)
temp = as.data.frame(meth[3,])
colnames(temp) = info$Biosample.term.name[(2*(1:7))]
plot(temp, pch = 19, cex = 0.2, main = "1kbp")
dev.off()
png("fig/luad/DNAmethylation/methylationCrossplot_1Mbp.png",
    height = 800, width = 800)
temp = as.data.frame(meth[5,])
colnames(temp) = info$Biosample.term.name[(2*(1:7))]
plot(temp, pch = 19, cex = 0.2, main = "1Mbp")
dev.off()

load("data/rdata/luad/Muts.RData")

library(corrplot)
corrplot(cor(cbind(methylation, Muts$mutated), use = "pair"))
methylation = apply(meth,1,function(x){
   rowMeans(as.data.frame(x), na.rm = T)
})
methylation[is.na(methylation)] = 0

colnames(methylation) = paste0("methylation_",c("local",  "1kb","1Mb"))
save(methylation, file = "data/rdata/luad/methylation.RData")
