# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)

# ENCODE files
files = list.files(paste0("data/procData/", tissue, "/methylation/"),
                   pattern = "out")
ids = unique(do.call(rbind,strsplit(files,split = "_"))[,1])
meta = read.table(paste0("data/rawdata/methylation/",
                         tissue, "/metadata.tsv"), 
                  sep = "\t", header = T, as.is = T)
rownames(meta) = meta$File.accession
meta = meta[ids,]
meth = sapply(unique(meta$Experiment.accession), function(e){
   f = meta$File.accession[meta$Experiment.accession == e]
   print(f)
   local = sapply(f,function(i){
      tab = read.table(paste0("data/procData/", tissue,
                              "/methylation/",i,"_local.out"),
                 sep = "\t", as.is = T, na.strings = ".")
      tab$V6[tab$V4 < 10]  = NA
      return(tab$V6)
   })
   kb = sapply(f,function(i){
      tab = read.table(paste0("data/procData/", tissue,
                              "/methylation/",i,"_1kbp.out"),
                       sep = "\t", as.is = T, na.strings = ".")
      tab$V6[tab$V4 < 10]  = NA
      return(tab$V6)
   })
   Mb = sapply(f,function(i){
      tab = read.table(paste0("data/procData/", tissue,
                              "/methylation/",i,"_1Mbp.out"),
                       sep = "\t", as.is = T, na.strings = ".")
      tab$V6[tab$V4 < 10]  = NA
      return(tab$V6)
   })

   res = data.frame(methylation_local = rowMeans(local, na.rm = T),
                    methylation_1kb = rowMeans(kb, na.rm = T),
                    methylation_1Mb = rowMeans(Mb, na.rm = T))
   return(res)
})


if(ncol(meth)>1){
   dir.create(paste0("fig/", tissue, "/DNAmethylation"), showWarnings=F)
   png(paste0("fig/", tissue, "/DNAmethylation/methylationCrossplot_local.png"),
       height = 800, width = 800)
   temp = as.data.frame(meth[1,])
   colnames(temp) = meta$Biosample.term.name[(2*(1:(nrow(meta)/2)))]
   dev.off()
   png(paste0("fig/", tissue, "/DNAmethylation/methylationCrossplot_1kbp.png"),
       height = 800, width = 800)
   temp = as.data.frame(meth[2,])
   colnames(temp) = meta$Biosample.term.name[(2*(1:(nrow(meta)/2)))]
   plot(temp, pch = 19, cex = 0.2, main = "1kbp")
   dev.off()
   png(paste0("fig/", tissue, "/DNAmethylation/methylationCrossplot_1Mbp.png"),
       height = 800, width = 800)
   temp = as.data.frame(meth[3,])
   colnames(temp) = meta$Biosample.term.name[(2*(1:(nrow(meta)/2)))]
   plot(temp, pch = 19, cex = 0.2, main = "1Mbp")
   dev.off()
} else{
   # png(paste0("fig/", tissue, "/DNAmethylation/methylationCrossplot.png"),
   #     height = 800, width = 800)
   # colnames(temp) = meta$Biosample.term.name[(2*(1:7))]
   # plot(temp, pch = 19, cex = 0.2, main = "1Mbp")
   # dev.off()
}

methylation = apply(meth,1,function(x){
   rowMeans(as.data.frame(x), na.rm = T)
})
plot(methylation)
methylation[is.na(methylation)] = 0

colnames(methylation) = paste0("methylation_",c("local",  "1kb","1Mb"))
save(methylation, file = paste0("data/rdata/", tissue,"/methylation.RData"))


tissue2File = c("luad" = "Lung", "breast" = "Breast", "skin" = "Skin", 
            "colon" = "Colon", "kidney" = "Kidney", "prostate" = "Prostate")
for (tissue in names(tissue2File)){
   meth = sapply(c(".out.bed", "_1kbp.out.bed",
                   "_1Mbp.out.bed"), function(i){
      tab = read.table(paste0("data/procData/", tissue,
                              "/methylation_methbank/",
                              tissue2File[tissue],i),
                       sep = "\t", as.is = T, na.strings = ".")
      res = tab$V4
      res[is.na(res)] = 0
      return(res)
   })
   colnames(meth) = paste("methbank",
                          c("local", "1kbp", "1Mbp"), sep = "_")
   save(meth, file = paste0("data/rdata/", tissue, "/methbank.RData"))
}
