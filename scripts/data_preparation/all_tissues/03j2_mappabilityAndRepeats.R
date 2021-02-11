# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


# mappability 
mappability = sapply(c("24mer", "40mer", "100mer"),function(x){
   t = read.table(paste0("data/procData/", tissue, 
                         "/UCSC_tracks/wgEncodeCrgMapabilityAlign",
                         x,".out.bed"))
   t = t[order(t$V1, t$V2),]
   return(t[,5])
})
colnames(mappability) = paste0("mappability_", c("24mer", "40mer", "100mer"))

# TandemRepeatsFinder
Trf = read.table(paste0("data/procData/", tissue,
                        "/UCSC_tracks/hg19.TandemRepeatsFinder.UCSC.out"),
                 as.is = T, na.strings = ".")
Trf = Trf[order(Trf$V1, Trf$V2),]
Trf = as.integer(!is.na(Trf[,4]))

# repeatMasker
repeatMasker = read.table(paste0("data/procData/", tissue,
                                 "/UCSC_tracks/hg19.repeatMasker.UCSC.out"),
                          as.is = T, na.strings = ".")
repeatMasker = repeatMasker[order(repeatMasker$V1, repeatMasker$V2),]
repeatMasker = as.integer(is.na(repeatMasker[,4]))

# save
save(mappability, Trf,repeatMasker, 
     file = paste0("data/rdata/", tissue, "/mappability_repeats.RData"))