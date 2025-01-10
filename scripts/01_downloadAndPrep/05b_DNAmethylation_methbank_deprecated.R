library(GenomicRanges)
dir.create("data/predictors/methylation")
files = list.files("data/rawdata/methylation_methbank/",
                   pattern=".*_liftover.bed")
tissueNames = tolower(sapply(files, function(x){
   strsplit(x,split="_", fixed = T)[[1]][1]
}))
files = setNames(files, tissueNames)
dumpVar = sapply(names(files), function(tissue){
   print(tissue)
   file = files[tissue]
   temp = read.table(paste0("data/rawdata/methylation_methbank/",
                            file))
   gr = GRanges(seqnames=temp$V1,
                ranges=IRanges(start=temp$V2+1, end=temp$V3), 
                val = temp$V5)
   gr = sort(gr)
   save(gr, file = paste0("data/predictors/methylation/", 
                          tissue, ".RData"))
})
