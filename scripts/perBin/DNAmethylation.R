tissues = c( "colon", "luad", "breast", "skin", # "ovary",
              "kidney", "prostate")
sizes = c( "bins100bp", "bins1kb",# "bins10bp",
           "bins10kb", "bins100kb", "bins1Mb")
library(IRanges)
for(tissue in tissues){
   print(tissue)
   for(size in sizes){
      # complete bins:
      print(size)
      methylation = read.table(paste0("data/procData/bins/", size, "/",
                              tissue,"/methylation.out.bed"),
                              as.is = T,sep = "\t", na.strings="NAN")[,5]
      save(methylation,
           file = paste0("data/procData/bins/", size, "/", tissue,
                         "/methylation/methylation.RData"))
      
      # covered portions:
      print("covered")
      methylation_covered = read.table(paste0("data/procData/bins/", size, "/",
                                              tissue,"/methylation_covered.out.bed"),
                                       as.is = T,sep = "\t", na.strings="NAN")[,5]
      load(paste0("data/procData/bins/", size, "/", size, ".RData"))
      bin2Cov = read.table(paste0("data/procData/bins/", size, "/",
                                  size, "_", tissue, "_covered.bed"))
      colnames(bin2Cov) = c("chr", "start", "end", "id")
      bin2Cov$start = bin2Cov$start+1
      binRanges = IRanges(bins$start, bins$end)
      covariateRanges = IRanges(bin2Cov$start, bin2Cov$end)
      overlapIndex = findOverlaps(binRanges, covariateRanges)
      overlapRanges = overlapsRanges(binRanges, covariateRanges, overlapIndex)
      vals = methylation_covered[subjectHits(overlapIndex)] / width(overlapRanges)
      vals = sapply(split(vals, queryHits(overlapIndex)), mean, na.rm = T)
      res = rep(NA, nrow(bins))
      res[as.numeric(names(vals))] = vals
      methylation_covered = res
      save(methylation_covered,
           file = paste0("data/procData/bins/", size, "/", tissue,
                         "/methylation/methylation_covered.RData"))
   }
   gc()
}


