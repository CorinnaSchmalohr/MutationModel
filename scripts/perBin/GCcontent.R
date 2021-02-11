tissues = c("luad", "breast", "skin", "colon",
            "ovary", "kidney", "prostate")
sizes = c("bins1Mb" = 1000000, "bins100kb" = 100000,
          "bins10kb" = 10000, "bins1kb" = 1000,
          "bins100bp" = 100, "bins10bp" = 10)
library(IRanges)
for(bin in names(sizes)){
   print(bin)
   load(paste0("data/procData/bins/", bin, "/", bin, ".RData"))
   GCcontent = read.table(paste0("data/procData/bins/",bin,
                                 "/GCcontent.out"))
   missing = which(!1:nrow(bins) %in% GCcontent[,2])
   if (length(missing) > 0) {
      missing = cbind(V1 = NA, V2 = missing)
      GCcontent = rbind(GCcontent, missing)
      GCcontent = GCcontent[order(GCcontent[,2]),]
   }
   GCcontent = GCcontent[,1]/sizes[bin]
   save(GCcontent, file = paste0("data/procData/bins/",bin,
                                 "/GCcontent.RData"))
   
   load(paste0("data/procData/bins/", bin, "/", bin, ".RData"))
   for (tissue in tissues){
      print(tissue)
      GCcontent = read.table(paste0("data/procData/bins/", bin, "/",
                                    tissue, "_covered_GCcontent.out"))
      bin2Cov = read.table(paste0("data/procData/bins/", bin, "/",
                                  bin, "_", tissue, "_covered.bed"))
      colnames(bin2Cov) = c("chr", "start", "end", "id")
      bin2Cov$start = bin2Cov$start+1
      missing = which(!1:nrow(bin2Cov) %in% GCcontent[,2])
      if (length(missing) > 0) {
         missing = cbind(V1 = NA, V2 = missing)
         GCcontent = rbind(GCcontent, missing)
         GCcontent = GCcontent[order(GCcontent[,2]),]
      }
      binRanges = IRanges(bins$start, bins$end)
      covariateRanges = IRanges(bin2Cov$start, bin2Cov$end)
      overlapIndex = findOverlaps(binRanges, covariateRanges)
      overlapRanges = overlapsRanges(binRanges, covariateRanges, overlapIndex)
      vals = GCcontent[subjectHits(overlapIndex),1] / width(overlapRanges)
      vals = sapply(split(vals, queryHits(overlapIndex)), mean, na.rm = T)
      res = rep(NA, nrow(bins))
      res[as.numeric(names(vals))] = vals
      GCcontent_covered = res
      save(GCcontent_covered, file = paste0("data/procData/bins/", bin, "/",
                                    tissue, "_covered_GCcontent.RData"))
   }
}
#####

