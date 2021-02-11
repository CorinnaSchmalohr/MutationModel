tissues = c( "colon", "luad", "breast", "skin",
             "ovary", "kidney", "prostate")
sizes = c( "bins100bp", "bins1kb",# "bins10bp",
           "bins10kb", "bins100kb", "bins1Mb")
library(IRanges)
for(tissue in tissues){
   for(size in sizes){
      # complete bins:
      print(size)
      NsomeGm12878 = read.table(paste0("data/procData/bins/", size,
                                       "/",tissue,
                                       "/wgEncodeSydhNsomeGm12878Sig.out.bed"))
      NsomeK562 = read.table(paste0("data/procData/bins/", size,
                                    "/",tissue,
                                    "/wgEncodeSydhNsomeK562Sig.out.bed"))
      Nucleosome = cbind(NsomeGm12878 = NsomeGm12878[,5],
                         NsomeK562 = NsomeK562[,5])
      save(Nucleosome, file = paste0("data/procData/bins/", size, "/", tissue,
                                     "/Nucleosome.RData"))
      
      # covered portions:
      print("covered")
      load(paste0("data/procData/bins/", size, "/", size, ".RData"))
      bin2Cov = read.table(paste0("data/procData/bins/", size, "/",
                                  size, "_", tissue, "_covered.bed"))
      colnames(bin2Cov) = c("chr", "start", "end", "id")
      bin2Cov$start = bin2Cov$start+1
      NsomeGm12878 = read.table(paste0("data/procData/bins/", size,
                                       "/",tissue,
                                       "/wgEncodeSydhNsomeGm12878Sig_covered.out.bed"))
      NsomeK562 = read.table(paste0("data/procData/bins/", size,
                                    "/",tissue,
                                    "/wgEncodeSydhNsomeK562Sig_covered.out.bed"))
      Nucleosome_covered = cbind(NsomeGm12878 = NsomeGm12878[,5],
                                     NsomeK562 = NsomeK562[,5])
      
      # unify portions that belong in the same bin.
      print("unify")
      binRanges = IRanges(bins$start, bins$end)
      covariateRanges = IRanges(bin2Cov$start, bin2Cov$end)
      overlapIndex = findOverlaps(binRanges, covariateRanges)
      overlapRanges = overlapsRanges(binRanges, covariateRanges, overlapIndex)
      vals = Nucleosome_covered[subjectHits(overlapIndex),] * width(overlapRanges)/sum(width(overlapRanges))
      vals = apply(vals,2,function(x){
         sapply(split(x, queryHits(overlapIndex)), mean, na.rm = T)
      })
      res = matrix(NA, ncol=ncol(vals), nrow=nrow(bins),
                   dimnames=list(NULL, colnames(vals)))
      res[as.numeric(rownames(vals)),] = vals
      Nucleosome_covered = res
      save(Nucleosome_covered,
           file = paste0("data/procData/bins/", size, "/", tissue,
                         "/Nucleosome_covered.RData"))
   }
   gc()
}
