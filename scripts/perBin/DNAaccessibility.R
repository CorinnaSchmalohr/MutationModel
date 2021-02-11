tissues = c( "colon", "luad", "breast", "skin",
            "ovary", "kidney", "prostate")
sizes = c( "bins100bp", "bins1kb",# "bins10bp",
          "bins10kb", "bins100kb", "bins1Mb")
library(IRanges)
for(tissue in tissues){
   print(tissue)
   meta = read.table(paste0("data/rawdata/DNAaccessibility/", 
                            tissue, "/metadata.tsv"), 
                     header = T, sep = "\t", stringsAsFactors = F)
   meta = meta[meta$File.Status == "released" & 
                  meta$Output.type %in% c("fold change over control", 
                                          "read-depth normalized signal", 
                                          "base overlap signal"),]
   for(size in sizes){
      # complete bins:
      print(size)
      DNAaccessibility_tissue = sapply(unique(meta$Experiment.accession), function(e){
         sub = meta[meta$Experiment.accession == e,]
         if(nrow(sub) > 1){
            sub = sub[sub$Assembly == "hg19",]
         }
         f = sub$File.accession
         t = read.table(file = paste0("data/procData/bins/",size,"/",
                                      tissue, "/DNAaccessibility/",
                                      f,".out.bed"),
                        as.is = T,sep = "\t", na.strings="NAN")
         return(t[,5])
      })
      DNAaccessibility_tissue = rowMeans(scale(DNAaccessibility_tissue),
                                         na.rm = T)
      
      # DNAse from UCSC, from multiple cell lines combined
      print("UCSC")
      DNAaccessibility_UCSC = read.table(paste0("data/procData/bins/",size, "/", tissue,
                                                "/UCSC_tracks/wgEncodeRegDnaseClusteredV3.out"),
                                         as.is = T, na.strings = "NAN")[,5]
      DNAaccessibility = cbind(DNAaccessibility_tissue,
                               DNAaccessibility_UCSC)
      save(DNAaccessibility,
           file = paste0("data/procData/bins/", size, "/", tissue,
                         "/DNAaccessibility/DNAaccessibility.RData"))
      rm(DNAaccessibility, DNAaccessibility_UCSC, DNAaccessibility_tissue)
      # covered portions:
      print("covered")
      load(paste0("data/procData/bins/", size, "/", size, ".RData"))
      bin2Cov = read.table(paste0("data/procData/bins/", size, "/",
                                  size, "_", tissue, "_covered.bed"))
      colnames(bin2Cov) = c("chr", "start", "end", "id")
      bin2Cov$start = bin2Cov$start+1
      DNAaccessibility_tissue_covered = sapply(unique(meta$Experiment.accession), function(e){
         sub = meta[meta$Experiment.accession == e,]
         if(nrow(sub) > 1){
            sub = sub[sub$Assembly == "hg19",]
         }
         f = sub$File.accession
         t = read.table(paste0("data/procData/bins/",size,"/",
                               tissue, "/DNAaccessibility/covered_",
                               f,".out.bed"),
                        as.is = T,sep = "\t", na.strings="NAN")
         return(t[,5])
      })
      DNAaccessibility_tissue_covered = rowMeans(scale(DNAaccessibility_tissue_covered),
                                         na.rm = T)
      DNAaccessibility_UCSC_covered = read.table(
         paste0("data/procData/bins/",size,"/",
                tissue, "/UCSC_tracks/covered_wgEncodeRegDnaseClusteredV3.out"),
         as.is = T, na.strings = "NAN")[,5]
      # unify portions that belong in the same bin.
      print("unify")
      DNAaccessibility_covered = data.frame(bin2Cov,
                                            DNAaccessibility_tissue_covered,
                                            DNAaccessibility_UCSC_covered)
      binRanges = IRanges(bins$start, bins$end)
      covariateRanges = IRanges(bin2Cov$start, bin2Cov$end)
      overlapIndex = findOverlaps(binRanges, covariateRanges)
      overlapRanges = overlapsRanges(binRanges, covariateRanges, overlapIndex)
      vals = DNAaccessibility_covered[subjectHits(overlapIndex),5:6] * width(overlapRanges)/sum(width(overlapRanges))
      vals = apply(vals,2,function(x){
         sapply(split(x, queryHits(overlapIndex)), mean, na.rm = T)
      })
      res = matrix(NA, ncol=ncol(vals), nrow=nrow(bins),
                   dimnames=list(NULL, colnames(vals)))
      res[as.numeric(rownames(vals)),] = vals
      DNAaccessibility_covered = res
      save(DNAaccessibility_covered,
           file = paste0("data/procData/bins/", size, "/", tissue,
                         "/DNAaccessibility/DNAaccessibility_covered.RData"))
      rm(DNAaccessibility_covered, vals, binRanges, 
         DNAaccessibility_tissue_covered, DNAaccessibility_UCSC_covered)
   }
   gc()
}


