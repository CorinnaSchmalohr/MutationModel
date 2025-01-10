library(GenomicRanges)
load("data/processedData/chrLengths.RData")
t = read.table("data/rawdata/gap.txt.gz")
# there is an error for GRCh37: telomeres are missing for chr17. Add them manually:
telomeres = t[t$V8 == "telomere",]
telomereGR = GRanges(seqnames=c(telomeres$V2, "chr17", "chr17"), 
                     ranges=IRanges(start=c(telomeres$V3+1,
                                            1,
                                            max(t$V4[t$V2 == "chr17"])+1), 
                                    end=c(telomeres$V4,
                                          min(t$V3[t$V2 == "chr17"]),
                                          chrLengths["chr17"])))
centromeres = t[t$V8 == "centromere",]
centromereGR = GRanges(seqnames=centromeres$V2, 
                       ranges=IRanges(start=centromeres$V3, end=centromeres$V4))
save(telomereGR, file="data/predictors/telomeres.RData")
save(centromereGR, file="data/predictors/centromeres.RData")
