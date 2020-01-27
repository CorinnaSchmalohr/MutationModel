# prepare chromosome lengths #####
fai = read.table("data/rawdata/GRCh37.p11.genome.fa.fai")
chrLengths = fai[,2]
names(chrLengths) = fai[,1]
chrLengths = chrLengths[paste0("chr", c(1:22))]
#####


# create Bins for chr1 and get number of mutations per bin #####
sizes = c("bins1Mb" = 1000000, "bins100kb" = 100000, "bins10kb" = 10000, 
          "bins1kb" = 1000, "bins100bp" = 100, "bins10bp" = 10)
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
library(IRanges)
for(size in names(sizes)){
   print(size)
   start = as.integer(seq(1,chrLengths["chr1"],by=sizes[size]))
   end = as.integer(c(start[2:length(start)]-1, chrLengths["chr1"]))
   bins = data.frame(chr = "chr1", 
                     start = start, 
                     end = end)
   MutsPerBins = sapply(tissues, function(tissue){
      print(tissue)
      
      # get coverage per bin
      covered = read.table(paste0("data/procData/", tissue,
                                  "/covered_regions_withsequences.bed"),as.is=T)
      colnames(covered) = c("chr", "start", "end", "sequence")
      covered = covered[covered$chr == "chr1",2:3]
      binRanges = IRanges(bins$start, bins$end)
      covRanges = IRanges(covered$start, covered$end)
      overlapIndex = findOverlaps(binRanges, covRanges)
      overlapRanges = overlapsRanges(binRanges, covRanges, overlapIndex)
      hash = split(width(overlapRanges), queryHits(overlapIndex))
      hash = sapply(hash, sum)
      coverage = rep(0, nrow(bins))
      coverage[as.numeric(names(hash))] = hash
      
      # get nMuts per bin
      load(paste0("data/rdata/", tissue, "/Muts.RData"))
      Muts = Muts[Muts$mutated == 1 & Muts$chr == "chr1",]
      binRanges = IRanges(bins$start, bins$end)
      mutRanges = IRanges(Muts$pos, Muts$pos)
      overlapIndex = findOverlaps(binRanges, mutRanges)
      overlapRanges = overlapsRanges(binRanges, mutRanges, overlapIndex)
      hash = split(width(overlapRanges), queryHits(overlapIndex))
      hash = sapply(hash, sum)
      muts = rep(0, nrow(bins))
      muts[as.numeric(names(hash))] = hash
      
      nMuts = muts/coverage
      return(nMuts)
   })
   MutsPerBin = cbind(bins, MutsPerBins)
   save(MutsPerBin, file = paste0("data/rdata/MutsPerBin_", size, ".RData"))
}
# MutsPerBin[is.na(MutsPerBin)] = 0
#####

# plot it #####
library(RColorBrewer)
cols = brewer.pal(length(tissues), "Set1")
names(cols) = tissues
for(size in names(sizes)[1:4]){
   print(size)
   load(paste0("data/rdata/MutsPerBin_", size, ".RData"))
   MutsPerBin = MutsPerBin[1:250,]
   # png(paste0("fig/nMutsPerMBbin_chr1_", size",.png"),
   #     width=1000,height=300)
   plot(NULL, xlim=c(0,max(MutsPerBin$end/1000000)),
        ylim=c(0,max(MutsPerBin[,4:10], na.rm=T)),  
        xlab = "Mb", ylab = "mutation density", main = size)
   for(i in tissues){
      lines(MutsPerBin$end/1000000,
            MutsPerBin[,i], lwd=2,
            type = "l", col = cols[i])
   }
   legend("topright", legend=tissues, col=cols,
          lty=1, lwd=2, ncol=2, bty = "n")
   # dev.off() 
}
#####

