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
   save(bins, file = paste0("data/procData/bins/", size, "/", size, ".RData"))
   Bins = bins
   Bins$start = as.integer(Bins$start - 1)
   Bins$id = 1:nrow(Bins)
   write.table(Bins, file=paste0("data/procData/bins/", size, "/", size, ".bed"), 
               quote=F, sep="\t", col.names=F, row.names=F)

   MutsPerBins = sapply(tissues, function(tissue){
      print(tissue)
      
      # get coverage per bin
      covered = read.table(paste0("data/procData/", tissue,
                                  "/covered_regions.bed"),as.is=T)
      colnames(covered) = c("chr", "start", "end")
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


# WGS #####
library(IRanges)
for(size in names(sizes)){
   load(paste0("data/procData/bins/", size, "/", size, ".RData"))
   WGSMutsPerBins = sapply(c("skin", "ovary", "kidney", "prostate"), function(tissue){
      print(tissue)
      # get nMuts per bin
      load(paste0("data/rdata/pancanWGS_icgc_muts_", 
                  tissue, ".RData"))
      muts = muts[muts$Chromosome == "chr1",]
      binRanges = IRanges(bins$start, bins$end)
      mutRanges = IRanges(muts$Start_position, muts$End_position)
      overlapIndex = findOverlaps(binRanges, mutRanges)
      overlapRanges = overlapsRanges(binRanges, mutRanges, overlapIndex)
      hash = split(width(overlapRanges), queryHits(overlapIndex))
      hash = sapply(hash, sum)
      nMuts = rep(0, nrow(bins))
      nMuts[as.numeric(names(hash))] = hash
      return(nMuts)
   })
   WGSMutsPerBins = cbind(bins, WGSMutsPerBins)
   save(WGSMutsPerBins, file = paste0("data/rdata/MutsPerBin_WGS_", size, ".RData"))
}

#####


