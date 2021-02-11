# prepare chromosome lengths #####
fai = read.table("data/rawdata/GRCh37.p11.genome.fa.fai")
chrLengths = fai[,2]
names(chrLengths) = fai[,1]
chrLengths = chrLengths[paste0("chr", c(1:22))]
chrLengths = chrLengths[order(names(chrLengths))]
#####


# create Bins for chr1 and get number of mutations per bin #####
sizes = c("bins1Mb" = 1000000, "bins100kb" = 100000, 
          "bins10kb" = 10000)# "bins1kb" = 1000, "bins100bp" = 100
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
library(IRanges)
for(size in names(sizes)){
   print(size)
   bins = do.call(rbind,lapply(names(chrLengths), function(chr){
      start = as.integer(seq(1,chrLengths[chr],
                             by=sizes[size]))
      end = as.integer(c(start[2:length(start)]-1,
                         chrLengths[chr]))
      bins = data.frame(chr = chr, 
                        start = start, 
                        end = end)
   }))
   save(bins, file = paste0("data/rdata/", size, "_allChrs.RData"))
   
   # save bed file
   Bins = bins
   Bins$start = as.integer(Bins$start - 1)
   Bins$id = 1:nrow(Bins)
   write.table(Bins, file=paste0("data/procData/", size, "_allChrs.bed"), 
               quote=F, sep="\t", col.names=F, row.names=F)
   
   MutsPerBins = sapply(tissues, function(tissue){
      print(tissue)
      
      # get coverage per bin
      covered = read.table(paste0("data/procData/", tissue,
                                  "/covered_regions.bed"),as.is=T)
      colnames(covered) = c("chr", "start", "end")
      coverage = unlist(sapply(names(chrLengths), function(chr){
         subBins = bins[bins$chr == chr,]
         subcovered = covered[covered$chr == chr,2:3]
         binRanges = IRanges(subBins$start, subBins$end)
         covRanges = IRanges(subcovered$start, subcovered$end)
         overlapIndex = findOverlaps(binRanges, covRanges)
         overlapRanges = overlapsRanges(binRanges, covRanges, overlapIndex)
         hash = split(width(overlapRanges), queryHits(overlapIndex))
         hash = sapply(hash, sum)
         cov = rep(0, nrow(subBins))
         cov[as.numeric(names(hash))] = hash
         return(cov)
      }))
      
      
      # get nMuts per bin
      load(paste0("data/rdata/", tissue, "/Muts.RData"))
      nMuts = unlist(sapply(names(chrLengths), function(chr){
         subMuts = Muts[Muts$mutated == 1 & Muts$chr == chr,]
         subBins = bins[bins$chr == chr,]
         binRanges = IRanges(subBins$start, subBins$end)
         mutRanges = IRanges(subMuts$pos, subMuts$pos)
         overlapIndex = findOverlaps(binRanges, mutRanges)
         overlapRanges = overlapsRanges(binRanges, mutRanges, overlapIndex)
         hash = split(width(overlapRanges), queryHits(overlapIndex))
         hash = sapply(hash, sum)
         muts = rep(0, nrow(subBins))
         muts[as.numeric(names(hash))] = hash
         return(muts)
      }))
      nMuts = nMuts/coverage
      return(nMuts)
   })
   MutsPerBin = cbind(bins, MutsPerBins)
   save(MutsPerBin, file = paste0("data/rdata/MutsPerBin_", size, "_allChrs.RData"))
}
#####


# WGS #####
library(IRanges)
for(size in names(sizes)){
   load(paste0("data/rdata/", size, "_allChrs.RData"))
   WGSMutsPerBins = sapply(c("skin", "ovary", "kidney", "prostate"), function(tissue){
      print(tissue)
      # get nMuts per bin
      load(paste0("data/rdata/pancanWGS_icgc_muts_", 
                  tissue, ".RData"))
      nMuts = unlist(sapply(names(chrLengths), function(chr){
         subMuts = muts[muts$Chromosome == chr,]
         subBins = bins[bins$chr == chr,]
         binRanges = IRanges(subBins$start, subBins$end)
         mutRanges = IRanges(subMuts$Start_position, subMuts$End_position)
         overlapIndex = findOverlaps(binRanges, mutRanges)
         overlapRanges = overlapsRanges(binRanges, mutRanges, overlapIndex)
         hash = split(width(overlapRanges), queryHits(overlapIndex))
         hash = sapply(hash, sum)
         muts = rep(0, nrow(subBins))
         muts[as.numeric(names(hash))] = hash
         return(muts)
      }))
      return(nMuts)
   })
   WGSMutsPerBins = cbind(bins, WGSMutsPerBins)
   save(WGSMutsPerBins, file = paste0("data/rdata/MutsPerBin_WGS_", 
                                      size, "_allChrs.RData"))
}

#####


