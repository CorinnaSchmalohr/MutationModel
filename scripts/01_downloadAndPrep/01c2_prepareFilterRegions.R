library(GenomicRanges)
load("data/processedData/chrLengths.RData")

temp = read.table("data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.bed")
temp = temp[temp$V1 %in% names(chrLengths),]

gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
# gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/hg19.repeatMasker.UCSC.RData")

temp = read.table("data/rawdata/UCSC_tracks/hg19.trf.bed.gz")
temp = temp[temp$V1 %in% names(chrLengths),]
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
# gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/hg19.trf.RData")


temp = read.table("data/rawdata/UCSC_tracks/wgEncodeDacMapabilityConsensusExcludable.bed.gz")
temp = temp[temp$V1 %in% names(chrLengths),]
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
# gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData")


# temp = read.table("data/processedData/MappabilityAlign24mer_regionsToKeep.bedGraph")
# gr = GRanges(seqnames=temp$V1, 
#              ranges=IRanges(start=temp$V2+1, end=temp$V3))
# gr = gr[seqnames(gr) %in% names(chrLengths)]
# save(gr, file = "data/predictors/MappabilityAlign24mer_regionsToKeep.RData")

temp = read.table("data/processedData/MappabilityAlign40mer_regionsToKeep.bedGraph")
temp = temp[temp$V1 %in% names(chrLengths),]
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
# gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/MappabilityAlign40mer_regionsToKeep.RData")

temp = read.table("data/processedData/MappabilityAlign100mer_regionsToKeep.bedGraph")
temp = temp[temp$V1 %in% names(chrLengths),]
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
# gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/MappabilityAlign100mer_regionsToKeep.RData")
