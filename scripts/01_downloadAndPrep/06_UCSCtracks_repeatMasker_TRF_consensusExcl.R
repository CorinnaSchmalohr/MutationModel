library(GenomicRanges)
load("data/procData/chrLengths.RData")

temp = read.table("data/rawdata/UCSC_tracks/hg19.repeatMasker.UCSC.bed")
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/hg19.repeatMasker.UCSC.RData")

temp = read.table("data/rawdata/UCSC_tracks/hg19.trf.bed.gz")
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/hg19.trf.RData")


temp = read.table("data/rawdata/UCSC_tracks/wgEncodeDacMapabilityConsensusExcludable.bed.gz")
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3))
gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData")



temp = read.table("data/rawdata/UCSC_tracks/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz")
nCells = sapply(strsplit(temp$V6, ","), length)
gr = GRanges(seqnames=temp$V1, 
             ranges=IRanges(start=temp$V2+1, end=temp$V3),
             TF = temp$V4, nCells = nCells)
gr = gr[seqnames(gr) %in% names(chrLengths)]
save(gr, file = "data/predictors/wgEncodeRegTfbsClusteredWithCellsV3.RData")
ETS = c("ELF", "ELF1", "ELF2", "NERF", "ELF4", "MEF", 
        "ELG", "GABPA", "GABP", 
        "ERG", "FLI1", "FEV",
        "ERF", "PE-2", "ETV3", "PE1", 
        "ESE", "ELF3", "ESE1", "ESX", "ELF5", "ESE2", "ESE3", "EHF",
        "ETS", "ETS1", "ETS2", 
        "PDEF", "SPDEF", "PSE", 
        "PEA3", "ETV4", "E1AF", "ETV5", "ERM", "ETV1", "ER81",
        "ER71", "ETV2",
        "SPI", "SPI1", "PU.1", "SPIB", "SPIC",
        "TCF", "ELK1", "ELK4", "SAP1", "SAP-1", "ELK3", "NET", "SAP2",
        "TEL", "ETV6", "ETV7", "TEL2") 
ETS_temp = temp[temp$V4 %in% ETS,]
gr = GRanges(seqnames=ETS_temp$V1, 
             ranges=IRanges(start=ETS_temp$V2+1, end=ETS_temp$V3))
save(gr, file = "data/predictors/wgEncodeRegTfbsClusteredWithCellsV3_ETS.RData")


# UCSC histones was removed