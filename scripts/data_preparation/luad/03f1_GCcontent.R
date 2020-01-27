LungMutsBed = read.table("data/procData/luad/Muts.bed",col.names = c("chr", "start", "end"), as.is = T)
chrLength = read.table("data/rawdata/GRCh37.p11.genome.fa.fai", sep = "\t")
chrL = chrLength[,2]
names(chrL) = chrLength[,1]
chrL = chrL[1:25]

dir.create("data/procData/luad/GCcontent")
LungMutsBed10bp = LungMutsBed
LungMutsBed10bp$start = as.integer(pmax(LungMutsBed10bp$start - 5,0))
LungMutsBed10bp$end = as.integer(pmin(LungMutsBed10bp$end + 5, chrL[LungMutsBed10bp$chr]))
write.table(LungMutsBed10bp, file = "data/procData/luad/GCcontent/MutsBed10bp.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)

LungMutsBed100bp = LungMutsBed
LungMutsBed100bp$start = as.integer(pmax(LungMutsBed100bp$start - 50,0))
LungMutsBed100bp$end = as.integer(pmin(LungMutsBed100bp$end + 50, chrL[LungMutsBed100bp$chr]))
write.table(LungMutsBed100bp, file = "data/procData/luad/GCcontent/MutsBed100bp.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)

LungMutsBed1kb = LungMutsBed
LungMutsBed1kb$start = as.integer(pmax(LungMutsBed1kb$start - 500,0))
LungMutsBed1kb$end = as.integer(pmin(LungMutsBed1kb$end + 500, chrL[LungMutsBed1kb$chr]))
write.table(LungMutsBed1kb, file = "data/procData/luad/GCcontent/MutsBed1kb.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)

LungMutsBed1Mb = LungMutsBed
LungMutsBed1Mb$start = as.integer(pmax(LungMutsBed1Mb$start - 500,0))
LungMutsBed1Mb$end = as.integer(pmin(LungMutsBed1Mb$end + 500, chrL[LungMutsBed1Mb$chr]))
write.table(LungMutsBed1Mb, file = "data/procData/luad/GCcontent/MutsBed1Mb.bed",
            sep = "\t", col.names = F, row.names = F, quote = F)

