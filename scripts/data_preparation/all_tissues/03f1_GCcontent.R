# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon",
            "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)

MutsBed = read.table(paste0("data/procData/", tissue, "/Muts.bed"),
                     col.names = c("chr", "start", "end"), as.is = T)
chrLength = read.table("data/rawdata/GRCh37.p11.genome.fa.fai",
                       sep = "\t")
chrL = chrLength[,2]
names(chrL) = chrLength[,1]
chrL = chrL[1:25]

dir.create(paste0("data/procData/", tissue, "/GCcontent"), showWarnings = F)

MutsBed5bp = MutsBed
MutsBed5bp$start = as.integer(pmax(MutsBed5bp$start - 2,0))
MutsBed5bp$end = as.integer(pmin(MutsBed5bp$end + 2, 
                                  chrL[MutsBed5bp$chr]))
write.table(MutsBed5bp, file = paste0("data/procData/", tissue, 
                                       "/GCcontent/MutsBed5bp.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)

MutsBed10bp = MutsBed
MutsBed10bp$start = as.integer(pmax(MutsBed10bp$start - 5,0))
MutsBed10bp$end = as.integer(pmin(MutsBed10bp$end + 5, 
                                  chrL[MutsBed10bp$chr]))
write.table(MutsBed10bp, file = paste0("data/procData/", tissue, 
                                       "/GCcontent/MutsBed10bp.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)

MutsBed100bp = MutsBed
MutsBed100bp$start = as.integer(pmax(MutsBed100bp$start - 50,0))
MutsBed100bp$end = as.integer(pmin(MutsBed100bp$end + 50,
                                   chrL[MutsBed100bp$chr]))
write.table(MutsBed100bp, file = paste0("data/procData/", tissue,
                                        "/GCcontent/MutsBed100bp.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)

MutsBed1kb = MutsBed
MutsBed1kb$start = as.integer(pmax(MutsBed1kb$start - 500,0))
MutsBed1kb$end = as.integer(pmin(MutsBed1kb$end + 500, 
                                 chrL[MutsBed1kb$chr]))
write.table(MutsBed1kb, file = paste0("data/procData/", tissue,
                                      "/GCcontent/MutsBed1kb.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)

MutsBed1Mb = MutsBed
MutsBed1Mb$start = as.integer(pmax(MutsBed1Mb$start - 500000,0))
MutsBed1Mb$end = as.integer(pmin(MutsBed1Mb$end + 500000, 
                                 chrL[MutsBed1Mb$chr]))
write.table(MutsBed1Mb, file = paste0("data/procData/", tissue,
                                      "/GCcontent/MutsBed1Mb.bed"),
            sep = "\t", col.names = F, row.names = F, quote = F)

