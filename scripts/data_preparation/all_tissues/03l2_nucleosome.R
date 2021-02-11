# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


NsomeGm12878 = read.table(paste0("data/procData/", tissue,
                          "/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.out.bed"))
NsomeGm12878 = NsomeGm12878[(order(NsomeGm12878$V1, NsomeGm12878$V2)),]
NsomeK562 = read.table(paste0("data/procData/", tissue,
                                 "/UCSC_tracks/wgEncodeSydhNsomeK562Sig.out.bed"))
NsomeK562 = NsomeK562[(order(NsomeK562$V1, NsomeK562$V2)),]

Nucleosome = cbind(NsomeGm12878 = NsomeGm12878[,5], NsomeK562 = NsomeK562[,5])
save(Nucleosome, file = paste0("data/rdata/", tissue, "/Nucleosome.RData"))

NsomeGm12878_100bp = read.table(paste0("data/procData/", tissue,
                                 "/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig_100bp.out.bed"))
NsomeGm12878_100bp = NsomeGm12878_100bp[(order(NsomeGm12878_100bp$V1, NsomeGm12878_100bp$V2)),]
NsomeK562_100bp = read.table(paste0("data/procData/", tissue,
                              "/UCSC_tracks/wgEncodeSydhNsomeK562Sig_100bp.out.bed"))
NsomeK562_100bp = NsomeK562_100bp[(order(NsomeK562_100bp$V1, NsomeK562_100bp$V2)),]
Nucleosome_100bp = cbind(NsomeGm12878_100bp = NsomeGm12878_100bp[,5],
                         NsomeK562_100bp = NsomeK562_100bp[,5])
save(Nucleosome_100bp, file = paste0("data/rdata/", tissue, "/Nucleosome_100bp.RData"))