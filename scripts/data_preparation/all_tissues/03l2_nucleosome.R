# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


NsomeGm12878 = read.table(paste0("data/procData/", tissue,
                          "/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.out.bed"))
NsomeK562 = read.table(paste0("data/procData/", tissue,
                                 "/UCSC_tracks/wgEncodeSydhNsomeK562Sig.out.bed"))
Nucleosome = cbind(NsomeGm12878 = NsomeGm12878[,5], NsomeK562 = NsomeK562[,5])
save(Nucleosome, file = paste0("data/rdata/", tissue, "/Nucleosome.RData"))