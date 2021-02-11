tissues = c( "colon", "luad", "breast", "skin",
             "ovary", "kidney", "prostate")
sizes = c("bins10kb", "bins100kb", "bins1Mb") # "bins100bp", "bins1kb",

for(tissue in tissues){
   for(size in sizes){
      print(size)
      NsomeGm12878 = read.table(paste0("data/procData/", tissue,
                                       "/UCSC_tracks/bins/", size,
                                       "_wgEncodeSydhNsomeGm12878Sig.out.bed"))
      NsomeGm12878 = NsomeGm12878[order(NsomeGm12878$V1, NsomeGm12878$V2),]
      NsomeK562 = read.table(paste0("data/procData/", tissue,
                                    "/UCSC_tracks/bins/", size,
                                    "_wgEncodeSydhNsomeK562Sig.out.bed"))
      NsomeK562 = NsomeK562[order(NsomeK562$V1, NsomeK562$V2),]
      Nucleosome = cbind(NsomeGm12878 = NsomeGm12878[,5],
                         NsomeK562 = NsomeK562[,5])
      save(Nucleosome, file = paste0("data/procData/", tissue,
                                     "/UCSC_tracks/bins/", size,
                                     "_Nucleosome.RData"))
   }
}
