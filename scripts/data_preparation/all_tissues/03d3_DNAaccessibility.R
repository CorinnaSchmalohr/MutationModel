# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


# tissue-specific:
# I checked that there are no merged replicates.
meta = read.table(paste0("data/rawdata/DNAaccessibility/", tissue, "/metadata.tsv"), 
                  header = T, sep = "\t", stringsAsFactors = F)
meta = meta[meta$File.Status == "released" & 
               meta$Output.type %in% c("fold change over control", 
                                       "read-depth normalized signal", 
                                       "base overlap signal"),]
muts = read.table(paste0("data/procData/", tissue, "/MutsWithIDs.bed"))
print("local")
DNAaccessibility_tissue = sapply(unique(meta$Experiment.accession), function(e){
   print(e)
   sub = meta[meta$Experiment.accession == e,]
   if(nrow(sub) > 1){
      sub = sub[sub$Assembly == "hg19",]
   }
   f = sub$File.accession
   t = read.table(file = paste0("data/procData/", tissue, "/DNAaccessibility/",f,".out.bed"),
                  as.is = T,sep = "\t", na.strings=".")
   t = t[order(t$V1,t$V2),]
   t[is.na(t[,ncol(t)]),ncol(t)] = 0
   return(t[,ncol(t)])
})
DNAaccessibility_tissue = rowMeans(log(DNAaccessibility_tissue + 1), na.rm = T)

print("1kb")
DNAaccessibility_tissue_100bp = sapply(unique(meta$Experiment.accession), function(e){
   print(e)
   sub = meta[meta$Experiment.accession == e,]
   if(nrow(sub) > 1){
      sub = sub[sub$Assembly == "hg19",]
   }
   f = sub$File.accession
   t = read.table(file = paste0("data/procData/", tissue, 
                                "/DNAaccessibility/",f,"_100bp.out.bed"),
                  as.is = T,sep = "\t", na.strings=".")
   return(t[,ncol(t)])
})
DNAaccessibility_tissue_100bp = rowMeans(log(DNAaccessibility_tissue_100bp + 1), na.rm = T)


# DNAse from UCSC, from multiple cell lines combined
DNAaccessibility_UCSC = read.table(paste0("data/procData/", tissue,
                                          "/UCSC_tracks/wgEncodeRegDnaseClusteredV3.out"),
                                   as.is = T, na.strings = ".")[,4]
DNAaccessibility_UCSC[is.na(DNAaccessibility_UCSC)] = 0
save(DNAaccessibility_tissue, DNAaccessibility_UCSC, DNAaccessibility_tissue_100bp,
     file = paste0("data/rdata/", tissue, "/DNAaccessibility.RData"))

