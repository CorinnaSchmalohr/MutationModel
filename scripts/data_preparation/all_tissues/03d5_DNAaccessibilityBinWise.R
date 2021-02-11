tissues = c(  "luad","colon", "breast", "skin",
             "ovary", "kidney", "prostate")
sizes = c("bins10kb", "bins100kb", "bins1Mb") # "bins100bp", "bins1kb",
for(tissue in tissues){
   print(tissue)
   meta = read.table(paste0("data/rawdata/DNAaccessibility/", 
                            tissue, "/metadata.tsv"), 
                     header = T, sep = "\t", stringsAsFactors = F)
   meta = meta[meta$File.Status == "released" & 
                  meta$Output.type %in% c("fold change over control", 
                                          "read-depth normalized signal", 
                                          "base overlap signal"),]
   for(size in sizes){
      print(size)
      DNAaccessibility_tissue = sapply(unique(meta$Experiment.accession), function(e){
         sub = meta[meta$Experiment.accession == e,]
         if(nrow(sub) > 1){
            sub = sub[sub$Assembly == "hg19",]
         }
         f = sub$File.accession
         t = read.table(file = paste0("data/procData/", tissue, 
                                      "/DNAaccessibility/bins/",
                                      size, "_", f,".out.bed"),
                        as.is = T,sep = "\t", na.strings="NAN")
         t = t[order(t$V1, t$V2),5]
         t[is.na(t)] = 0
         return(t)
      })
      DNAaccessibility_tissue = rowMeans(scale(DNAaccessibility_tissue),
                                         na.rm = T)
      
      # DNAse from UCSC, from multiple cell lines combined
      print("UCSC")
      DNAaccessibility_UCSC = read.table(paste0("data/procData/", tissue, 
                                                "/UCSC_tracks/bins/",
                                                size, "_wgEncodeRegDnaseClusteredV3.out.bed"),
                                         as.is = T, na.strings = "NAN")
      DNAaccessibility_UCSC = DNAaccessibility_UCSC[order(DNAaccessibility_UCSC$V1, 
                                                          DNAaccessibility_UCSC$V2),5]
      DNAaccessibility_UCSC[is.na(DNAaccessibility_UCSC)] = 0
      DNAaccessibility = cbind(DNAaccessibility_tissue,
                               DNAaccessibility_UCSC)
      save(DNAaccessibility,
           file = paste0("data/procData/", tissue,
                         "/DNAaccessibility/bins/", size,
                         "DNAaccessibility.RData"))
   }
}


