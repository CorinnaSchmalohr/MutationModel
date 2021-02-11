tissues = c("luad", "breast", "skin", "colon",
            "ovary", "kidney", "prostate")
bins = c( "bins100bp", "bins1kb", #"bins10bp",
          "bins10kb", "bins100kb", "bins1Mb")

for(tissue in tissues){
   print(tissue)
   meta = read.table(paste0("data/rawdata/DNAaccessibility/", 
                            tissue, "/metadata.tsv"), 
                     header = T, sep = "\t", stringsAsFactors = F)
   meta = meta[meta$File.Status == "released" & 
                  meta$Output.type %in% c("fold change over control", 
                                          "read-depth normalized signal", 
                                          "base overlap signal"),]
   for(bin in bins){
      print(bin)
      DNAaccessibility_tissue = sapply(unique(meta$Experiment.accession), function(e){
         cat(e, ' ')
         sub = meta[meta$Experiment.accession == e,]
         if(nrow(sub) > 1){
            sub = sub[sub$Assembly == "hg19",]
         }
         f = sub$File.accession
         t = read.table(file = paste0("data/procData/bins/",bin,"/",
                                      tissue, "/DNAaccessibility/",
                                      f,".out.bed"),
                        as.is = T,sep = "\t", na.strings="NAN")
         return(t[,5])
      })
      DNAaccessibility_tissue = rowMeans(scale(DNAaccessibility_tissue),
                                         na.rm = T)
      
      # DNAse from UCSC, from multiple cell lines combined
      DNAaccessibility_UCSC = read.table(paste0("data/procData/bins/",bin, "/", tissue,
                                                "/UCSC_tracks/wgEncodeRegDnaseClusteredV3.out"),
                                         as.is = T, na.strings = ".")[,5]
      DNAaccessibility = data.frame(DNAaccessibility_tissue,DNAaccessibility_UCSC)
      save(DNAaccessibility,
           file = paste0("data/procData/bins/", bin, "/", tissue,
                         "/DNAaccessibility/DNAaccessibility.RData"))
   }
   DNAaccessibility_tissue = sapply(unique(meta$Experiment.accession), function(e){
      cat(e,' ')
      sub = meta[meta$Experiment.accession == e,]
      if(nrow(sub) > 1){
         sub = sub[sub$Assembly == "hg19",]
      }
      f = sub$File.accession
      t = read.table(file = paste0("data/procData/",tissue, 
                                   "/DNAaccessibility/covered_",
                                   f,".out.bed"),
                     as.is = T,sep = "\t", na.strings="NAN")
      return(t[,5])
   })
   DNAaccessibility_tissue = rowMeans(scale(DNAaccessibility_tissue),
                                      na.rm = T)
   # DNAse from UCSC, from multiple cell lines combined
   DNAaccessibility_UCSC = read.table(paste0("data/procData/", tissue,
                                             "/UCSC_tracks/covered_wgEncodeRegDnaseClusteredV3.out"),
                                      as.is = T, na.strings = ".")[,5]
   DNAaccessibility = data.frame(DNAaccessibility_tissue,DNAaccessibility_UCSC)
   save(DNAaccessibility,
        file = paste0("data/procData/",tissue, 
                      "/DNAaccessibility/covered_DNAaccessibility.RData"))
}


