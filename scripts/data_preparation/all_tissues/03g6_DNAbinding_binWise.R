tissues = c("luad", "breast", "skin", "colon",
            "ovary", "kidney", "prostate")
sizes = c("bins10kb", "bins100kb", "bins1Mb") # "bins100bp", "bins1kb",

marks = c("H3k27ac", "H3k4me1", "H3k4me3")
lines = c("Gm12878", "H1hesc", "Hsmm", "Huvec", 
          "K562", "Nhek", "Nhlf")
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

for(tissue in tissues){
   print(tissue)
   meta = read.table(paste0("data/rawdata/DNAbinding/",
                            tissue, "/metadata.tsv"),
                     header = T, sep = "\t", as.is = T)
   meta = meta[meta$Output.type %in% c("fold change over control",
                                       "signal") & 
                  meta$Assembly == "hg19" & 
                  meta$File.Status == "released" , ]
   meta = meta[grep(x = meta$Audit.ERROR,
                    pattern = "extreme",invert = T),]
   meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,
                    pattern = "severe",invert = T),]
   meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,
                    pattern = "insufficient",invert = T),]
   rownames(meta) = meta$File.accession
   for(size in sizes){
      print(size)
      
      # DNAbinding
      print("DNAbinding")
      files = list.files(paste0("data/procData/bins/",size, 
                                "/", tissue, "/DNAbinding/"), pattern="^ENC")
      files = substring(files, 1,11)
      filesByTarget = split(files, f = meta[files, "Experiment.target"])
      
      DNAbinding = sapply(names(filesByTarget) , function(target){
         ids=filesByTarget[[target]]
         res = sapply(ids, function(id){
            t = read.table(paste0("data/procData/", tissue, 
                                  "/DNAbinding/bins/", size, "_", 
                                  id, ".out.bed"))
            t = t[order(t$V1, t$V2),]
            return(t[,5])
         })
         return(rowMeans(res, na.rm = T))
      }) 
      colnames(DNAbinding) = do.call(rbind,strsplit(colnames(DNAbinding), 
                                                    split = "-"))[,1]
      save(DNAbinding,
           file = paste0("data/procData/",tissue,
                         "/DNAbinding/bins/", size, "_DNAbinding.RData"))
      
      # UCSC histone marks
      print("histone marks")
      histones = sapply(marks, function(m){
         res = sapply(lines, function(l){
            t = read.table(paste0("data/procData/", tissue, 
                                  "/UCSC_tracks/bins/", size, 
                                  "_wgEncodeBroadHistone", l,m, 
                                  "StdSig.out.bed"))
            t = t[order(t$V1, t$V2),]
            return(t[,5])
         })
         res = rowMeans(res)
         return(res)
      })
      colnames(histones) = paste("UCSC", colnames(histones), sep = "_")
      save(histones, 
           file = paste0("data/procData/",tissue,
                         "/UCSC_tracks/bins/", size, 
                         "_UCSC_histones.RData"))
      # TFBS
      print("TFBS")
      tfbs = read.table(paste0("data/procData/", tissue, 
                               "/UCSC_tracks/bins/", size, 
                               "_wgEncodeRegTfbsClusteredV3.out"), 
                        as.is = T, na.strings = ".")
      TF_BS = cbind(ETS_BS = as.integer(sapply(strsplit(tfbs[,5],split = ",", fixed = T), 
                                               function(x){sum(x %in% ETS)})),
                    TF_BS = as.integer(sapply(strsplit(tfbs[,5],split = ",", fixed = T), 
                                              function(x){sum(!is.na(x))}))) 
      save(TF_BS, file = paste0("data/procData/",tissue,
                                "/UCSC_tracks/bins/", size, 
                                "_Tfbs.RData"))
   }
}


