tissues = c("luad", "breast", "skin", "colon",
            "ovary", "kidney", "prostate")
binsizes = c("bins100bp", "bins1kb",  #"bins10bp",
         "bins10kb", "bins100kb", "bins1Mb")
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
library(BiocGenerics,quietly=T)
library(IRanges,quietly=T)

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
   for(bin in binsizes){
      print(bin)
      # prepare bin 2 range assignment
      load(paste0("data/procData/bins/", bin, "/", bin, ".RData"))
      bin2Cov = read.table(paste0("data/procData/bins/", bin, "/",
                                  bin, "_", tissue, "_covered.bed"))
      colnames(bin2Cov) = c("chr", "start", "end", "id")
      bin2Cov$start = bin2Cov$start+1
      binRanges = IRanges(bins$start, bins$end)
      covariateRanges = IRanges(bin2Cov$start, bin2Cov$end)
      overlapIndex = findOverlaps(binRanges, covariateRanges)
      overlapRanges = overlapsRanges(binRanges, covariateRanges, overlapIndex)

      # DNAbinding
      print("DNAbinding")
      files = list.files(paste0("data/procData/bins/",bin, 
                                "/", tissue, "/DNAbinding/"), pattern="^ENC")
      files = substring(files, 1,11)
      filesByTarget = split(files, f = meta[files, "Experiment.target"])
      
      DNAbinding = sapply(names(filesByTarget) , function(target){
         ids=filesByTarget[[target]]
         res = sapply(ids, function(id){
            t = read.table(paste0("data/procData/bins/",bin, "/", tissue,
                                  "/DNAbinding/", id,".out.bed"))
            return(t[,5])
         })
         return(rowMeans(res, na.rm = T))
      }) 
      colnames(DNAbinding) = do.call(rbind,strsplit(colnames(DNAbinding), 
                                                    split = "-"))[,1]
      save(DNAbinding,
           file = paste0("data/procData/bins/", bin,"/", tissue,
                         "/DNAbinding/DNAbinding.RData"))
      
      DNAbinding_covered = sapply(names(filesByTarget) , function(target){
         ids=filesByTarget[[target]]
         res = sapply(ids, function(id){
            t = read.table(paste0("data/procData/bins/",bin, "/", tissue,
                                  "/DNAbinding/covered_", id,".out.bed"))
            return(t[,5])
         })
         return(rowMeans(res, na.rm = T))
      }) 
      colnames(DNAbinding_covered) = do.call(rbind,strsplit(colnames(DNAbinding_covered), 
                                                            split = "-"))[,1]
      vals = DNAbinding_covered[subjectHits(overlapIndex),] * width(overlapRanges)/sum(width(overlapRanges))
      vals = apply(vals,2,function(x){
         sapply(split(x, queryHits(overlapIndex)), mean, na.rm = T)
      })
      res = matrix(NA, ncol=ncol(vals), nrow=nrow(bins),
                   dimnames=list(NULL, colnames(vals)))
      res[as.numeric(rownames(vals)),] = vals
      DNAbinding_covered = res
      save(DNAbinding_covered,
           file = paste0("data/procData/bins/", bin,"/", tissue,
                         "/DNAbinding/covered_DNAbinding.RData"))
      
      # histone marks
      print("histone marks")
      histones = sapply(marks, function(m){
         res = sapply(lines, function(l){
            t = read.table(paste0("data/procData/bins/", bin,"/", tissue, 
                                  "/UCSC_tracks/wgEncodeBroadHistone",
                                  l,m,"StdSig.out.bed"))
            return(t[,5])
         })
         res = rowMeans(res)
         return(res)
      })
      colnames(histones) = paste("UCSC", colnames(histones), sep = "_")
      save(histones, 
           file = paste0("data/procData/bins/", bin, "/", tissue,
                         "/UCSC_tracks/UCSC_histones.RData"))
      histones_covered = sapply(marks, function(m){
         res = sapply(lines, function(l){
            t = read.table(paste0("data/procData/bins/", bin,"/", tissue, 
                                  "/UCSC_tracks/covered_wgEncodeBroadHistone",
                                  l,m,"StdSig.out.bed"))
            return(t[,5])
         })
         res = rowMeans(res)
         return(res)
      })
      colnames(histones_covered) = paste("UCSC", colnames(histones_covered), sep = "_")
      vals = histones_covered[subjectHits(overlapIndex),] * width(overlapRanges)/sum(width(overlapRanges))
      vals = apply(vals,2,function(x){
         sapply(split(x, queryHits(overlapIndex)), mean, na.rm = T)
      })
      res = matrix(NA, ncol=ncol(vals), nrow=nrow(bins),
                   dimnames=list(NULL, colnames(vals)))
      res[as.numeric(rownames(vals)),] = vals
      histones_covered =res
      save(histones_covered, 
           file = paste0("data/procData/bins/", bin, "/", tissue,
                         "/UCSC_tracks/covered_UCSC_histones.RData"))
      
      print("TFBS")
      tfbs = read.table(paste0("data/procData/bins/", bin, "/", tissue, 
                               "/UCSC_tracks/wgEncodeRegTfbsClusteredV3.out"),
                        as.is = T, na.strings = ".")
      TF_BS = cbind(ETS_BS = as.integer(sapply(strsplit(tfbs[,5],split = ",", fixed = T), 
                                               function(x){sum(x %in% ETS)})),
                    TF_BS = as.integer(sapply(strsplit(tfbs[,5],split = ",", fixed = T), 
                                              function(x){sum(!is.na(x))}))) 
      save(TF_BS, file = paste0("data/procData/bins/", bin,"/",
                                tissue, "/UCSC_tracks/Tfbs.RData"))
      
      tfbs_covered = read.table(paste0("data/procData/bins/", bin, "/", tissue, 
                                       "/UCSC_tracks/covered_wgEncodeRegTfbsClusteredV3.out"),
                                as.is = T, na.strings = ".")
      TF_BS_covered = cbind(ETS_BS_covered = as.integer(sapply(strsplit(tfbs_covered[,5],split = ",", fixed = T), 
                                                               function(x){sum(x %in% ETS)})), 
                            TF_BS_covered =as.integer(sapply(strsplit(tfbs_covered[,5],split = ",", fixed = T), 
                                                             function(x){sum(!is.na(x))})))
      vals = TF_BS_covered[subjectHits(overlapIndex),] * width(overlapRanges)/sum(width(overlapRanges))
      vals = apply(vals,2,function(x){
         sapply(split(x, queryHits(overlapIndex)), mean, na.rm = T)
      })
      res = matrix(NA, ncol=ncol(vals), nrow=nrow(bins),
                   dimnames=list(NULL, colnames(vals)))
      res[as.numeric(rownames(vals)),] = vals
      TF_BS_covered = res
      save(TF_BS_covered,
           file = paste0("data/procData/bins/", bin,"/",
                         tissue, "/UCSC_tracks/covered_Tfbs.RData"))
   }
}


