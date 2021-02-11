# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "ovary", "breast", "colon", "kidney","prostate",  "skin" )
tissue = tissues[args]
print(tissue)

load(paste0("data/rdata/", tissue, "/Muts.RData"))

if(tissue %in% c("luad", "ovary")){
   tissue2file = c("luad" = "LG", "ovary" = "OV")
   # PC1 used for compartment A/B prediction
   HiC_compPCA_table = read.table("data/rawdata/HiC_Schmitt2016/HiC_ABcompartments_hg19_Schmitt2016_PCA.csv",
                                  header = T,sep = ",",as.is = T)
   HiC_compPCA_table$chr = paste0("chr", HiC_compPCA_table$chr)
   # A/B labels of compartment prediction
   HiC_compLabels_table = read.table("data/rawdata/HiC_Schmitt2016/HiC_ABcompartments_hg19_Schmitt2016_ABlabels.csv",
                                     header = T,sep = ",",as.is = T)
   HiC_compLabels_table$chr = paste0("chr", HiC_compLabels_table$chr)
   # FIRE scores (frequently interacting regions)
   HiC_FIRE_table = read.table("data/rawdata/HiC_Schmitt2016/HiC_FIREscores_hg19_Schmitt2016_primaryCohort.csv",
                               header = T,sep = ",",as.is = T)
   HiC_FIRE_table$chr = paste0("chr", HiC_FIRE_table$chr)
   HiC_TADbound_table = read.table(paste0("data/rawdata/HiC_Schmitt2016/HiC_TADboundaries_hg19_Schmitt2016_", 
                                   tissue2file[tissue], ".csv"),
                                   header = F,sep = ",",as.is = T)
   colnames(HiC_TADbound_table) = c("chr", "start", "end")
   HiC = t(apply(Muts,1, function(x){
      c = x["chr"]
      p = as.numeric(x["pos"])
      HiC_compPCA = HiC_compPCA_table[HiC_compPCA_table$chr == c & 
                                         HiC_compPCA_table$start < p & 
                                         HiC_compPCA_table$end > p, tissue2file[tissue]]
      if (length(HiC_compPCA) == 0) {
         print(x)
         HiC_compPCA = NA
      }
      HiC_compLabels = HiC_compLabels_table[HiC_compLabels_table$chr == c & 
                                               HiC_compLabels_table$start < p & 
                                               HiC_compLabels_table$end > p, tissue2file[tissue]]
      if (length(HiC_compLabels) == 0) {HiC_compLabels = NA}
      HiC_FIRE =  HiC_FIRE_table[HiC_FIRE_table$chr == c & 
                                    HiC_FIRE_table$start < p & 
                                    HiC_FIRE_table$end > p, tissue2file[tissue]]
      if (length(HiC_FIRE) == 0) {HiC_FIRE = NA}
      HiC_TADbound = nrow(HiC_TADbound_table[HiC_TADbound_table$chr == c & 
                                                HiC_TADbound_table$start < p & 
                                                HiC_TADbound_table$end > p, ])
      if (length(c(HiC_compPCA, HiC_compLabels, HiC_FIRE, HiC_TADbound)) > 4)
         print(x)
      return(c(HiC_compPCA = HiC_compPCA, 
               HiC_compLabels = HiC_compLabels,
               HiC_FIRE = HiC_FIRE, 
               HiC_TADbound = HiC_TADbound))
   }))
   HiC = data.frame(HiC_compPCA = as.integer(HiC[,"HiC_compPCA"]),
                    HiC_compLabels = as.integer(HiC[,"HiC_compLabels"] == "A"), 
                    HiC_FIRE = as.integer(HiC[,"HiC_FIRE"]),
                    HiC_TADbound = as.integer(HiC[,"HiC_TADbound"]))
   save(HiC, file = paste0("data/rdata/", tissue, "/HiC.RData"))
   rm(HiC_compPCA_table, HiC_compLabels_table, HiC_FIRE_table, HiC_TADbound_table, HiC) 
}

if(tissue == "prostate"){
   library(parallel)
   ints = read.table("data/rawdata/HiC_prostate/GSE118629_22Rv1_HiC_40k.normalized.matrix.txt.gz")
   bins = read.table("data/rawdata/HiC_prostate/GSE118629_hg19_40k.bed.gz", as.is = T)
   colnames(bins) = c("chr", "start", "end", "id")
   id2chr = bins$chr
   names(id2chr) = bins$id
   col1chr = id2chr[ints$V1]
   col2chr = id2chr[ints$V2]
   intsPerChr = lapply(unique(Muts$chr), function(cr){
      cat(cr, ' ')
      subMuts = Muts[Muts$chr == cr,]
      toSubset = which(bins$chr == cr)
      subInts = ints[col1chr == cr | col2chr == cr,]
      subBins = bins[toSubset,]
      cl = makeCluster(4,"SOCK")
      res = parSapply(cl, subMuts$pos, function(m_pos, subBins, subInts){
         hit = subBins[subBins$start <= m_pos & subBins$end > m_pos,"id"]
         intVals = subInts[subInts$V1 == hit | subInts$V2 == hit,]
         intVals = intVals[intVals$V1 != intVals$V2,]
         return(mean(intVals$V3, na.rm = T))
      }, subBins = subBins, subInts = subInts)
      stopCluster(cl)
      return(res)
   })
   HiCints = unlist(intsPerChr)
   save(HiCints, file = paste0("data/rdata/", tissue, "/HiC_interactions.RData"))
}


if(tissue %in% c("breast", "skin", "colon", "kidney")){
   fun = function(m,subBins, subInts){
      binHits = which(subBins$start <= m & subBins$end > m)
      binVals = c(subInts[binHits,-binHits])
      return(mean(binVals, na.rm = T))
   }
   meta = read.table(paste0("data/rawdata/HiC_Encode/", 
                     tissue, "/metadata.tsv"), 
                     header = T, as.is = T, sep = "\t")
   meta = meta[meta$Output.type %in% c("chromatin interactions",
                                       "topologically associated domains",
                                       "genome compartments"),] 
   exp = unique(meta$Experiment.accession)
   type2ending = c("hdf5" = ".h5", 
                   "bed bed3+" = ".bed.gz", 
                   "bigWig" = ".bigWig",
                   "hic" = ".hic", 
                   "bedpe" = ".bedpe.gz")
   library(methods)
   library(h5)
   library(parallel)
   library(strawr)
   print("ints")
   HiCints = lapply(exp, function(e){
      print(e)
      sub = meta[meta$Output.type == "chromatin interactions" & 
                    meta$Experiment.accession == e,]
      ints = apply(sub, 1,function(row){
         f = row["File.accession"]
         print(f)
         type = row["File.format"]
         genome = row["Assembly"]
         if(type == "hdf5"){
            cat("hdf5",' ')
            dat = h5file(paste0("data/rawdata/HiC_Encode/", tissue,
                                "/", f,type2ending[type]))
            bins = as.data.frame(dat["bin_positions"][])
            colnames(bins) = c("chr", "start", "end")
            chrs = dat["chrs"][]
            bins$chr = chrs[bins$ch+1]
            intsPerChr = sapply(unique(Muts$chr), function(cr){
               toSelect = which(bins$chr == cr)
               subBins = bins[toSelect,]
               subInts = dat["interactions"][toSelect,toSelect]
               subMuts = Muts[Muts$chr == cr,]
               cl = makeCluster(8,"SOCK")
               res = parSapply(cl,subMuts$pos, fun, 
                               subBins = subBins, subInts = subInts)
               stopCluster(cl)
               return(res)
            })
            h5close(dat)
            intsPerChr = lapply(intsPerChr, scale)
            intsPerChr = unlist(intsPerChr)
            return(intsPerChr)
         } else if (type == "hic"){
            cat("hic",' ')
            ending = type2ending[type]
            res = apply(Muts,1, function(m){
               m_chr = substring(m["chr"],4,10)
               m_pos = as.integer(m["pos"])
               pos = paste(m_chr, m_pos, m_pos+10000, sep = ":")
               ints = straw(norm = "NONE", fname = paste0("data/rawdata/HiC_Encode/", tissue,
                                                          "/", f,ending), 
                            unit = "BP", binsize=10000, chr1loc=pos, chr2loc=m_chr)
               return(sum(ints[,3]))
            })
            return(res)
         }
      })
      return(rowMeans(ints, na.rm=T))
   })
   HiCints = do.call(cbind,HiCints)
   HiCints = scale(HiCints)
   HiCints = rowMeans(HiCints, na.rm=T)
   save(HiCints, file = paste0("data/rdata/", tissue, "/HiCints.RData"))
   rm(HiCints)
   
   print("Tads")
   HiCTADs = lapply(exp, function(e){
      print(e)
      sub = meta[meta$Output.type == "topologically associated domains" & 
                    meta$Experiment.accession == e,]
      if(nrow(sub) == 0){
         return(list(inTAD = rep(NA, nrow(Muts)), onTADbound = rep(NA, nrow(Muts))))
      }
      tads = apply(sub, 1,function(row){
         f = row["File.accession"]
         print(f)
         type = row["File.format"]
         genome = row["Assembly"]
         dat = read.table(paste0("data/rawdata/HiC_Encode/", 
                                 tissue, "/", f,type2ending[type]))
         colnames(dat) = c("chr", "start", "end", "description", "someValue")
         res = t(apply(Muts,1,function(m){
            cr = m["chr"]
            pos = as.numeric(m["pos"])
            intad = any(dat$chr == cr & dat$start <= pos & dat$end > pos)
            onbound = any(dat$chr == cr & 
                             (abs(dat$start - pos) <= 5000 | 
                                 abs(dat$end - pos) <= 5000))
            return(as.integer(c(intad, onbound)))
         }))
         return(as.data.frame(res))
      })
      inTAD = sapply(tads,function(x){return(x$V1)})
      inTAD = rowMeans(inTAD, na.rm = T)
      onTADbound = sapply(tads,function(x){return(x$V2)})
      onTADbound = rowMeans(onTADbound, na.rm = T)
      return(list(inTAD = inTAD, onTADbound = onTADbound))
   })
   HiC_inTAD = sapply(HiCTADs, function(x){x$inTAD})
   HiC_inTAD = round(rowMeans(HiC_inTAD, na.rm = T))
   HiC_onTADboundary = sapply(HiCTADs, function(x){x$onTADbound})
   HiC_onTADboundary = ceiling(rowMeans(HiC_onTADboundary, na.rm = T))
   save(HiC_inTAD, HiC_onTADboundary, 
        file = paste0("data/rdata/", tissue, "/HiC_TADs.RData"))
   rm(HiCTADs, HiC_inTAD, HiC_onTADboundary)
   
   library(rtracklayer)
   print("comp")
   HiCcomp = sapply(exp, function(e){
      print(e)
      sub = meta[meta$Output.type == "genome compartments" & 
                    meta$Experiment.accession == e,]
      if(nrow(sub) == 0){
         return(rep(NA, nrow(Muts)))
      }
      comp = apply(sub, 1,function(row){
         f = row["File.accession"]
         print(f)
         type = row["File.format"]
         genome = row["Assembly"]
         dat = as.data.frame(import.bw(paste0("data/rawdata/HiC_Encode/",
                                              tissue, "/", f, type2ending[type])))
         res = apply(Muts,1,function(m){
            cr = m["chr"]
            pos = as.numeric(m["pos"])
            s = dat[dat$seqnames == cr & 
                           dat$start <= pos & 
                           dat$end > pos,"score"]
            if(length(s) < 1){
               s = NA
            }
            return(s)
         })
         return(res)
      })
      comp = scale(comp)
      return(rowMeans(comp, na.rm = T))
   })
   HiCcomp = rowMeans(HiCcomp, na.rm = T)
   save(HiCcomp, file = paste0("data/rdata/", tissue, "/HiCcomp.RData"))
   
   load(paste0("data/rdata/", tissue, "/HiCints.RData"))
   load(paste0("data/rdata/", tissue, "/HiC_TADs.RData"))
   load(paste0("data/rdata/", tissue, "/HiCcomp.RData"))
   HiC = cbind(HiC_ints = HiCints, HiC_inTAD, 
               HiC_TADbound = HiC_onTADboundary, HiC_compPCA = HiCcomp)
   HiC = HiC[,apply(HiC,2,function(x){all(is.na(x))})]
   save(HiC, file = paste0("data/rdata/", tissue, "/HiC.RData"))
}
