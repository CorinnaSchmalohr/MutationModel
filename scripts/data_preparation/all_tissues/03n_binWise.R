tissues = c(  "luad","colon", "breast", "skin",
              "ovary", "kidney", "prostate")
sizes = c("bins10kb", "bins100kb", "bins1Mb")
chrs = paste0("chr", c(1:22))

temp = sapply(tissues, function(tissue){
   load(paste0("data/rdata/", tissue, "/Muts.RData"))
   
   binWiseValues = lapply(sizes, function(size){
      # create table of all features that were represented binWise
      load(paste0("data/rdata/", size, "_allChrs.RData"))
      filenames = c(load(paste0("data/procData/", tissue,
                                "/DNAaccessibility/bins/", size,
                                "DNAaccessibility.RData")),
                    load(paste0("data/procData/",tissue,
                                "/DNAbinding/bins/", size, 
                                "_DNAbinding.RData")),
                    load(paste0("data/procData/",tissue,
                                "/UCSC_tracks/bins/", size, 
                                "_UCSC_histones.RData")),
                    load(paste0("data/procData/",tissue,
                                "/UCSC_tracks/bins/", size, 
                                "_Tfbs.RData")),
                    if(tissue != "ovary"){
                       load(paste0("data/procData/", 
                                   tissue, "/methylation/bins/",
                                   size, "_methylation.RData"))},
                    load(paste0("data/procData/", tissue,
                                "/UCSC_tracks/bins/", size,
                                "_Nucleosome.RData")),
                    load(paste0("data/procData/",
                                size,"_GCcontent.RData")))
      binWise = data.frame(bins,do.call(cbind,sapply(filenames, function(x){
         get(x)
      })))
      res = do.call(rbind,lapply(unique(Muts$chr), function(chr){
         subBinWise = binWise[binWise$chr == chr,]
         subMuts = Muts[Muts$chr == chr,]
         fun = approxfun(c(subBinWise$start, subBinWise$end[nrow(subBinWise)]), 
                         seq(nrow(subBinWise)+1), 
                         method = "constant")
         inds = fun(subMuts$pos)
         subBinWise[inds, -(1:3)]
      }))
      colnames(res) = paste(colnames(res), size, sep = "_")
      n = paste("binWiseValue", size, sep="_")
      assign(n,res)
      save(list=n, file = paste0("data/rdata/", tissue,
                                 "/binwisePreds_allChrs_", size, ".RData"))
      return(res)
   })
   binWiseValues = data.frame(binWiseValues)
   # temp = binWiseValues[order(names(binWiseValues))]
   save(binWiseValues, file = paste0("data/rdata/", tissue,
                                     "/binwisePreds_allChrs.RData"))
})
