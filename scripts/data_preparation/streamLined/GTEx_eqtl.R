gtexFiles = list.files("data/rawdata/GTEx_Analysis_v7_eQTL/", 
                       pattern=".v7.signif_variant_gene_pairs.txt.gz")
FilesByTissue = split(gtexFiles,
                      sapply(strsplit(gtexFiles,split="[_.]"), function(x){x[1]}))
# set chromosome lengths
chrs = sort(paste0("chr", c(1:22)))
chrLengths = read.table("data/rawdata/GRCh37.p11.genome.fa.fai")
chrLengths = setNames(chrLengths[,2],chrLengths[,1])
chrLengths = chrLengths[chrs]

# tissue = "Adipose"
for(tissue in names(FilesByTissue)){
   print(tissue)
   files = FilesByTissue[[tissue]]
   # load GTEx eqtl results
   eqtlDat = do.call(rbind,lapply(files, function(file){
      temp = read.table(paste0("data/rawdata/GTEx_Analysis_v7_eQTL/",file),
                 sep = "\t", as.is = T, header = T)
      return(temp)
   }))
   
   # reformat to get SNP positions
   temp = t(sapply(eqtlDat[,1], function(x){
      t = strsplit(x,split = "_")[[1]]
      c(paste0("chr", t[1]), t[2])
   }))
   eQTL = data.frame(chr = temp[,1], 
                     start = as.integer(temp[,2]),
                     pval = eqtlDat$pval_nominal)
   
   eQTLapprox = sapply(chrs, function(cr){
      # subset to each chromosome
      subEqtl = eQTL[eQTL$chr == cr,]
      # sort and remove duplicate SNPs (keeping the smallest pval)
      subEqtl = subEqtl[order(subEqtl[,"start"]),]
      duplPos = unique(subEqtl[duplicated(subEqtl[,"start"]),"start"])
      toRemove = unlist(sapply(duplPos, function(pos){
         inds = which(subEqtl[,"start"] == pos)
         inds[-which.min(subEqtl[inds,"pval"])]
      }))
      subEqtl = subEqtl[-toRemove,]
      # reformat to cover all positions (non-eQTL --> p-value = 1) 
      i = 1; res = rbind(do.call(rbind,apply(subEqtl, 1,function(x){
         s = as.integer(x["start"])
         p = as.numeric(x["pval"])
         if(i<s){
            out = rbind(c(i,1),
                        c(s,p))
         } else{ out = c(s,p)}
         i <<- s+1
         return(out)
      })),c(subEqtl[nrow(subEqtl),"start"]+1,1))
      # make approxfunction
      logPvals = -log(res[,2])
      aFun = stepfun(x = res[-1,1], y = logPvals, f=0)
      # plot(aFun)
      # points(x=res[,1], y = -log(res[,2]), col = "red")
      # plot(aFun, xlim = c(722600,722610))
      # points(x=res[,1], y = -log(res[,2]), col = "red")
      return(aFun)
   }, simplify=F)
   save(eQTLapprox, file = paste0("data/streamLined/eQTLpositions_",
                                  tissue, "_stepfun.RData"))
}



