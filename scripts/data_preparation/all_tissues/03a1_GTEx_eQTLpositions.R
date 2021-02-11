args = as.numeric(commandArgs(trailingOnly=T)) # 1-6
tissues = c("luad", "breast", "skin", "colon", "ovary", "prostate")
# prostate not possible!
tissue = tissues[args]
print(tissue)
gtexTissue = list("luad"   = "Lung",
                  "breast" = "Breast_Mammary_Tissue",
                  "skin"   = c("Skin_Sun_Exposed_Lower_leg",
                               "Skin_Not_Sun_Exposed_Suprapubic"),
                  "colon"  = c("Colon_Sigmoid",
                               "Colon_Transverse"),
                  "ovary"  = "Ovary",
                  #"kidney" = NA,
                  "prostate" = "Prostate")

eqtl = lapply(gtexTissue[[tissue]], function(i){
   eqtl_tab = read.table(paste0("data/rawdata/GTEx_Analysis_v7_eQTL/",
                                i,
                                ".v7.signif_variant_gene_pairs.txt.gz"),
                         sep = "\t", as.is = T, header = T)
   pos = t(sapply(eqtl_tab[,1], function(x){
      t = strsplit(x,split = "_")[[1]]
      cr = paste0("chr", t[1])
      s = as.numeric(t[2]) - 11 # bed is 0-based!!
      e = s + nchar(t[3]) + 10
      r = substr(t[3],start = 1, stop = 1)
      return(c(cr,s,e,r))
   }))
   cbind(pos,eqtl_tab$pval_nominal)
})
eqtl  = do.call(rbind,eqtl)
eqtl = data.frame(eqtl, stringsAsFactors = F, row.names = NULL)
colnames(eqtl) = c("chr", "start", "end", "ref", "pval")
eqtl$start = as.integer(eqtl$start)
eqtl$end = as.integer(eqtl$end)
eqtl$pval = as.numeric(eqtl$pval)
eqtl = eqtl[order(eqtl[,1],eqtl[,2]),]
eqtl = eqtl[!duplicated(eqtl[,1:3]),]
write.table(eqtl, file = paste0("data/procData/",tissue,"/GTEx_eqtls.bed"),
            sep = "\t", quote = F, row.names = F, col.names = F)

