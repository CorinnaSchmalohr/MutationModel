eqtl_tab = read.table("data/rawdata/GTEx_Analysis_v7_eQTL/Lung.v7.signif_variant_gene_pairs.txt.gz",
                  sep = "\t", as.is = T, header = T)
pos = t(sapply(eqtl_tab[,1], function(x){
   t = strsplit(x,split = "_")[[1]]
   c = paste0("chr", t[1])
   s = as.numeric(t[2]) - 11 # bed is 0-based!!
   e = s + nchar(t[3]) + 10
   r = substr(t[3],start = 1, stop = 1)
   return(c(c,s,e,r))
}))
eqtl = data.frame(cbind(pos,eqtl_tab$pval_nominal), stringsAsFactors = F, row.names = NULL)
colnames(eqtl) = c("chr", "start", "end", "ref", "pval")
eqtl$start = as.integer(eqtl$start)
eqtl$end = as.integer(eqtl$end)
eqtl$pval = as.numeric(eqtl$pval)
eqtl = eqtl[order(eqtl[,1],eqtl[,2]),]
write.table(eqtl, file = "data/procData/lung/GTEx_lung_eqtls.bed", sep = "\t", quote = F, row.names = F, col.names = F)

