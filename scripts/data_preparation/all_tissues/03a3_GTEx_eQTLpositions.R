args = as.numeric(commandArgs(trailingOnly=T)) # 1-7
tissues = c("luad", "breast", "skin", "colon", "ovary", "prostate", "kidney")
# kidney not possible!
tissue = tissues[args]

# whether variant lies within 10 bp of an eQTL in this tissue. column 4 is minimum p-value of all hits.
if(tissue == "kidney"){
   GTEx_eqtl = NA
} else{
   GTEx_eqtl = read.table(paste0("data/procData/",tissue,"/GTEx_eqtls.out"),
                          sep = "\t", na.strings = ".", as.is = T)
   GTEx_eqtl = GTEx_eqtl[,5]
   GTEx_eqtl = as.integer(!is.na(GTEx_eqtl))
   save(GTEx_eqtl, file = paste0("data/rdata/",tissue,"/GTEx_eqtl.RData"))
}


     