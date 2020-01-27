# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)

load("data/rdata/gtf_exons.RData")

library(parallel)


load(paste0("data/rdata/", tissue, "/Muts.RData"))
# new stuff
cl = makeCluster(8, type = "SOCK")
inexon = parApply(cl,Muts,1,function(m, gtf_exons){
   return(sum(gtf_exons$chr == m[2] & 
                 gtf_exons$start <= as.numeric(m[3]) &
                 gtf_exons$end >= as.numeric(m[3])) > 1)
}, gtf_exons = gtf_exons)
stopCluster(cl)
save(inexon,file = paste0("data/rdata/", tissue, "/inexon.RData"))



