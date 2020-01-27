# make sure to use files that combined replicates 1 and 2 over each individual file.
# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)

conservation = read.table(paste0("data/procData/", tissue, 
                          "/conservation/hg19.100way.phyloP100way.out.bed"), as.is = T)
conservation = conservation[,5]
save(conservation, file = paste0("data/rdata/", tissue,
                                 "/conservation.RData"))
