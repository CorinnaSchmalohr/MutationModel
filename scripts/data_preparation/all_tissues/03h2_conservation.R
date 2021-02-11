# make sure to use files that combined replicates 1 and 2 over each individual file.
# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)

conservation = read.table(paste0("data/procData/", tissue, 
                          "/conservation/hg19.100way.phyloP100way.out.bed"), as.is = T)
conservation = conservation[order(conservation$V1, conservation$V2),]
conservation = conservation[,5]
save(conservation, file = paste0("data/rdata/", tissue,
                                 "/conservation.RData"))


conservation_100bp = read.table(paste0("data/procData/", tissue, 
                                 "/conservation/hg19.100way.phyloP100way_100bp.out.bed"), as.is = T)
conservation_100bp = conservation_100bp[order(conservation_100bp$V1, conservation_100bp$V2),]
conservation_100bp = conservation_100bp[,5]
save(conservation_100bp, file = paste0("data/rdata/", tissue,
                                       "/conservation_100bp.RData"))
