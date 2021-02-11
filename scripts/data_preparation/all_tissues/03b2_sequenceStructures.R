# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


# abnormal sequence structures  #####
# first, execute script 03b_sequenceStructures.sh
aPhasedRepeats = read.table(paste0("data/procData/", tissue, "/structures/aPhasedRepeats.tsv"), 
                            sep = "\t", stringsAsFactors = F)
directRepeats = read.table(paste0("data/procData/", tissue, "/structures/directRepeats.tsv"), 
                           sep = "\t", stringsAsFactors = F)
gQuadruplex = read.table(paste0("data/procData/", tissue, "/structures/gQuadruplex.tsv"), 
                         sep = "\t", stringsAsFactors = F)
mirrorRepeats = read.table(paste0("data/procData/", tissue, "/structures/mirrorRepeats.tsv"), 
                           sep = "\t", stringsAsFactors = F)
shortTandemRepeats = read.table(paste0("data/procData/", tissue, "/structures/shortTandemRepeats.tsv"), 
                                sep = "\t", stringsAsFactors = F)
zDNAmotifs = read.table(paste0("data/procData/", tissue, "/structures/zDNAmotifs.tsv"), 
                        sep = "\t", stringsAsFactors = F)
structures = data.frame("aPhasedRepeats" = as.integer(aPhasedRepeats[,4] != 0),
                        "directRepeats" = as.integer(directRepeats[,4] != 0),
                        "gQuadruplex" = as.integer(gQuadruplex[,4] != 0),
                        "mirrorRepeats" = as.integer(mirrorRepeats[,4] != 0),
                        "shortTandemRepeats" = as.integer(shortTandemRepeats[,4] != 0),
                        "zDNAmotifs" = as.integer(zDNAmotifs[,4] != 0))
save(structures, file = paste0("data/rdata/", tissue, "/structures.RData"))

# 100 bp window
aPhasedRepeats = read.table(paste0("data/procData/", tissue, "/structures/aPhasedRepeats_100bp.tsv"), 
                            sep = "\t", stringsAsFactors = F)
directRepeats = read.table(paste0("data/procData/", tissue, "/structures/directRepeats_100bp.tsv"), 
                           sep = "\t", stringsAsFactors = F)
gQuadruplex = read.table(paste0("data/procData/", tissue, "/structures/gQuadruplex_100bp.tsv"), 
                         sep = "\t", stringsAsFactors = F)
mirrorRepeats = read.table(paste0("data/procData/", tissue, "/structures/mirrorRepeats_100bp.tsv"), 
                           sep = "\t", stringsAsFactors = F)
shortTandemRepeats = read.table(paste0("data/procData/", tissue, "/structures/shortTandemRepeats_100bp.tsv"), 
                                sep = "\t", stringsAsFactors = F)
zDNAmotifs = read.table(paste0("data/procData/", tissue, "/structures/zDNAmotifs_100bp.tsv"), 
                        sep = "\t", stringsAsFactors = F)
structures_100bp = data.frame("aPhasedRepeats" = as.integer(aPhasedRepeats[,4] != 0),
                              "directRepeats" = as.integer(directRepeats[,4] != 0),
                              "gQuadruplex" = as.integer(gQuadruplex[,4] != 0),
                              "mirrorRepeats" = as.integer(mirrorRepeats[,4] != 0),
                              "shortTandemRepeats" = as.integer(shortTandemRepeats[,4] != 0),
                              "zDNAmotifs" = as.integer(zDNAmotifs[,4] != 0))
colnames(structures_100bp) = paste(colnames(structures_100bp),
                                   "100bp", sep = "_")
save(structures_100bp, file = paste0("data/rdata/", tissue, "/structures_100bp.RData"))