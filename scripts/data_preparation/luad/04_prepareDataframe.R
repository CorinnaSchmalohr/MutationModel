# combine all data #####
load("data/rdata/luad/Muts.RData")
load("data/rdata/luad/mers.RData")
load("data/rdata/luad/GCcontent.RData")
load("data/rdata/luad/methylation.RData")
load("data/rdata/luad/conservation.RData")
load("data/rdata/luad/GTEx_eqtl.RData")
load("data/rdata/luad/DNAbinding.RData")
load("data/rdata/luad/UCSC_histones.RData")
load("data/rdata/luad/Nucleosome.RData")
load("data/rdata/luad/Tfbs.RData")
load("data/rdata/luad/structures.RData")
load("data/rdata/luad/mappability_repeats.RData")
load("data/rdata/luad/healthyExpr.RData")
load("data/rdata/luad/cancerExpr.RData")
load("data/rdata/luad/HiC.RData")
# load("data/rdata/luad/replTiming.RData")
# "data/rdata/luad/replDirection.RData"
load("data/rdata/luad/replDirection_Koren.RData")
load("data/rdata/luad/DNAaccessibility.RData")
data = data.frame(Muts,
                  trimer, pentamer, septamer,
                  precedingBase, followingBase, precedingDimer, followingDimer,
                  GCcontent,
                  methylation, conservation,GTEx_eqtl,
                  DNAbinding, Nhlf_histoneMarks, Nucleosome, 
                  ETS_BS, TF_BS,
                  structures, mappability, repeatMasker, Trf,
                  healthyExpr, cancerExpr,
                  HiC, replDirection_Koren, 
                  DNAaccessibility_lung, DNAaccessibility_UCSC,
                  stringsAsFactors = F)
data$ref = factor(data$ref)

# add label of how often we see a mutation at each position and remove duplicates
pos = paste(data$chr , data$pos, sep = "_")
noccur = table(pos)
multip = noccur[noccur > 1]
torm = NULL
for (i in names(multip)) {
   index = which(pos == i)
   if (length(unique(data$mutated[index])) > 1) {
      print(i)
   }
   if (all(data$mutated[index] == 0)) {
      print(i)
   }
   data$mutated[index[1]] <- sum(data$mutated[index])
   torm = c(torm,index[-1])
}
data = data[-torm,]

# remove positions with missing data
miss = apply(data[,-(which(colnames(data) == "alt"))], 1, function(x){sum(is.na(x))})
data = data[miss == 0,]
toRM = c("sample", "chr", "pos", "alt", "geneID", "baseGeneID")
removed = data[,toRM]
data = data[,(!colnames(data) %in% toRM)]
save(data, removed, file = "data/rdata/luad/completeData.RData")
#####
