# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)



# combine all data #####
load(paste0("data/rdata/", tissue, "/Muts.RData"))
load(paste0("data/rdata/", tissue, "/mers.RData"))
if(tissue != "kidney"){
   load(paste0("data/rdata/", tissue, "/GTEx_eqtl.RData"))
   GTEx_eqtl = as.integer(GTEx_eqtl)
} else
   GTEx_eqtl = rep(NA, nrow(Muts))
load(paste0("data/rdata/", tissue, "/structures.RData"))
load(paste0("data/rdata/", tissue, "/healthyExpr.RData"))
load(paste0("data/rdata/", tissue, "/cancerExpr.RData"))
load(paste0("data/rdata/", tissue, "/DNAaccessibility.RData"))
load(paste0("data/rdata/", tissue, "/methylation.RData"))
if(tissue != "ovary"){
   load(paste0("data/rdata/", tissue, "/methbank.RData"))
} else
   meth = rep(NA, nrow(Muts))
load(paste0("data/rdata/", tissue, "/GCcontent.RData"))
load(paste0("data/rdata/", tissue, "/DNAbinding.RData"))
load(paste0("data/rdata/", tissue, "/UCSC_histones.RData"))
load(paste0("data/rdata/", tissue, "/Tfbs.RData"))
load(paste0("data/rdata/", tissue, "/conservation.RData"))
load(paste0("data/rdata/", tissue, "/replication.RData"))
replication$replDirection = as.factor(replication$replDirection)
load(paste0("data/rdata/", tissue, "/mappability_repeats.RData"))
load(paste0("data/rdata/", tissue, "/inexon.RData"))
inexon = as.integer(inexon)
load(paste0("data/rdata/", tissue, "/Nucleosome.RData"))
if(tissue %in% c("luad", "ovary")){
   load(paste0("data/rdata/", tissue, "/HiC.RData"))
} else if (tissue == "prostate"){
   load(paste0("data/rdata/", tissue, "/HiC_interactions.RData"))
   HiC = cbind(HiC_ints = HiCints)
} else if (tissue %in% c("breast", "kidney", "colon", "skin")){
   load(paste0("data/rdata/", tissue, "/HiCints.RData"))
   load(paste0("data/rdata/", tissue, "/HiC_TADs.RData"))
   load(paste0("data/rdata/", tissue, "/HiCcomp.RData"))
   HiC = cbind(HiC_ints = HiCints, HiC_inTAD, 
               HiC_TADbound = HiC_onTADboundary, HiC_compPCA = HiCcomp)
}
data = data.frame(Muts,
                  trimer, pentamer, septamer,
                  precedingBase, followingBase, 
                  GCcontent, inexon,
                  methylation, meth,
                  conservation, GTEx_eqtl,
                  DNAbinding, histones, Nucleosome, 
                  ETS_BS, TF_BS,
                  structures, mappability, repeatMasker, Trf,
                  healthyExpr, cancerExpr,
                  HiC,
                  replication, 
                  DNAaccessibility_tissue, DNAaccessibility_tissue_1kb, DNAaccessibility_UCSC,
                  stringsAsFactors = F)
data$ref = factor(data$ref)
if(tissue == "kidney"){
   data$GTEx_eqtl = NULL
}
if(tissue == "ovary"){
   data$meth = NULL
}
if(tissue == "breast"){
   data$HiC_compPCA = NULL
}

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
if(length(torm)>0){
   data = data[-torm,]
}

# remove positions with missing data
data$methbank_local = NULL
data$mutated = as.factor(data$mutated)
data$strand  = as.integer(data$strand == "+")
miss = apply(data[,-(which(colnames(data) == "alt"))],
             1, function(x){sum(is.na(x))})
data = data[miss == 0,]
toRM = c("chr", "pos", "alt", "geneID", "baseGeneID")
removed = data[,toRM]
data = data[,(!colnames(data) %in% toRM)]
save(data, removed, file = paste0("data/rdata/", tissue, 
                                  "/completeData_unscaled.RData"))

toScale = sapply(data,is.numeric) & !sapply(data,is.integer)
data[,toScale] = scale(data[,toScale])
save(data, removed, file = paste0("data/rdata/", tissue, "/completeData.RData"))
#####
