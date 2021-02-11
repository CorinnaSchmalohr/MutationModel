# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


# combine all data #####
filenames = c(load(paste0("data/rdata/", tissue, "/Muts.RData")),
              load(paste0("data/rdata/", tissue, "/mers.RData")),
              if(tissue != "kidney"){
                 load(paste0("data/rdata/", tissue, "/GTEx_eqtl.RData"))
              },
              # load(paste0("data/rdata/", tissue, "/structures.RData")),
              load(paste0("data/rdata/", tissue, "/structures_100bp.RData")),
              load(paste0("data/rdata/", tissue, "/healthyExpr.RData")),
              load(paste0("data/rdata/", tissue, "/cancerExpr.RData")),
              load(paste0("data/rdata/", tissue, "/DNAaccessibility.RData"))[3],
              # load(paste0("data/rdata/", tissue, "/methylation.RData")),
              if(tissue != "ovary"){
                 load(paste0("data/rdata/", tissue, "/methbank.RData"))
              },
              load(paste0("data/rdata/", tissue, "/GCcontent.RData")),
              # load(paste0("data/rdata/", tissue, "/DNAbinding.RData")),
              load(paste0("data/rdata/", tissue, "/DNAbinding_100bp.RData")),
              # load(paste0("data/rdata/", tissue, "/UCSC_histones.RData")),
              # load(paste0("data/rdata/", tissue, "/UCSC_histones_100bp.RData")),
              # load(paste0("data/rdata/", tissue, "/Tfbs.RData")),
              load(paste0("data/rdata/", tissue, "/Tfbs_100bp.RData")),
              # load(paste0("data/rdata/", tissue, "/conservation.RData")),
              load(paste0("data/rdata/", tissue, "/conservation_100bp.RData")),
              load(paste0("data/rdata/", tissue, "/replication.RData")),
              load(paste0("data/rdata/", tissue, "/mappability_repeats.RData")),
              load(paste0("data/rdata/", tissue, "/inexon.RData")),
              # load(paste0("data/rdata/", tissue, "/Nucleosome.RData")),
              load(paste0("data/rdata/", tissue, "/Nucleosome_100bp.RData")),
              if(tissue %in% c("luad", "ovary")){
                 load(paste0("data/rdata/", tissue, "/HiC.RData"))
              } else if (tissue == "prostate"){
                 load(paste0("data/rdata/", tissue, "/HiC_interactions.RData"))
              } else if (tissue %in% c("breast", "kidney", "colon", "skin")){
                 load(paste0("data/rdata/", tissue, "/HiC.RData"))
              },
              load(paste0("data/rdata/", tissue,
                          "/binwisePreds_allChrs_bins10kb.RData")))
if(tissue != "ovary"){meth = meth[,"methbank_100bp"]}
GCcontent = GCcontent[,"GCcontent_100bp"]
data = data.frame(lapply(filenames, function(f){
   x = get(f)
   if(is.null(dim(x))){
      x = data.frame(x)
      colnames(x) = f
   }
   return(x)
})) #, check.names=F
toSort = names(data)[14:length(data)]
newOrder = c(names(data)[1:13], toSort[order(toSort)])
data = data[newOrder]


# add label of how often we see a mutation at each position and remove duplicates
# afaik there should be no duplicates, but leave this step for security
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
data$replDirection = as.factor(data$replDirection)
data$inexon = as.integer(data$inexon)
data$ref = factor(data$ref)
data$mutated = as.factor(data$mutated)
data$strand  = as.integer(data$strand == "+")
miss = apply(data[,-(which(colnames(data) == "alt"))],
             1, function(x){sum(is.na(x))})
data = data[miss == 0,]
toRM = c("chr", "pos", "alt", "geneID", "baseGeneID")
removed = data[,toRM]
data = data[,(!colnames(data) %in% toRM)]
save(data, removed, file = paste0("data/rdata/", tissue, 
                                  "/completeData_unscaled_withBinwise.RData"))

toScale = sapply(data,is.numeric) & !sapply(data,is.integer)
data[,toScale] = scale(data[,toScale])
save(data, removed, file = paste0("data/rdata/", tissue, "/completeData_withBinwise.RData"))
#####
