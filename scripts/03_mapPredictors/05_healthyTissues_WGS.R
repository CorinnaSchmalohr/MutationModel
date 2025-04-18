args = as.numeric(commandArgs(trailingOnly=T))
tissues =c("adipocytes", "adrenal_gland", "bladder", 
           "bonemarrow", "brain", "breast", "colon", 
           "esophagus", "fibroblast", "heart",  "kidney",
           "liver", "lung", "pancreas", "placenta", "prostate", "rectum", 
           "skeletal_muscle", "skin", "small_intestine", "spleen", 
           "stomach", "testicle", "thyroid", "tonsil", "ureter") #"iPSC","blood", "embryonic_stem_cell", "endometrium"(has only indels)
tissue = tissues[args]
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
source("lib/dataMapping.R")
library(readxl)
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)
print(tissue)


# for each tissue, prepare corresponding data
print("mapping data")
tab = read_xlsx("data/rawdata/dataMappingAlltissues_WGS.xlsx", 
                sheet="allTissues", col_names=T)
tab$NA. = NULL
tab[tab == "NA"] = NA
if(tissue %in% colnames(tab)){
  tab = tab[,c(colnames(tab)[1:9],tissue)]
  tab = tab[!is.na(tab[,tissue]),]
} else{
  tab = tab[,c(colnames(tab)[1:9],"allTissues")]
  tab = tab[!is.na(tab[,"allTissues"]),]
}
# for predictors where we want multiple ranges, expand table
tab = apply(tab,1,function(x){
  if(is.na(x["range"])){
    return(x)
  } else{
    rangeWindow = strsplit(x["range"],";")[[1]]
    rangeInds = which(names(ranges) == rangeWindow[2]):which(names(ranges) == rangeWindow[1])
    subRanges = ranges[rangeInds]
    t(sapply(names(subRanges), function(r,y){
      y["range"] = subRanges[r]
      if(length(subRanges)>1){
        y["Name"] = paste0(y["Name"]," ",r)
        y["abbreviation"] = paste0(y["abbreviation"],"_",r)
      }
      return(y)
    },y=x))
  }
})
tab = do.call(rbind, tab)
tab = as.data.frame(tab)
pred = mapPredictors(x=tab, 
                     posFile=paste0("data/MutTables/healthyTissuesWGS/", 
                                    tissue, "_SomaMutDB.bed"))
load(paste0("data/MutTables/healthyTissuesWGS/", tissue, "_SomaMutDB.RData"))
data = list(meta = tab, pred = pred, muts = Muts)
save(data, file = paste0("data/MutTables/healthyTissuesWGS/", tissue, "_Muts_mapped.RData"))

dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
datchroms = data$muts$chr
dat = as.data.frame(dat)
save(dat,datchroms, file=paste0("data/MutTables/healthyTissuesWGS/", 
                                tissue, "_Muts_mapped_processed.RData"))
cat("\n")
