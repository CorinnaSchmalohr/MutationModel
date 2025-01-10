source("lib/dataMapping.R")
library(readxl)
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "luad",
            "ovary", "prostate", "skin")
tissue = tissues[args]

print("mapping data")
print(tissue)
tab = read_xlsx("data/rawdata/dataMapping.xlsx",
                sheet=tissue, col_names = T)
tab = as.data.frame(tab)
tab[tab == "NA"] = NA

for(i in 1:4){
  pred = mapPredictors(x=tab, 
                       posFile=paste0("data/procData/Chr1arm/Chr1arm_allPos_part",i,".bed"))
  Muts = get(load(paste0("data/rdata/Chr1arm/Chr1arm_allPos_part",i,".RData")))
  data = list(meta = tab, pred = pred, muts = Muts)
  save(data, file = paste0("data/procData/Chr1arm/Chr1arm_allPos_part",i,"_", tissue, ".RData"))
  
  cat("\n")
  dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
  datchroms = data$muts$chr
  
  # remove datapoints that could be problematic
  #toExclude = dat$ConsensusExcludable == 1 | dat$repeatMasker ==1 | dat$tandemRepeatFinder == 1 | 
  #  dat$mappability_100mer < 1 |dat$mappabiliy_40mer < 1 | dat$mappability_24mer < 1
  #dat = dat[!toExclude,]
  #datchroms = datchroms[!toExclude]
  #dat$ConsensusExcludable = NULL
  #dat$repeatMasker = NULL
  #dat$tandemRepeatFinder = NULL
  #dat$mappability_100mer = NULL
  #dat$mappabiliy_40mer = NULL
  #dat$mappability_24mer = NULL
  
  # center and log and replace missing values with mean:
  dat = sapply(colnames(dat),FUN=function(j){ 
    x = dat[,j]
    if(!is.numeric(x)){return(x)}
    if(max(x, na.rm = T)-min(x, na.rm = T) == 0){
      x = x-min(x, na.rm = T)
    } else{
      x = (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)) # scale from 0 to 1
    }
    x = log(x +1) # log
    x[is.na(x)] = mean(x,na.rm = T) # replace missing
    # x = scale(x)
    return(x) 
  }, simplify=F)
  dat = as.data.frame(dat)
  
  # Edit bases to CG or not CG as we are not looking at the strand 
  dat[c("precedingBase_CG", "ref_CG", "followingBase_CG")] = 0
  dat[dat$context_ref == "C" | dat$context_ref == "G",]["ref_CG"] = 1
  dat$ref_CG = as.factor(dat$ref_CG)
  dat[dat$context_precedingBase == "C" | dat$context_precedingBase == "G",]["precedingBase_CG"] = 1
  dat$precedingBase_CG = as.factor(dat$precedingBase_CG)
  dat[dat$context_followingBase == "C" | dat$context_followingBase == "G",]["followingBase_CG"] = 1
  dat$followingBase_CG = as.factor(dat$followingBase_CG)
  
  dat = dat[, !(colnames(dat) %in% c("context_precedingBase","context_ref","context_followingBase"))]
  
  save(dat,datchroms, file=paste0("data/procData/Chr1arm/Chr1arm_allPos_part",i,"_",
                                  tissue,"_processed.RData"))
}


