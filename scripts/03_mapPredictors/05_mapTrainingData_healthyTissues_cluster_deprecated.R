args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus", "kidney", "liver", "luad","prostate", "skin")
tissue = tissues[args]
print(tissue)

source("lib/dataMapping.R")
library(readxl)


# Moore et al. #####
tissues_Moore = c("colon","esophagus", "kidney","liver", "luad", 
                  "prostate", "skin")
if(tissue %in% tissues_Moore){
  print("Moore et al.")
  #load table of predictor files to use
  tab = read_xlsx("data/rawdata/dataMapping.xlsx",
                  sheet=tissue, col_names = T)
  tab = as.data.frame(tab)
  tab[tab == "NA"] = NA
  
  # map predictors
  pred = mapPredictors(x=tab,
                       posFile=paste0("data/procData/", tissue,
                                      "/mutData_validation_Moore.bed"))
  
  # save raw mapped data
  load(paste0("data/rdata/", tissue, "/mutData_validation_Moore.RData"))
  data = list(meta = tab, pred = pred, muts = Muts)
  save(data, file = paste0("data/procData/validationData/validationData_Moore_", tissue, ".RData"))
  # load(paste0("data/procData/validationData/validationData_Moore_", tissue, ".RData"))
  cat("\n")
  
  # process data further
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
  
  # save
  save(dat,datchroms, file=paste0("data/procData/validationData/validationData_Moore_processed_", tissue, ".RData"))
}  
#####

# Li et al. #####
tissues_Li = c("colon","esophagus", "liver", "luad")
if(tissue %in% tissues_Li){
  print("Li et al")
  #load table of predictor files to use
  tab = read_xlsx("data/rawdata/dataMapping.xlsx",
                  sheet=tissue, col_names = T)
  tab = as.data.frame(tab)
  tab[tab == "NA"] = NA
  
  # map predictors
  pred = mapPredictors(x=tab,
                       posFile=paste0("data/procData/", tissue,
                                      "/mutData_validation_Li.bed"))
  
  # save raw mapped data
  load(paste0("data/rdata/", tissue, "/mutData_validation_Li.RData"))
  data = list(meta = tab, pred = pred, muts = Muts)
  save(data, file = paste0("data/procData/validationData/validationData_Li_", tissue, ".RData"))
  # load(paste0("data/procData/validationData/validationData_Li_", tissue, ".RData"))
  cat("\n")
  
  # process data further
  dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
  datchroms = data$muts$chr
  
  # remove datapoints that could be problematic
  toExclude = dat$ConsensusExcludable == 1 | dat$repeatMasker ==1 | dat$tandemRepeatFinder == 1 | 
    dat$mappability_100mer < 1 |dat$mappabiliy_40mer < 1 | dat$mappability_24mer < 1
  dat = dat[!toExclude,]
  datchroms = datchroms[!toExclude]
  dat$ConsensusExcludable = NULL
  dat$repeatMasker = NULL
  dat$tandemRepeatFinder = NULL
  dat$mappability_100mer = NULL
  dat$mappabiliy_40mer = NULL
  dat$mappability_24mer = NULL
  
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
  
  # save
  save(dat,datchroms, file=paste0("data/procData/validationData/validationData_Li_processed_", tissue, ".RData"))
  #####
}



# SomaMutDB #####
print("SomaMutDB")
#load table of predictor files to use
tab = read_xlsx("data/rawdata/dataMapping.xlsx",
                 sheet=tissue, col_names = T)
tab = as.data.frame(tab)
tab[tab == "NA"] = NA

# map predictors
pred = mapPredictors(x=tab,
                     posFile=paste0("data/procData/", tissue,
                                    "/mutData_validation_SomaMutDB.bed"))

# save raw mapped data
load(paste0("data/rdata/", tissue, "/mutData_validation_SomaMutDB.RData"))
data = list(meta = tab, pred = pred, muts = Muts)
save(data, file = paste0("data/procData/validationData/validationData_SomaMutDB_", tissue, ".RData"))

cat("\n")

# process data further
dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
datchroms = data$muts$chr

# remove datapoints that could be problematic
toExclude = dat$ConsensusExcludable == 1 | dat$repeatMasker ==1 | dat$tandemRepeatFinder == 1 | 
  dat$mappability_100mer < 1 |dat$mappabiliy_40mer < 1 | dat$mappability_24mer < 1
dat = dat[!toExclude,]
datchroms = datchroms[!toExclude]
dat$ConsensusExcludable = NULL
dat$repeatMasker = NULL
dat$tandemRepeatFinder = NULL
dat$mappability_100mer = NULL
dat$mappabiliy_40mer = NULL
dat$mappability_24mer = NULL

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

# save
save(dat,datchroms, file=paste0("data/procData/validationData/validationData_SomaMutDB_processed_", tissue, ".RData"))



# Choudhury et al. #####
print("Choudhury")
#load table of predictor files to use
tab = read_xlsx("data/rawdata/dataMapping.xlsx",
                sheet="breast", col_names = T)  # use breast instead of heart?
tab = as.data.frame(tab)
tab[tab == "NA"] = NA

# map predictors
pred = mapPredictors(x=tab,
                     posFile=paste0("data/procData/heart/mutData_validation_Choudhury.bed"))

# save raw mapped data
load(paste0("data/rdata/heart/mutData_validation_Choudhury.RData"))
data = list(meta = tab, pred = pred, muts = Muts)
save(data, file = paste0("data/procData/validationData/validationData_Choudhury.RData"))

cat("\n")

# process data further
dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
datchroms = data$muts$chr

# remove datapoints that could be problematic
toExclude = dat$ConsensusExcludable == 1 | dat$repeatMasker ==1 | dat$tandemRepeatFinder == 1 | 
  dat$mappability_100mer < 1 |dat$mappabiliy_40mer < 1 | dat$mappability_24mer < 1
dat = dat[!toExclude,]
datchroms = datchroms[!toExclude]
dat$ConsensusExcludable = NULL
dat$repeatMasker = NULL
dat$tandemRepeatFinder = NULL
dat$mappability_100mer = NULL
dat$mappabiliy_40mer = NULL
dat$mappability_24mer = NULL

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

# save
save(dat,datchroms, file=paste0("data/procData/validationData/validationData_Choudhury_processed.RData"))
