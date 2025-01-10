.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
source("lib/dataMapping.R")
library(readxl)
dir.create("data/procData/validationData/", showWarnings = F)

# SomaMutDB #####
print("SomaMutDB")
tissues = c("brain","breast", "colon","esophagus", 
            "kidney", "liver", "lung","prostate", "skin")
for(tissue in tissues){
  print(tissue)
  # load table of predictor files to use
  tab = read_xlsx("data/rawdata/dataMapping.xlsx", 
                  sheet=tissue, col_names=T)
  tab$NA. = NULL
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
  
  # save
  save(dat,datchroms, file=paste0("data/procData/validationData/validationData_SomaMutDB_processed_", tissue, ".RData"))
}
#####
