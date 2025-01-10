.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
source("lib/dataMapping.R")
library(readxl)
tissues = c("brain","breast", "colon","kidney", "liver", "lung", #"esophagus",
            "ovary", "prostate", "skin")
# for each tissue, prepare corresponding data
print("mapping data")
for(tissue in tissues){
   print(tissue)
  tab = read_xlsx("data/rawdata/dataMapping.xlsx", 
                  sheet=tissue, col_names=T)
   tab$NA. = NULL
   tab[tab == "NA"] = NA
   pred = mapPredictors(x=tab, 
                        posFile=paste0("data/procData/", tissue,
                                       "/Muts_subsampledNsamples.bed"))
   load(paste0("data/rdata/", tissue, "/Muts_subsampledNsamples.RData"))
   data = list(meta = tab, pred = pred, muts = Muts)
   save(data, file = paste0("data/procData/traindata/traindata_", tissue,
                            "_subsampledNsamples.RData"))
   
   cat("\n")
   
   # process data further
   dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
   datchroms = data$muts$chr
   
   
   # center and log and replace missing values with mean:
   dat = sapply(colnames(dat),FUN=function(j){ 
      x = dat[,j]
      if(!is.numeric(x)){return(x)}
      if((max(x, na.rm = T)-min(x, na.rm = T)) != 0){
        x = (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)) # scale from 0 to 1
      }
      x = log(x +1) # log
      x[is.na(x)] = mean(x,na.rm = T) # replace missing
      return(x) 
   }, simplify=F)
   dat = as.data.frame(dat)
   
   save(dat,datchroms, file=paste0("data/procData/traindata/traindata_processed_",
                                   tissue, "_subsampledNsamples.RData"))
}

