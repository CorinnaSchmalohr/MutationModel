args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
t = tissues[args]
print(t)
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")


source("lib/dataMapping.R")
library(readxl)


ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0) # Range to the left and right around the pos


# Map data over various ranges
dumpVar =   lapply(seq(length(ranges)), function(r){
  print(ranges[r])
  
  # preprocess data table for different ranges of predictors 
  tab = read_xlsx("data/rawdata/dataMapping_RangeAnalysis.xlsx", 
                  sheet=t, col_names = T)
  features_noRange = tab$abbreviation[tab$range != "x"]
  
  
  tab = tab[tab$range == "x",]# Only map the features with varying ranges 
  features = tab$abbreviation
  tab$abbreviation = paste(tab$abbreviation, names(ranges[r]), sep = "_")
  tab$Name = paste(tab$Name, names(ranges[r]), sep = " ")
  tab$range = as.integer(ranges[r])
  
  tab$NA. = NULL
  tab[tab == "NA"] = NA
  
  
  # map predictors 
  print("Mapping predictors")
  pred = mapPredictors(x=tab, 
                       posFile=paste0("data/procData/", t, "/Muts.bed"))
  load(paste0("data/rdata/", t, "/Muts.RData"))
  data = list(meta = tab, pred = pred, muts = Muts)
  
  
  # preprocessing further"
  print("preprocessing further")
  
  dat_bin = cbind(data$pred, mutated = as.factor(data$muts$mutated))
  datchroms_bin = data$muts$chr
  
  
  # center and log and replace missing values with mean:
  dat_bin = sapply(colnames(dat_bin),FUN=function(j){ 
    x = dat_bin[,j]
    if(!is.numeric(x)){return(x)}
    if((max(x, na.rm = T)-min(x, na.rm = T)) != 0){
      x = (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)) # scale from 0 to 1
    }
    x = log(x +1) # log
    x[is.na(x)] = mean(x,na.rm = T) # replace missing
    return(x) 
  }, simplify=F)
  dat_bin = as.data.frame(dat_bin)
  
  # Merge mapped features w/o ranges and the new mapped features with a different range
  load(file=paste0("data/procData/traindata/traindata_processed_",
                   t, ".RData"))
  dat_noRange = dat[,na.omit(match(features_noRange,colnames(dat)))]
  
  dat_bin = cbind(dat_bin, dat_noRange)
  
  save(dat_bin,datchroms_bin,features, 
       file=paste0("data/procData/traindata/traindata_processed_",
                   t, "_BinAnalysis_",names(ranges[r]),".RData"))
})
