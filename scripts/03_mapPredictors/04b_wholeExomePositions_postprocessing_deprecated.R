args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
tissue = tissues[args]
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")

load(paste0("data/procData/traindata/traindata_processed_", 
            tissue, ".RData"))
x = colnames(dat)
return(x)

for(i in 1:31){
  print(i)
  
  # list of meta, pred, muts
  for(tissue in tissues){
    cat(tissue, ' ')
    load(paste0("data/procData/exomeMuts/exomeMuts_part",i,
                "_", tissue, ".RData")) # data
    
    dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
    datpos = as.data.frame(cbind(chr = data$muts$chr, context = data$muts$context))
    
    # 6. center and log and replace missing values with mean:
    dat = sapply(colnames(dat),FUN=function(j){
      x = dat[,j]
      if(!is.numeric(x)){return(x)}
      if((max(x, na.rm = T)-min(x, na.rm = T)) != 0){
        x = (x-min(x, na.rm = T))/(max(x, na.rm = T)-min(x, na.rm = T)) # scale from 0 to 1
      }
      x = log(x +1) # log
      x[is.na(x)] = mean(x,na.rm = T) # replace missing
      # x = scale(x)
      return(x)
    }, simplify=F)
    dat = as.data.frame(dat)
    
    
    save(dat,datpos, file=paste0("data/procData/exomeMuts/exomeMuts_part",i,"_",
                                 tissue,"_processed.RData"))
  }
}


