.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
source("lib/dataMapping.R")
library(readxl)
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)


# prepare mapping table
print("mapping data")
tab = read_xlsx("data/rawdata/dataMappingAlltissues.xlsx", 
                sheet="allTissues", col_names=T)
tab$NA. = NULL
tab[tab == "NA"] = NA
tab = tab[,c(colnames(tab)[1:9],"allTissues")]
tab = tab[!is.na(tab[,"allTissues"]),]
# for predictors where we want multiple ranges, expand table
tab = apply(tab,1,function(x){
  # print( x[["Name"]])
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
#####


# iterate through the three parts that the data was saved in and map them
dumpVar = sapply(1:3, function(i){
  cat("part ", i, "\n")
  if(file.exists(paste0("data/MutTables/exomeTrainData/allTissues_Muts_part", i, "_mapped.RData"))){
    return(NA)
  }
  pred = mapPredictors(x=tab, 
                       posFile=paste0("data/MutTables/exomeTrainData/allTissues_Muts_part", i, ".bed"))
  load(paste0("data/MutTables/exomeTrainData/allTissues_Muts_part", i,".RData"))
  tempData = list(meta = tab, pred = pred, muts = Muts)
  save(tempData, file = paste0("data/MutTables/exomeTrainData/allTissues_Muts_part", i, "_mapped.RData"))
  rm(pred, Muts, tempData); gc()
  return(NA)
})
dataList = sapply(1:3, function(i){
  load(paste0("data/MutTables/exomeTrainData/allTissues_Muts_part", i, "_mapped.RData"))
  return(tempData)
}, simplify = F)
# combine the parts
data = list(meta = do.call(rbind, sapply(dataList, function(x)  {x$meta} , simplify = F)),
            pred = do.call(rbind, sapply(dataList, function(x)  {x$pred}  , simplify = F)), 
            muts = do.call(rbind, sapply(dataList, function(x)  {x$muts}  , simplify = F)))
save(data, file = paste0("data/MutTables/exomeTrainData/allTissues_Muts_mapped.RData"))
dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
datchroms = data$muts$chr
dat = as.data.frame(dat)
save(dat,datchroms, file=paste0("data/MutTables/exomeTrainData/allTissues_Muts_mapped_processed.RData"))
cat("\n")

