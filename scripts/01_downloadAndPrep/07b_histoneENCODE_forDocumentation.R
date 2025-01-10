tissues = list.files("data/rawdata/DNAbinding/", full.names=F)
dumpVar = sapply(tissues, function(tissue){
   print(tissue)
   meta = read.table(paste0("data/rawdata/DNAbinding/", tissue, "/metadata.tsv"),
                     header = T, sep = "\t", as.is = T)
   fc = meta[meta$Output.type %in% c("fold change over control", "signal") & 
                meta$File.assembly == "hg19" & 
                meta$File.Status == "released" , ]
   fc = fc[grep(x = fc$Audit.ERROR,pattern = "extreme",invert = T),]
   fc = fc[grep(x = fc$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
   fc = fc[grep(x = fc$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
   
   experiments = unique(fc$Experiment.accession)
   urls = sapply(experiments, function(i){
      temp = fc[fc$Experiment.accession == i,]
      if(nrow(temp) == 1){
         return(temp$File.download.URL)
      } else {
         replLength = sapply(strsplit(temp$Biological.replicate.s., ","), length)
         return(temp$File.download.URL[replLength == max(replLength)])
      }
   })
   
   files = unname(do.call(rbind,strsplit(urls, split="/"))[,5])
   
   rownames(fc) = fc$File.accession
   targets = fc[files,"Experiment.target"]
   filesByTarget = split(files, targets)
   cbind(tissue, 
         target = names(filesByTarget),
         IDs = sapply(filesByTarget, function(x){paste(x, collapse = ", ")} ))
})
source("lib/general_function.R")

targets = table(unlist(sapply(dumpVar, function(x)x[,2])))
TFtargets = names(targets)[targets == 1]
restTargets = names(targets)[targets > 1]

temp = sapply(tissues, function(tissue){
   sapply(restTargets, function(target){
      x = dumpVar[[tissue]]
      x = x[x[,"target"] == target,3]
      if(length(x) <1) {return(NA)} else{return(x)}
   })
})
temp = cbind(target = rownames(temp), temp)
cat.table.redmine(temp)

temp2 = sapply(TFtargets, function(target){
   x = dumpVar[["liver"]]
   x = x[x[,"target"] == target,3]
   if(length(x) <1) {return(NA)} else{return(x)}
})
cat.table.redmine(cbind(TF = names(temp2), IDs = temp2))
