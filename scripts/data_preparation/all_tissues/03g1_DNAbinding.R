# give 1 to 7 to get the different tissues
tissue = commandArgs(trailingOnly=T)
# args = as.numeric(commandArgs(trailingOnly=T))
# tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
# tissue = tissues[args]
print(tissue)


dir.create(paste0("data/procData/", tissue, "/DNAbinding/"), showWarnings=F)
meta = read.table(paste0("data/rawdata/DNAbinding/", tissue, "/metadata.tsv"),
                         header = T, sep = "\t", as.is = T)
fc = meta[meta$Output.type %in% c("fold change over control", "signal") & 
             meta$Assembly == "hg19" & 
             meta$File.Status == "released" , ]
fc = fc[grep(x = fc$Audit.ERROR,pattern = "extreme",invert = T),]
fc = fc[grep(x = fc$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
fc = fc[grep(x = fc$Audit.NOT_COMPLIANT,pattern = "insufficient",invert = T),]

experiments = unique(fc$Experiment.accession)
urls = sapply(experiments, function(i){
   temp = fc[fc$Experiment.accession == i,]
   if(nrow(temp) == 1){
      return(temp$File.download.URL)
   } else {
      return(temp$File.download.URL[temp$Biological.replicate.s. == "1, 2"])
   }
})
write.table(urls,file = paste0("data/rawdata/DNAbinding/", tissue, "/filesToDownload.txt"), 
            quote = F,row.names = F, col.names = F)
