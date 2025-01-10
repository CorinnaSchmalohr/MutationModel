tissues = c("brain","breast",  "colon", "esophagus", "kidney", "liver", "lung", "ovary", "prostate","skin" )#  
dumpVar = lapply(tissues, function(tissue){
  print(tissue)
  meta = read.table(paste0("data/rawdata/HiCENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$Output.type %in% c("genome compartments",
                                      "mapping quality thresholded contact matrix"),]
  write.table(meta$File.download.URL,
              file = paste0("data/rawdata/HiCENCODE/", tissue, "/filesToDownload.txt"),
              quote = F,row.names = F, col.names = F)
})