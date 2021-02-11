tissues = c( "colon", "luad", "breast", "skin", # "ovary",
             "kidney", "prostate")
sizes = c("bins10kb", "bins100kb", "bins1Mb") # "bins100bp", "bins1kb",
for(tissue in tissues){
   print(tissue)
   for(size in sizes){
      print(size)
      methylation = read.table(paste0("data/procData/",
                                      tissue,"/methylation/bins/", 
                                      size, "_methylation.out.bed"),
                               as.is = T,sep = "\t", na.strings="NAN")[,5]
      methylation[is.na(methylation)] = 0
      save(methylation,
           file = paste0("data/procData/", 
                         tissue, "/methylation/bins/",
                         size, "_methylation.RData"))
   }
}
