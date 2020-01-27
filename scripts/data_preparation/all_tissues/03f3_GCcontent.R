# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


GCcontent = sapply(c("5bp", "10bp", "100bp", "1kb"), #, "1Mb"
                   function(f){
   print(f)
   temp = read.table(paste0("data/procData/", tissue, 
                         "/GCcontent/MutsBed",f,".out"),
                  as.is = T)
   ref = read.table(paste0("data/procData/", tissue,
                           "/GCcontent/MutsBed",f,".bed"), as.is = T)
   missing = which(!1:nrow(ref) %in% temp[,2])
   if (length(missing) > 0) {
      temp = rbind(temp,cbind(0,which(!1:nrow(ref) %in% temp[,2])))
   }
   temp = temp[order(temp[,2]),]
   ref = read.table(paste0("data/procData/", tissue,
                           "/GCcontent/MutsBed",f,".bed"), as.is = T)
   temp[,1]/(ref[,3] - ref[,2])
})
colnames(GCcontent) = paste0("GCcontent_",
                             c("5bp", "10bp", "100bp", "1kb")) #, "1Mb"
save(GCcontent, file = paste0("data/rdata/", tissue, 
                              "/GCcontent.RData"))
#####