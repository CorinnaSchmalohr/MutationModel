# load all files from replicationdomain, scale and then save again with appendix "_scaled"
files = list.files("data/rawdata/replicationdomain/lung")
# files = files[grep("_scaled", files, invert = T, fixed = T)]
files = files[files != "metadata.tsv"]
dir.create("data/procData/lung/replicationdomain")
for(i in files){
   print(i)
   t = read.table(paste0("data/rawdata/replicationdomain/lung/",i), stringsAsFactors = F)
   if(length(table(t[,1]))<2){
      files = files[files != i]
      print(paste0(i, " was removed"))
      next
   }
   v = t[,4]
   v[is.na(v)] = 0
   v2 = v / (sqrt(sum(v^2)/(length(v)-1)))
   t[,4] = v2
   assign(i,t)
   write.table(t,file = paste0("data/procData/lung/replicationdomain/",i,"_scaled"), col.names = F, quote = F, row.names = F, sep = "\t")
}

