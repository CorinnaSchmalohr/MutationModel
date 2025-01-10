files =read.table("data/rawdata/replication/files.txt", sep = "\t")
temp = lapply(files$V2,function(x){
   x = strsplit(x,split="; ")[[1]]
   x = do.call(rbind,strsplit(x,split="="))
   setNames(x[,2], x[,1])
})
tags = unique(unlist(sapply(temp, names)))
meta = t(sapply(temp, function(x){
   x[tags]
}))
colnames(meta) = tags
meta = data.frame(fileName = files$V1, meta)
toUse = meta[meta$view %in% c("WaveSignal", "Peaks", "Valleys"),
             c("fileName", "view", "cell", "replicate", 
               "dccAccession", "subId", "labExpId", "type", 
               "md5sum", "geoSampleAccession", "tableName")]
write.table(toUse, file="data/rawdata/replication/metadata.tsv", 
            sep="\t", quote=F, col.names=T, row.names=F)
filesToDownload = paste0("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/",
                         toUse$fileName)
write.table(filesToDownload, file="data/rawdata/replication/filesToDownload.txt", 
            sep="\t", quote=F, col.names=F, row.names=F)
