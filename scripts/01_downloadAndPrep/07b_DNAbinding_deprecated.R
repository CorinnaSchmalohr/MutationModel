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
   write.table(urls,file = paste0("data/rawdata/DNAbinding/", tissue, "/filesToDownload.txt"), 
               quote = F,row.names = F, col.names = F)
   
   files = unname(do.call(rbind,strsplit(urls, split="/"))[,5])
   
   rownames(fc) = fc$File.accession
   targets = fc[files,"Experiment.target"]
   filesByTarget = split(files, targets)
   command = do.call(c,lapply(names(filesByTarget), function(target){
      toUnify = paste0("data/rawdata/DNAbinding/", tissue, "/",
                       filesByTarget[[target]], ".bigWig")
      printCmd = paste("echo", target, length(toUnify), "files", sep = " ")
      if(length(toUnify) < 2){
         cmd = paste0("cp ",
                      toUnify,
                      " data/predictors/DNAbinding/", target, "_", tissue, ".bigWig")
      } else{
         wigToolsCommand = paste0("wiggletools write_bg - mean \\\n",
                                  paste0("scale $(wiggletools stddevI ", toUnify,
                                  " | awk '{print 1/$1}') \\\n", toUnify, " \\\n", collapse = ""),
                                  "> temp/", target, "_", tissue, ".bedGraph")
         sortCommand = paste0("sort -k1,1 -k2,2n --buffer-size=60G \\\n",
                              "--temporary-directory=temp --parallel=16 \\\n",
                              "temp/", target, "_", tissue, ".bedGraph \\\n" ,
                              "> temp/", target, "_", tissue, "_sorted.bedGraph")
         toBWCommand = paste0("lib/bedGraphToBigWig \\\n",
                              "temp/", target, "_", tissue, "_sorted.bedGraph \\\n",
                              "data/rawdata/hg19.chrom.sizes \\\n",
                              "data/predictors/DNAbinding/", target, "_", tissue, ".bigWig")
         rmCommand = paste0("rm ",
                            " temp/", target, "_", tissue, ".bedGraph",
                            " temp/", target, "_", tissue, "_sorted.bedGraph")
         cmd = c(wigToolsCommand, sortCommand, toBWCommand, rmCommand)
      }
      return(c(printCmd,cmd))
   }))
   write.table(command,
               file = paste0("data/rawdata/DNAbinding/", tissue,
                             "/UnifyCommand.sh"), 
               quote = F,row.names = F, col.names = F)
})


