tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
dir.create("data/predictors/histoneENCODE/", showWarnings = F)
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  meta = read.table(paste0("data/rawdata/histoneENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$File.type %in% c("bigWig", "bed") , ]
  
  # bW = meta[meta$File.type == "bigWig" , ]
  # bed = meta[meta$File.type == "bed", ]
  # if(length(table(table(bW$Experiment.accession)))>1){
  #   print("more than one bW replicate present")
  # }
  # if(length(table(table(bed$Experiment.accession)))>1){
  #   print("more than one bedreplicate present")
  # }
  
  # narrowPeak: chrom, start, end, name, score, strand, signal value, pValue, qValue, peakPosition
  urls = meta$File.download.URL
  write.table(urls,file = paste0("data/rawdata/histoneENCODE/", tissue, "/filesToDownload.txt"),
              quote = F,row.names = F, col.names = F)
  
  # prepare code for the bigWigs
  filesByTargetBW = split(meta$File.accession[meta$File.format == "bigWig"],
                          meta$Experiment.target[meta$File.format == "bigWig"])
  command = do.call(c,lapply(names(filesByTargetBW), function(target){
    toUnify = paste0("data/rawdata/histoneENCODE/", tissue, "/",
                     filesByTargetBW[[target]], ".bigWig")
    printCmd = paste("echo", target, length(toUnify), "files", sep = " ")
    if(length(toUnify) < 2){
      wigToolsCommand = paste0("wiggletools write_bg -  ",
                               toUnify,
                               " > temp/", target, "_", tissue, ".bedGraph")
      # cmd = paste0("cp ", toUnify,
      #              " data/predictors/histoneENCODE/", target, "_", tissue, ".bigWig")
    } else{
      wigToolsCommand = paste0("wiggletools write_bg -  mean ",
                               paste0(toUnify, collapse = " "),
                               " > temp/", target, "_", tissue, ".bedGraph")
    }
    liftOverCommand= paste0("./lib/liftOver ",
                            "-bedPlus=4 temp/", target, "_", tissue, ".bedGraph ",
                            "data/rawdata/genome/hg38ToHg19.over.chain ",
                            "temp/", target, "_", tissue, "_liftOver.bedGraph ",
                            "temp/", target, "_", tissue, "_unmapped.bed")
    sortCommand = paste0("sort -k1,1 -k2,2n --buffer-size=60G ",
                         "--temporary-directory=temp --parallel=16 ",
                         "temp/", target, "_", tissue, "_liftOver.bedGraph " ,
                         "> temp/", target, "_", tissue, "_liftOver_sorted.bedGraph")
    fixOverlapping = paste0("bedtools merge  -c 4 -o mean -d -1 ",
                            "-i temp/", target, "_", tissue, "_liftOver_sorted.bedGraph ",
                            "> temp/", target, "_", tissue, "_liftOver_sorted_merged.bedGraph ")
    toBWCommand = paste0("lib/bedGraphToBigWig ",
                         "temp/", target, "_", tissue, "_liftOver_sorted_merged.bedGraph ",
                         "data/rawdata/genome/hg19.chrom.sizes ",
                         "data/predictors/histoneENCODE/", target, "_", tissue, ".bigWig")
    rmCommand = paste0("rm ",
                       " temp/", target, "_", tissue, ".bedGraph",
                       " temp/", target, "_", tissue, "_liftOver.bedGraph",
                       " temp/", target, "_", tissue, "_liftOver_sorted.bedGraph",
                       " temp/", target, "_", tissue, "_unmapped.bed",
                       " temp/", target, "_", tissue, "_liftOver_sorted_merged.bedGraph")
    cmd = c(wigToolsCommand, liftOverCommand, sortCommand, 
            fixOverlapping, toBWCommand, rmCommand)
    return(c(printCmd,cmd))
  }))
  write.table(command,
              file = paste0("data/rawdata/histoneENCODE/", tissue,
                            "/UnifyCommand.sh"),
              quote = F,row.names = F, col.names = F)
  
})


