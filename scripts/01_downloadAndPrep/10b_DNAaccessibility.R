# Nucleosome from UCSC
# use raw bigWig file.

# DNAse from UCSC, from multiple cell lines combined #####
library(GenomicRanges)
DNAse = read.table("data/rawdata/UCSC_tracks/wgEncodeRegDnaseClusteredV3.bed.gz", 
                   as.is = T)
colnames(DNAse) = c("chr", "start", "end", "name", "score", "nCelltypes", "whichCelltypes", "CellScores")
gr = GRanges(seqnames=DNAse$chr, 
             ranges=IRanges(start = DNAse$start+1, end=DNAse$end), 
             score = DNAse$score)
gr = gr[grep("[YXM]", seqnames(gr), invert=T)]
save(gr, file = "data/predictors/wgEncodeRegDnaseClusteredV3.RData")
##### 


# DNase-seq from ENCODE ##### 
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
dir.create("data/predictors/DNAseENCODE", showWarnings = F)
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  meta = read.table(paste0("data/rawdata/DNAseENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T, quote = "")
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
  #   print("more than one bed replicate present")
  # }
  
  # narrowPeak: chrom, start, end, name, score, strand, signal value, pValue, qValue, peakPosition
  urls = meta$File.download.URL
  write.table(urls,file = paste0("data/rawdata/DNAseENCODE/", tissue, "/filesToDownload.txt"),
              quote = F,row.names = F, col.names = F)
  
  # prepare code for the bigWigs
  files = meta$File.accession[meta$File.format == "bigWig"]
  toUnify = paste0("data/rawdata/DNAseENCODE/", tissue, "/",
                   files, ".bigWig")
  printCmd = paste("echo",  length(toUnify), "files", sep = " ")
  if(length(toUnify) < 2){
    wigToolsCommand = paste0("wiggletools write_bg -  ",
                             toUnify,
                             " > temp/DNAse_", tissue, ".bedGraph")
  } else{
    wigToolsCommand = paste0("wiggletools write_bg -  mean ",
                             paste0(toUnify, collapse = " "),
                             " > temp/DNAse_", tissue, ".bedGraph")
  }
  liftOverCommand= paste0("./lib/liftOver ",
                          "-bedPlus=4 temp/DNAse_", tissue, ".bedGraph ",
                          "data/rawdata/genome/hg38ToHg19.over.chain ",
                          "temp/DNAse_", tissue, "_liftOver.bedGraph ",
                          "temp/DNAse_", tissue, "_unmapped.bed")
  sortCommand = paste0("sort -k1,1 -k2,2n --buffer-size=60G ",
                       "--temporary-directory=temp --parallel=16 ",
                       "temp/DNAse_", tissue, "_liftOver.bedGraph " ,
                       "> temp/DNAse_", tissue, "_liftOver_sorted.bedGraph")
  fixOverlapping = paste0("bedtools merge  -c 4 -o mean -d -1 ",
                          "-i temp/DNAse_", tissue, "_liftOver_sorted.bedGraph ",
                          "> temp/DNAse_", tissue, "_liftOver_sorted_merged.bedGraph ")
  toBWCommand = paste0("lib/bedGraphToBigWig ",
                       "temp/DNAse_", tissue, "_liftOver_sorted_merged.bedGraph ",
                       "data/rawdata/genome/hg19.chrom.sizes ",
                       "data/predictors/DNAseENCODE/", tissue, ".bigWig")
  rmCommand = paste0("rm ",
                     " temp/DNAse_", tissue, ".bedGraph",
                     " temp/DNAse_", tissue, "_liftOver.bedGraph",
                     " temp/DNAse_", tissue, "_liftOver_sorted.bedGraph",
                     " temp/DNAse_", tissue, "_unmapped.bed",
                     " temp/DNAse_", tissue, "_liftOver_sorted_merged.bedGraph")
  cmd = c(wigToolsCommand, liftOverCommand, sortCommand, 
          fixOverlapping, toBWCommand, rmCommand)
  command = c(printCmd,cmd)
  write.table(command,
              file = paste0("data/rawdata/DNAseENCODE/", tissue,
                            "/UnifyCommand.sh"),
              quote = F,row.names = F, col.names = F)
})


##### 