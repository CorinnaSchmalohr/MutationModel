tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
library(GenomicRanges)
library(rtracklayer)
# library(BRGenomics)
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  meta = read.table(paste0("data/rawdata/histoneENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$File.type == "bed", ]
  filesByTargetBed = split(meta$File.accession,
                          meta$Experiment.target)
  dumpVar2 = sapply(names(filesByTargetBed), function(target){
    cat(target, ' ')
    ids = filesByTargetBed[[target]]
    if(length(ids)>1){
      bedToBedgraph = paste0("zcat data/rawdata/histoneENCODE/", 
                             tissue,"/", ids, ".bed.gz", 
                             " | cut -f1-3,8 | sort -k 1,1 -k2,2n ",
                             " > temp/", ids, ".bedgraph")
      fixOverlapping = paste0("bedtools merge  -c 4 -o mean -d -1 ",
                              "-i temp/", ids, ".bedgraph ",
                              "> temp/", ids, "_merged.bedgraph ")
      cmd=paste(c(bedToBedgraph, fixOverlapping), collapse = "\n")
      system(cmd)
      unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                          paste(" temp/", ids, "_merged.bedgraph", sep = "", collapse = " "))
      comb = read.table(text = system(unionbedg, intern = T))  
      file.remove(c(paste0("temp/", ids, ".bedgraph"), paste0("temp/", ids, "_merged.bedgraph")))
      comb = comb[comb$V1 %in% paste0("chr", 1:22),]
      val = apply(comb[,-(1:3)],1,function(x){sum(!is.na(x))})/(ncol(comb)-3)
      gr = GRanges(seqnames=comb$V1,
                   ranges=IRanges(start=comb$V2+1, end = comb$V3), 
                   val = val)
    } else{
      comb = read.table(paste0("data/rawdata/histoneENCODE/",tissue,"/", ids, ".bed.gz")) 
      comb = comb[comb$V1 %in% paste0("chr", 1:22),]
      gr = GRanges(seqnames=comb$V1,
                   ranges=IRanges(start=comb$V2+1, end = comb$V3), 
                   val = 1)
    }
    ch = import.chain("data/rawdata/genome/hg38ToHg19.over.chain")
    gr = liftOver(gr, ch)
    gr = unlist(gr[elementNROWS(gr)==1])
    gr = sort(gr)
    save(gr, file=paste0("data/predictors/histoneENCODE/",
                         target, "_", tissue,".RData"))
  }); cat('\n')
})