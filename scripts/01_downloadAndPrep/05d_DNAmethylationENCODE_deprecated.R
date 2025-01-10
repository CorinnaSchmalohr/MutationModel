dir.create("data/predictors/methylationENCODE/", showWarnings = F)
library(GenomicRanges)
library(rtracklayer)
library(parallel)
ch = import.chain("/cellnet/MutationModel/data/rawdata/genome/hg38ToHg19.over.chain")
chroms = paste0("chr", 1:22)
fileList = read.table("data/rawdata/methylationENCODE/files2Download.csv", sep = ",")
fileList$id = do.call(rbind,strsplit(fileList$V2, split = "[/.]"))[,9]
dumpVar = sapply(unique(fileList$V1), function(tissue){
  print(tissue)
  ids = fileList[fileList$V1 == tissue, 3]
  GRL = lapply(ids, function(id){
    print(id)
    cmd = paste0("zcat data/rawdata/methylationENCODE/", id, ".bed.gz  | ",
                 "awk '($10 >= 6 && $14 >= 20){print $1,$2,$11}'")
    temp = read.table(text = system(cmd, intern = T), sep = " ")
    temp = temp[temp$V1 %in% chroms,]
    gr = GRanges(seqnames=temp$V1,
                 ranges=IRanges(start=temp$V2+1, width = 1), 
                 score = temp$V3)
    gr = liftOver(gr, ch)
    gr = unlist(gr[elementNROWS(gr)==1])
    return(gr)
  })
  if(length(GRL) < 2){
    gr = sort(GRL[[1]])
    save(gr, file = paste0("data/predictors/methylationENCODE/",
                           tissue,".RData")) 
  }
  # GRL = GRangesList(GRL)
  names(GRL) = ids
  gr <- do.call(c, c(GRL, use.names = FALSE))
  mcols(gr) <- NULL
  gr <- unique(sort(gr))
  # get dataframe of signal value for each dataset in the output GRanges
  idx <- mclapply(GRL, function(x) which(gr %in% x), mc.cores = 8)
  counts <- mcmapply(function(dat, idx) {
    out <- rep.int(NA, length(gr))
    out[idx] <- mcols(dat)[["score"]]
    out
  }, GRL, idx, mc.cores = 8, SIMPLIFY = TRUE)
  gr$score = rowMeans(counts, na.rm = T)
  save(gr, file = paste0("data/predictors/methylationENCODE/",
                          tissue,".RData")) 
})



