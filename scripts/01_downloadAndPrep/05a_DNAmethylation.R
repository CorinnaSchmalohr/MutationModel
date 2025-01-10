# Data from this paper https://www.nature.com/articles/s41586-022-05580-6#MOESM1
# downloaded here: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186458
library(rtracklayer)
library(parallel)
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
# run only once:
# untarCmd = paste0("mkdir data/rawdata/methylationGSE186458/", tissues, "; ",
#                   "tar -xf  data/rawdata/methylationGSE186458/", tissues, 
#                   ".tar --directory  data/rawdata/methylationGSE186458/", tissues, "; ", collapse = "")
# system(untarCmd)
# dir.create("data/predictors/methylationGSE186458/")

dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  files = list.files(paste0("data/rawdata/methylationGSE186458/", tissue))
  # iterate over chromosomes, for memory reasons
  GRL = sapply(paste0("chr", 1:22), function(chrom){
    cat(chrom, ' ')
    GRLchrom  = lapply(files, function(x){
      gr = import(paste0("data/rawdata/methylationGSE186458/", tissue, "/", x))
      gr = gr[seqnames(gr) == chrom]
      start(gr) = start(gr)+1
      return(gr)
    }) 
    names(GRLchrom) = files
    gr <- do.call(c, c(GRLchrom, use.names = FALSE))
    mcols(gr) <- NULL
    gr <- unique(sort(gr))
    # get dataframe of signal value for each dataset in the output GRanges
    idx <- mclapply(GRLchrom, function(x) which(gr %in% x), mc.cores = 4)
    counts <- mcmapply(function(dat, idx) {
      out <- rep.int(NA, length(gr))
      out[idx] <- mcols(dat)[["score"]]
      out
    }, GRLchrom, idx, mc.cores = 8, SIMPLIFY = TRUE)
    gr$score = rowMeans(counts, na.rm = T)
    return(gr)
  })
  gr = do.call(c, c(GRL, use.names = FALSE))
  save(gr, file = paste0("data/predictors/methylationGSE186458/", tissue, ".RData"))
})

