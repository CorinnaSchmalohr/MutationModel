dir.create("data/predictors/HiCENCODE", showWarnings = F)
library(rtracklayer)
# library(rhdf5)
library(strawr) #remotes::install_github("aidenlab/straw/R")
library(Matrix)
library(GenomicRanges)
ch = import.chain("data/rawdata/genome/hg38ToHg19.over.chain")

# for each tissue
tissues =c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
           "ovary", "prostate", "skin")
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  ## load meta file
  meta = read.table(paste0("data/rawdata/HiCENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[meta$File.Status == "released",]
  
  # interactions
  intMeta  = meta[meta$Output.type == "mapping quality thresholded contact matrix",]
  print("ints")
  GRLints = lapply(intMeta$File.accession, function(id){
    cat(id, ' ')
    fname = paste0("data/rawdata/HiCENCODE/", 
                   tissue, "/", id, ".hic")
    # get chromosomes
    # chroms = readHicChroms(fname)
    # if(!all(paste0("chr", 1:22) %in% chroms$name)){
    #   error("check chromosomes")
    # }
    chroms = paste0("chr", 1:22)
    # get resolutions and normalization
    # res = readHicBpResolutions(fname)
    # norm = readHicNormTypes(fname) # SCALE replaced KR
    # if(!"SCALE" %in% norm){
    #   error("check normalization")
    # }
    GRperChr = lapply(chroms, function(chrom){
      hic = straw(matrix = "oe", # observed, oe
                  norm = "SCALE", 
                  fname = fname, 
                  chr1loc=chrom, chr2loc=chrom, 
                  unit="BP", binsize=50000)
      bins = sort(unique(c(hic$x, hic$y)))
      bins2Index = setNames(seq_len(length(bins)), bins)
      mat = sparseMatrix(i = bins2Index[as.character(hic$x)], 
                         j=bins2Index[as.character(hic$y)], 
                         x=as.numeric(log(hic$counts+1)), 
                         symmetric = T)
      # calculate mean over columns. cannot simply use rowMeans
      # because 0-entries of sparse matrix are used too.
      vals = rowSums(mat, na.rm = T)
      vals = scale(vals)[,1]
      gr = GRanges(seqnames=chrom,
                   ranges=IRanges(start=bins+1,width=50000),
                   score = vals)
      gr = gr[!is.na(score(gr))]
      gr = liftOver(gr, ch)
      gr = unlist(gr[elementNROWS(gr)==1])
      gr = sort(gr)
      return(gr)
    })
    GR = suppressWarnings(do.call(c,GRperChr))
    return(GR)
  }); cat('\n')
  
  GRLints = lapply(GRLints, function(x){
    s = score(x)
    score(x) = (s-min(s))/(max(s)-min(s))
    return(x)
  })
  GRLints = GRangesList(GRLints)
  names(GRLints) = intMeta$File.accession
  save(GRLints, file = paste0("data/predictors/HiCENCODE/", tissue,
                              "_ints.RData"))
  
  ## genome compartments
  print("compartments")
  compMeta = meta[meta$Output.type == "genome compartments",]
  if(nrow(compMeta) >=  1){
    GRLcomp = lapply(compMeta$File.accession, function(id){
      cat(id, ' ')
      gr = import.bw(con = paste0("data/rawdata/HiCENCODE/", tissue, "/",
                                  id, ".bigWig"))
      gr = gr[!is.na(score(gr))]
      gr = liftOver(gr, ch)
      gr = unlist(gr[elementNROWS(gr)==1])
      gr = sort(gr)
      return(gr)
    }); cat('\n')
    GRLcomp = GRangesList(GRLcomp)
    names(GRLcomp) = compMeta$File.accession
    scores = sapply(GRLcomp, function(x)score(x))
    save(GRLcomp, file = paste0("data/predictors/HiCENCODE/", tissue,
                                "_comp.RData"))
  }
  
})



