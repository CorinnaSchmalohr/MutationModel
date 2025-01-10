library(rtracklayer)
library(GenomicRanges)
library(matrixStats)
ch = import.chain("data/rawdata/genome/hg38ToHg19.over.chain")
tissues = c("brain","breast", "colon","esophagus", 
            "kidney", "liver", "lung","ovary", 
            "prostate", "skin") 

# load annotation #####
#450k
annot450k = read.table("data/rawdata/methbank/humanmethylation450_15017482_v1-2.csv", sep = ",", 
                       skip = 7, quote = "", fill = T, header = T)
annot450k = annot450k[annot450k$CHR %in% as.character(1:22),c("Name","CHR", "MAPINFO"),] #hg37
gr450k = GRanges(seqnames=paste0("chr",annot450k$CHR), 
             ranges=IRanges(start=annot450k$MAPINFO,width = 1))
names(gr450k) = annot450k$Name
# 850k
annot850k = read.table("data/rawdata/methbank/MethylationEPIC v2.0 Files/EPIC-8v2-0_A1.csv", sep = ",", 
                       skip = 7, quot = "", fill = T, header = T)
annot850k = annot850k[annot850k$CHR %in% paste0("chr",1:22),c("Name", "CHR", "MAPINFO")] #hg38
annot850k = unique(annot850k)

annot850k2 = read.table("data/rawdata/methbank/MethylationEPIC v2.0 Files/Removed probes_EPICv1 to EPICv2.csv", sep = ",", header = T, quote = "")
annot850k2 = annot850k2[,c("EPICv1_Loci", "CHR", "MAPINFO")] #hg38
annot850k2 = unique(annot850k2)
colnames(annot850k2) = c("Name", "CHR", "MAPINFO")
annot850k2 = annot850k2[!annot850k2$Name %in% annot850k$Name, ]

gr850k = GRanges(seqnames=c(annot850k$CHR, annot850k2$CHR), 
             ranges=IRanges(start=c(annot850k$MAPINFO, annot850k2$MAPINFO),width = 1))
names(gr850k) = c(annot850k$Name, annot850k2$Name)
genome(gr850k) = "hg38"
gr850k = liftOver(gr850k, ch)
gr850k = unlist(gr850k)
genome(gr850k) = "hg19"
gr850k = sort(gr850k)

annot = c(gr450k[!names(gr450k) %in% names(gr850k)], gr850k)
# missing probes in annot correspond to control probes, so no need to try to fix those

rm(annot450k, annot850k, annot850k2, gr450k, gr850k, ch)
#####


# load samples
fileEndings = c(".txt.gz","_F.txt.gz", "_M.txt.gz", 
                "_850K.txt.gz",  "_850K_F.txt.gz", "_850K_M.txt.gz")
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  data = lapply(fileEndings, function(f){
    if(!file.exists(paste0("data/rawdata/methbank/", tissue, f))){
      return(NULL)
    }
    dat = read.table(paste0("data/rawdata/methbank/", tissue, f), header = T)
    dat = dat[dat$probe %in% names(annot),]
    n = setNames((ncol(dat)-1)-(rowCounts(as.matrix(dat[,(-1)]), value = NA)), dat$probe) # get the number of non-na samples per probe
    dat = setNames(rowSums(dat[,(-1)], na.rm= T), dat$probe) # add all non-na values per probe
    # map to positions in annot
    list(n = n[names(annot)], res= dat[names(annot)])
  })
  data = data[!sapply(data, is.null)]
  if(length(data) <1){return(NULL)}
  # manually calculate mean
  totN = rowSums(sapply(data, function(x){x$n}), na.rm = T)
  totSum = rowSums(sapply(data, function(x){x$res}), na.rm = T)
  totMean = totSum/totN
  gr = annot
  gr$val = totMean
  gr = sort(gr)
  save(gr, file = paste0("data/predictors/methylation/", 
                         tissue, ".RData"))
})

