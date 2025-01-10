library(GenomicRanges)
load("data/processedData/chrLengths.RData")
dir.create("data/predictors/nonB", showWarnings=F)
#get a list of non-B DNA structures available
types = list.files("data/rawdata/non-b_db_human_hg19_tsv")
# iterate through structures 
# load, create a GRanges object and save
for(type in types){
   print(type)
   files = list.files(paste0("data/rawdata/non-b_db_human_hg19_tsv/", 
                     type, "/tsv"), pattern="chr[0-9]", full.names=T)
   gr = do.call(c,lapply(files, function(file){
      tab = read.table(file, header = T, stringsAsFactors=F)
      rang = GRanges(seqnames=tab$Sequence_name, 
                     ranges=IRanges(start=tab$Start, end=tab$Stop),
                     seqlengths=chrLengths)
      return(rang)
   }))
   save(gr, file = paste0("data/predictors/nonB/", type, ".RData"))
}

