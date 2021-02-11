# cd data/rawdata
# wget -O non-b_db_human_hg19_tsv.tar.gz \
#    https://ncifrederick.cancer.gov/bids/ftp/actions/download/?resource=/bioinfo/static/nonb_dwnld/human_hg19/human_hg19.tsv.tar.gz
# tar -xz --one-top-level -f non-b_db_human_hg19_tsv.tar.gz
# cd ../..


library(GenomicRanges)
# set chromosome lengths
chrs = sort(paste0("chr", c(1:22)))
chrLengths = read.table("data/rawdata/GRCh37.p11.genome.fa.fai")
chrLengths = setNames(chrLengths[,2],chrLengths[,1])
chrLengths = chrLengths[chrs]

structures = list.files("data/rawdata/non-b_db_human_hg19_tsv")

apr = do.call(c,lapply(chrs, function(cr){
   temp = read.table(paste0("data/rawdata/non-b_db_human_hg19_tsv/a-phased_repeats/tsv/", 
                            cr, "_APR.tsv"), header = T, stringsAsFactors=F)
   gr = GRanges(seqnames=temp$Sequence_name, 
                ranges=IRanges(start=temp$Start, end=temp$Stop),
                seqlengths=chrLengths)
   return(gr)
}))

dr = do.call(c,lapply(chrs, function(cr){
   temp = read.table(paste0("data/rawdata/non-b_db_human_hg19_tsv/direct_repeats/tsv/", 
                            cr, "_DR.tsv"), header = T, stringsAsFactors=F)
   gr = GRanges(seqnames=temp$Sequence_name, 
                ranges=IRanges(start=temp$Start, end=temp$Stop),
                seqlengths=chrLengths)
   return(gr)
}))

# this is how to find the values for the mut positions:
# load("data/rdata/luad/Muts.RData")
# MutsR = GRanges(seqnames=Muts$chr, ranges=IRanges(Muts$pos, width=1))
# res = findOverlaps(query=MutsR, subject=grl, select = "first")
# res = as.integer(!is.na(findOverlaps(query=MutsR, subject=apr, select = "first")))
# temp = read.table(paste0("data/procData/luad/structures/aPhasedRepeats.tsv"),
#            sep = "\t", stringsAsFactors = F)

