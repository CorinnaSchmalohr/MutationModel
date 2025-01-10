args = as.numeric(commandArgs(trailingOnly=T))
cell2ID = c("cell1" =  "MMEB9_E1_CKDL230009768-1A_H3TG2DSX7_L4.final.snv.vcf",
                  "cell2" =  "MMEB9_E2_CKDL230009768-1A_H3TG2DSX7_L4.final.snv.vcf",
                  "cell3" =  "MMEB9_E3_CKDL230009768-1A_H3TG2DSX7_L4.final.snv.vcf",
                  "cell4" =  "MMEA9_E1_CKDL230009767-1A_H3TG2DSX7_L3.final.snv.vcf",
                  "cell5" =  "MMEA9_E2_CKDL230009768-1A_H3TG2DSX7_L4.final.snv.vcf",
                  "cell6" =  "MMEA9_E3_CKDL230009768-1A_H3TG2DSX7_L4.final.snv.vcf" )
cell = names(cell2ID)[args]
print(cell)

library(GenomicRanges)
library(parallel)
library(rtracklayer)
library(stringr)
setwd("/cellfile/cellnet/MutationModel/")

chrs = paste0("chr", c(1:22))
load("data/processedData/pentamers.RData")
rownames(fivemers) = pent2context[fivemers[,1]]

# load vcf file
muts <- read.table(paste0("data/rawdata/scWGS/mutations/",cell2ID[cell]))
muts <- read.table("data/rawdata/scWGS/mutations/MMEB9_E1_CKDL230009768-1A_H3TG2DSX7_L4.final.snv.vcf")
colnames(muts) <- c("chr","pos","ID","ref","alt","qual","filter","info","format","cell")

# extract values from column "cell001" 
stringr::str_split_fixed

mutssplit = str_split_fixed(muts$cell, ":", 6)

'> unique(cell1$format)
[1] "GT:SO:AD:BI:GQ:PL"'
colnames(mutssplit) = c("GT","SO","AD","BI","GQ","PL")

# Add columns with individual cell values to original data frame
muts <- data.frame(as.matrix(muts),mutssplit)
muts = subset(muts, select = -c(ID, qual, filter, info))

# GT = Genotype
# SO = Somatic mutation
# AD = Allelic depths for ref and alt alleles
# BI = Amplification Bias
# GQ = Genotype Quality
# PL = sequencing noise, amplification artifact, heterozygous SNV and homozygous SNV respectively

# mapping positions from HG38 genome to HG19
ch = import.chain("data/rawdata/genome/hg38ToHg19.over.chain")
muts.gr <- makeGRangesFromDataFrame(muts, keep.extra.columns=T, seqnames.field = "chr", start.field="pos",
                                     end.field="pos")
muts.gr = liftOver(muts.gr, ch)
muts.gr = unlist(muts.gr[elementNROWS(muts.gr)==1])
muts.gr = sort(muts.gr)
muts = as.data.frame(muts.gr)
muts = subset(muts, select = -c(end, width, strand, format, cell001))
colnames(muts)[1] <- "chr"
colnames(muts)[2] <- "pos"

# extract sequence context for positions from genome
TPs_bed = data.frame(muts$chr, as.numeric(muts$pos)-3, as.numeric(muts$pos)+2)
write.table(TPs_bed, file="temp/scTPs_bed.bed", quote=F, col.names=F, row.names=F, sep="\t")
cmd = paste0("bedtools getfasta -fi data/rawdata/GRCh37.primary_assembly.genome.fa ", 
             "-bed temp/scTPs_bed.bed ",
             "-name -tab | cut -f 2")

context = system(cmd, intern = T)
muts$context = context

# extract positions to exclude during TN sampling
posToFilter = cbind(muts$chr,
                    paste(muts$chr, muts$pos, sep = "_"))

###  filter TNs for callable regions --> from coverage bed files
# import coverage bed files into one GRanges object
pathToChromFiles <-"data/rawdata/scWGS/coverage/MMEA9_E1_CKDL230009767-1A_H3TG2DSX7_L3/MMEA9_E1_CKDL230009767-1A_H3TG2DSX7_L3."
chrList_bed = lapply(chrs, function(cr){
  x = read.table(paste0(pathToChromFiles, cr, ".bed"),
                 col.names = c("chr","start","end"))
  return(x)
})

chrList_bed = do.call("rbind", chrList_bed)
coverage.gr = makeGRangesFromDataFrame(chrList_bed, keep.extra.columns=T, seqnames.field = "chr", start.field="start",
                                       end.field="end")

# mapping positions from HG38 genome to HG19 for coverage files
coverage.gr = liftOver(coverage.gr, ch)
coverage.gr = unlist(coverage.gr[elementNROWS(coverage.gr)==1])
coverage.gr = sort(coverage.gr)

### get TNs #####

# sample from these positions. 
print("getting TNs")
cl <- makeCluster(8, type="FORK")
TNs = parLapply(cl = cl, unique(pent2context[muts$context]), function(pent){
  cat(pent, ' ')
  pent2search = fivemers[pent,]
  samplePos = do.call(c,lapply(pent2search, function(p){
    load(paste0("data/processedData/pentLocations/", p, ".RData"))
    ## filter pentLoc for coverage:
    # asked ChatGTP: In this code, subsetByOverlaps function is used to subset the pentLoc object 
    # using the intervals in the coverage.gr object. It returns a new GRanges object, 
    # which contains only the values from pentLoc that are covered by the intervals in coverage.gr.
  #  pentLoc = subsetByOverlaps(pentLoc, coverage.gr)
    pentLoc$context = p
    return(pentLoc)
    
  }))
  start(samplePos) = start(samplePos) + 2
  width(samplePos) = 1
  toFilter = paste(seqnames(samplePos), start(samplePos), sep = "_")
  samplePos = samplePos[!toFilter %in% posToFilter]
  TNchr = lapply(as.character(unique(muts$chr)), function(cr){ 
    n = sum(muts$context[muts$chr == cr] %in% pent2search)
    toSamp = which(seqnames(samplePos) == cr)
    samp = sample(toSamp, size=n)
    return(samplePos[samp])
  })
  TNchr = do.call(c, TNchr)
  TNchr = as.data.frame(TNchr)
  TNchr = TNchr[,c("seqnames", "start", "context")]
  colnames(TNchr) = c("chr", "pos", "context")
  TNchr$chr = as.character(TNchr$chr)
  return(TNchr)
}); cat('\n')
print("done with TNs")
stopCluster(cl)
TNs = do.call(rbind,TNs)
#####


############# copy&paste from Corinna's script #######################

# save #####
TNs$ref = substr(TNs$context, 3,3)
TNs$alt = NA
TNs$mutated = 0
TPs$mutated = 1
Muts = rbind(TNs,TPs)
Muts = Muts[order(Muts[,1], Muts[,2]),]
Muts = Muts[,c("chr", "pos", "ref", "alt", "context", "mutated")]
ids = do.call(paste,c(Muts, sep = "_"))
MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                ids)
save(Muts, file = paste0("data/MutTables/scWholeGenomeData/scWGSMuts_",cell2ID[cell],".RData"))
write.table(MutsBed, file = paste0("data/MutTables/scWholeGenomeData/scWGSMuts_",cell2ID[cell],"cell1.bed"),
            col.names = F, row.names = F, sep = "\t", quote = F)
print("end")
#####
# })



