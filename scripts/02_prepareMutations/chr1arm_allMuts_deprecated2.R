
# load exonic regions
#load("data/procData/codExons.RData")
load("data/procData/codExons_noUTR.RData")


# get all positions on the long arm of chr1
chr1exons = codExons_noUTR[codExons_noUTR$chr == "chr1",]
sum(chr1exons$end - chr1exons$start) 

# chr1armExons = chr1exons[chr1exons$end < 130000000,]
# sum(chr1armExons$end - chr1armExons$start) # just the arm is feasible

pos = do.call(c,apply(chr1exons,1, function(x){
   x["start"]:x["end"]
}))
length(pos)
pos = unique(pos)
length(pos)
Muts = data.frame(chr = "chr1", pos = pos, ref = NA, alt = NA)

# get context 
bed = data.frame(Muts$chr, as.integer(Muts$pos-3), as.integer(Muts$pos+2))
write.table(bed, file="temp/TPs_bed.bed", quote=F, col.names=F, row.names=F, sep="\t")
cmd = paste0("bedtools getfasta -fi data/rawdata/GRCh37.p11.genome.fa ", 
             "-bed temp/TPs_bed.bed ",
             "-name -tab | cut -f 2")
context = system(cmd, intern = T)
Muts$context = context
Muts$ref = substr(Muts$context, 3,3)
Muts$mutated = NA

# save 
ids = do.call(paste,c(Muts, sep = "_"))
MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                ids)
save(Muts, file = "data/rdata/Chr1_allPos.RData")
write.table(MutsBed, file = "data/procData/Chr1_allPos.bed",
            col.names = F, row.names = F, sep = "\t", quote = F)


# save in 4 parts because of memory problems otherwise
dir.create("data/rdata/Chr1arm/", showWarnings = F)
dir.create("data/procData/Chr1arm/", showWarnings = F)
n = nrow(Muts)
q = ceiling(n/4)
seqs = list(1:q, (q+1):(q*2), (q*2+1):(q*3),(q*3+1):n)
for(i in 1:4){
   subMuts = Muts[seqs[[i]],]
   subMutsBed = MutsBed[seqs[[i]],]
   save(subMuts, file = paste0("data/rdata/Chr1arm/Chr1arm_allPos_part",i,".RData"))
   write.table(subMutsBed, file = paste0("data/procData/Chr1arm/Chr1arm_allPos_part",i,".bed"),
               col.names = F, row.names = F, sep = "\t", quote = F)
}
