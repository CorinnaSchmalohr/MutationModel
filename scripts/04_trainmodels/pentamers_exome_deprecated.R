
# load exonic regions
load("data/procData/codExons.RData")


# get all positions in the exome
Muts = lapply(paste0("chr", 1:22), function(cr){
   print(cr)
   chrexons = codExons[codExons$chr == cr,]
   pos = do.call(c,apply(chrexons,1, function(x){
      x["start"]:x["end"]
   }))
   pos = unique(pos)
   m = data.frame(chr = cr, pos = pos, ref = NA, alt = NA)
   return(m)
})



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
n = nrow(Muts)
q = ceiling(n/4)
seqs = list(1:q, (q+1):(q*2), (q*2+1):(q*3),(q*3+1):n)
for(i in 1:4){
   subMuts = Muts[seqs[[i]],]
   subMutsBed = MutsBed[seqs[[i]],]
   save(subMuts, file = paste0("data/rdata/Chr1arm_allPos_part",i,".RData"))
   write.table(subMutsBed, file = paste0("data/procData/Chr1arm_allPos_part",i,".bed"),
               col.names = F, row.names = F, sep = "\t", quote = F)
}
