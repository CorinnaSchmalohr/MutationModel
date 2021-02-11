# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


# get TNs #####
library(stringr)
library(parallel)
load(paste0("data/rdata/", tissue, "/submuts.RData"))
load(paste0("data/rdata/", tissue, "/exomemuts.RData"))
chrs = paste0("chr", c(1:22))
covered_forTNs = read.table(paste0("data/procData/", tissue,
                             "/covered_regions_withsequences.bed"),as.is=T)
# in this file, covered regions are extended by two bases in 
# both direction, in order to be able to find the 5mers.
colnames(covered_forTNs) = c("chr", "start", "end", "sequence")
covered_forTNs = covered_forTNs[covered_forTNs$chr %in% chrs,]

posToFilter = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
set.seed(235)
print("starting")
pents = unique(exomemuts$CONTEXT)
TNs = sapply(pents, function(pent){
   nmatch = str_count(string = covered_forTNs$sequence, pattern = pent)
   weights = nmatch/sum(nmatch)
   TNs_sample = sample(1:nrow(covered_forTNs),
                       size = sum(exomemuts$CONTEXT == pent),
                       prob = weights, replace = T)
   # cl = makeCluster(10,type = "FORK")
   # pos = t(parSapply(cl,TNs_sample, function(i){
   pos = t(sapply(TNs_sample, function(i){
      gene = unlist(covered_forTNs[i,])
      As = gregexpr(pattern = pent, text = gene["sequence"])[[1]]
      if(length(As) == 1){
         samp = As[1]
      } else{
         samp = sample(As,size = 1)
      }
      c(ref = substr(pent,3,3),
        gene["chr"], 
        pos = (as.numeric(gene["start"]) + samp +2),
        context = pent)
   }))
   # stopCluster(cl)
   # make sure that none of the positions are actually
   # overlapping positions with SNPs
   temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
   toRm = which(temp %in% posToFilter | duplicated(temp))
   while (length(toRm) > 0) {
      TNs_sample = sample(1:nrow(covered_forTNs),
                          size = length(toRm),
                          prob = weights, replace = T)
      new = t(sapply(TNs_sample, function(i){
         gene = unlist(covered_forTNs[i,])
         As = gregexpr(pattern = pent, text = gene["sequence"])[[1]]
         if(length(As) == 1){
            samp = As[1]
         } else{
            samp = sample(As,size = 1)
         }
         c(ref = substr(pent,3,3),
           gene["chr"], 
           pos = (as.numeric(gene["start"]) + samp +2),
           context = pent)
      }))
      pos[toRm,] = new
      temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
      toRm = which(temp %in% posToFilter | duplicated(temp))
   }
   return(pos)
})
print("done")
TNs = do.call(rbind,TNs)
TNs = data.frame(chr=TNs[,"chr"], pos = TNs[,"pos"],
                 ref=TNs[,"ref"], alt = NA, context=TNs[,"context"], stringsAsFactors = F)
exomemuts$End_Position = NULL
colnames(exomemuts) = colnames(TNs)
save(TNs,file = paste0("data/rdata/", tissue, "/TNs.RData"))
write.table(TNs, file=paste0("data/procData/", tissue, "/TNs.tsv"),
            col.names = F, row.names = F, sep = "\t", quote = F)
#####

# combine TPs and TNs #####
TNs$mutated = 0
exomemuts$mutated = 1
Muts = rbind(TNs,exomemuts)
Muts$pos = as.numeric(Muts$pos)
Muts = Muts[order(Muts[,1], Muts[,2]),]
#####


# assign gene id and strand to each variant #####
print("assigning genes")
load("data/rdata/gtf_genes.RData")
library(parallel)
cl = makeCluster(8, type = "FORK")
ingene = parApply(cl,Muts,1,function(m){
   target = gtf_genes[gtf_genes$chr == m[1] & 
                         gtf_genes$start <= as.numeric(m[2]) &
                         gtf_genes$end >= as.numeric(m[2]),]
   g = unique(target[,"gene_id"])
   if (length(g) > 1) {
      lvls = strsplit(target$info,";")
      lvls = as.numeric(sapply(lvls, function(x) {
         strsplit(x[grep("level",x)], " ")[[1]][3] }))
      g = unique(target[lvls == min(lvls) &
                           target$transcript_type == "protein_coding","gene_id"])
   }
   if (length(g) == 0) {
      return(NA)
   }
   return(g)
})
stopCluster(cl)
save(ingene,file =paste0("data/rdata/", tissue, "/ingene.RData"))
ingene[sapply(ingene, length) != 1] = NA
ingene = unlist(ingene)
name2strand = gtf_genes[,"strand"]
names(name2strand) = gtf_genes[,"gene_id"]
strand = name2strand[ingene]
Muts$geneID = ingene
Muts$strand = strand
Muts$baseGeneID = do.call(rbind,strsplit(Muts$geneID,split = ".", fixed = T))[,1]
save(Muts,file = paste0("data/rdata/", tissue, "/Muts.RData"))
#####


# write Bed ####
MutsBed = Muts[,c(1,2,2)]
MutsBed[,2] = as.integer(MutsBed[,2] - 1)
MutsBed[,3] = as.integer(MutsBed[,3])
write.table(MutsBed, file = paste0("data/procData/", tissue, "/Muts.bed"),
            col.names = F, row.names = F, sep = "\t", quote = F)
MutsBedwithIDs = cbind(MutsBed, 1:nrow(MutsBed))
write.table(MutsBedwithIDs, file = paste0("data/procData/", tissue, "/MutsWithIDs.bed"),
            col.names = F, row.names = F, sep = "\t", quote = F)
#####