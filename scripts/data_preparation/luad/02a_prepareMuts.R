setwd("/cellnet/MutationModel")


# read in mutation information for LUAD #####
muts = read.table("data/rawdata/Alexandrov2013/somatic_mutation_data/Lung Adeno/Lung Adeno_clean_somatic_mutations_for_signature_analysis.txt",
                  as.is = T)
colnames(muts) = c("sample","type","chr", "start", "end", "ref", "alt", "comment")
# restrict to mutations that were created from WES
samples = read.table("data/rawdata/Alexandrov2013/samples_summary.txt", 
                     header = T, sep = "\t", as.is = T)
samples = samples[samples$Cancer.Type == "Lung Adeno",]
rownames(samples) = samples$Sample.Name
seqtype = samples[muts$sample,"Sequencing.Type"]
exomemuts = muts[muts$type == "subs" & seqtype == "Exome",]
table(table(exomemuts$sample) > 4000) # there is one sample with 
# exceptionally many mutations. Remove it?
# remove mitochondria and y chromosome
exomemuts = exomemuts[!exomemuts$chr  %in% c("MT","Y", "X"),]
exomemuts$chr = paste0("chr", exomemuts$chr)
colnames(exomemuts)[4] = "pos"
save(exomemuts, file = "data/rdata/luad/exomemuts.RData")
#####


# for each mutation, get a corresponding TN mutation (with same ref base) #####
# load("data/rdata/luad/exomemuts.RData")
load("data/rdata/ref_coding.RData")
library(stringr)
library(parallel)
unfiltered = read.table("data/rawdata/Alexandrov2013/somatic_mutation_data/Lung Adeno/Lung Adeno_raw_mutations_data.txt",
                               as.is = T)
posToFilter = paste("chr", unfiltered$V3, "_", unfiltered$V4, sep = "")
set.seed(235)
TNs = sapply(c("A", "C", "G", "T"), function(base){
   nmatch = str_count(string = ref_coding$sequence, pattern = base)
   weights = nmatch/sum(nmatch)
   TNs_A_sample = sample(1:nrow(ref_coding),
                         size = sum(exomemuts$ref == base),
                         prob = weights, replace = T)
   cl = makeCluster(8,type = "FORK")
   pos = t(parSapply(cl,TNs_A_sample, function(i){
      gene = unlist(ref_coding[i,])
      As = gregexpr(pattern = base, text = gene["sequence"])[[1]]
      c(ref = base,
        gene["chr"], 
        pos = (as.numeric(gene["start"]) + sample(As,size = 1)))
   }))
   stopCluster(cl)
   # make sure that none of the positions are actually
   #  overlapping positions with SNPs
   temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
   toRm = which(temp %in% posToFilter | duplicated(temp))
   while (length(toRm) > 0) {
      TNs_A_sample = sample(1:nrow(ref_coding),
                            size = length(toRm),
                            prob = weights, replace = T)
      new = t(sapply(TNs_A_sample, function(i){
         gene = unlist(ref_coding[i,])
         As = gregexpr(pattern = base, text = gene["sequence"])[[1]]
         c(ref = base,gene["chr"], 
           pos = (as.numeric(gene["start"]) + sample(As,size = 1)))
      }))
      pos[toRm,] = new
      temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
      toRm = which(temp %in% posToFilter | duplicated(temp))
   }
   return(pos)
})
TNs = do.call(rbind,TNs)
TNs = cbind(sample = "TN", TNs, alt = NA)
TNs = as.data.frame(TNs, stringsAsFactors = F)
save(TNs,file = "data/rdata/luad/TNs.RData")
#####


# combine TN and TP table and sort #####
load("data/rdata/luad/exomemuts.RData")
load("data/rdata/luad/TNs.RData")
TNs$mutated = 0
exomemuts$mutated = 1
Muts = rbind(TNs[,c(1,3,4,2,5,6)],exomemuts[,c(1,3,4,6,7,9)])
Muts$pos = as.numeric(Muts$pos)
Muts = Muts[order(Muts[,2], Muts[,3]),]
##### 



# assign gene id and strand to each variant #####
load("data/rdata/gtf_genes.RData")
library(parallel)
cl = makeCluster(8, type = "FORK")
ingene = parApply(cl,Muts,1,function(m){
   target = gtf_genes[gtf_genes$chr == m[2] & 
                         gtf_genes$start <= as.numeric(m[3]) &
                         gtf_genes$end >= as.numeric(m[3]),]
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
save(ingene,file = "data/rdata/luad/ingene.RData")
ingene[sapply(ingene, length) != 1] = NA
ingene = unlist(ingene)
name2strand = gtf_genes[,"strand"]
names(name2strand) = gtf_genes[,"gene_id"]
strand = name2strand[ingene]
Muts$geneID = ingene
Muts$strand = strand
Muts$baseGeneID = do.call(rbind,strsplit(Muts$geneID,split = ".", fixed = T))[,1]
save(Muts,file = "data/rdata/luad/Muts.RData")
#####


# write Bed ####
MutsBed = Muts[,c(2,3,3)]
MutsBed[,2] = as.integer(MutsBed[,2] - 1)
MutsBed[,3] = as.integer(MutsBed[,3])
write.table(MutsBed, file = "data/procData/luad/Muts.bed",
            col.names = F, row.names = F, sep = "\t", quote = F)
#####


