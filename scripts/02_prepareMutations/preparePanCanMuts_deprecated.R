
# preparation #####
tissue2Cancer = list("luad" = "Lung adenocarcinoma",
                     "breast" = "Breast invasive carcinoma",
                     "skin" = "Skin Cutaneous Melanoma",
                     "colon" = c("Colon adenocarcinoma", "Rectum adenocarcinoma"),
                     "ovary" = "Ovarian serous cystadenocarcinoma",
                     "kidney" = "Kidney renal clear cell carcinoma",
                     "prostate" = "Prostate adenocarcinoma",
                     "esophagus" = "Esophageal carcinoma",
                     "liver" = "Liver hepatocellular carcinoma",
                     "brain" = "Brain Lower Grade Glioma")
chrs = paste0("chr", c(1:22))
library(stringr)
library(GenomicRanges)
#####

# load the pancan mutations #####
muts = read.table("data/rawdata/pancan/mc3.v0.2.8.PUBLIC.maf.gz",
                  header = T, na.strings=".", quote="", as.is=T)
# get  cancer IDs
cancerIDs = read.table("data/rawdata/pancan/TCGA-CDR-SupplementalTableS1_TSScode_info.csv",
                       sep = ",", header = T, as.is = T, skip=2, na.strings="")
cancerID2cancer = setNames(cancerIDs$Disease.type, cancerIDs$TSS.Code)
cancerCode = substr(muts$Tumor_Sample_Barcode, 6,7)
muts$cancerType = cancerID2cancer[cancerCode]
muts = muts[muts$cancerType %in% unlist(tissue2Cancer) & muts$Variant_Type == "SNP",]
muts$Chromosome = paste0("chr", muts$Chromosome)
muts = muts[muts$Chromosome %in% chrs,]
#####

# filter out variants with too little depth #####
muts = muts[muts$n_depth>=8 & muts$t_depth >= 14,]
# these cutoffs are based on what was used for the "sufficient coverage" files.
#####


# filter out variants with >1% population frequency  #####
mafs = muts[,c("GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF",
               "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF")]
mafsFilter = apply(mafs,1, function(x){
   y = gsub("[T,C,G,A,-]:", replacement="",x=x)
   z = as.numeric(unlist(strsplit(split=",", x=y)))
   if(all(is.na(z))){
      return(T)
   }  else if(max(z, na.rm = T) < 0.01){
      return(T)
   } else{
      return(F)
   }
})
muts = muts[mafsFilter,]
#####


# subset to exome regions #####
load("data/procData/codExons_noUTR_GR.RData")
#load("data/procData/codExons.RData")
#exonGR = GRanges(seqnames=codExons$chr, 
#                 ranges = IRanges(start = codExons$start, end = codExons$end))
mutGR = GRanges(seqnames= muts$Chromosome, 
                ranges=IRanges(start = muts$Start_Position,
                               end=muts$End_Position))
ov = !is.na(findOverlaps(query=mutGR, subject=codExons_noUTR_GR, select = "first"))
print(table(ov))
muts = muts[ov,]
#####


# load sequences for TN generation #####
exonSeqs = read.table("data/procData/codExons_noUTR.withsequences.bed")
# in this file, covered regions are extended by two bases in 
# both direction, in order to be able to find the 5mers.
colnames(exonSeqs) = c("chr", "start", "end", "id","sequence")
exonSeqs = exonSeqs[exonSeqs$chr %in% chrs,]
#####



# for each tissue #####
set.seed(235)
for(tissue in names(tissue2Cancer)){
   print(tissue)
   dir.create(paste0("data/rdata/", tissue, "/"), showWarnings = F)
   # subset to tissue
   submuts = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]

   # exclude positions that were mutated more than once, 
   # to exclude possible selection
   positions = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
   multiple = duplicated(positions) | duplicated(positions, fromLast=T)
   submuts = submuts[!multiple,]
   
   # save intermediate 
   save(submuts, file = paste0("data/rdata/", tissue, "/submuts.RData"))
   
   # get TNs #####
   posToFilter = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
   context = submuts$CONTEXT
   contextL = unique(str_length(context))
   # wrong:
   # pents = substr(context,start=(contextL-1)/2-2, stop=(contextL-1)/2+2)
   pents = substr(context,start=(contextL-1)/2-1, stop=(contextL-1)/2+3)
   if(!all.equal(substring(pents, 3,3), submuts$Reference_Allele)){
      error("context is not matching reference base. check script")
   }
   TNs = sapply(unique(pents), function(pent){
      # cat(pent, ' ')
      nmatch = str_count(string = exonSeqs$sequence, pattern = pent)
      weights = nmatch/sum(nmatch)
      TNs_sample = sample(1:nrow(exonSeqs),
                          size = sum(pents == pent),
                          prob = weights, replace = T)
      pos = t(sapply(TNs_sample, function(i){
         gene = unlist(exonSeqs[i,])
         matches = gregexpr(pattern = pent, text = gene["sequence"])[[1]]
         if(length(matches) == 1){
            samp = matches[1]
         } else{
            samp = sample(matches,size = 1)
         }
         c(ref = substr(pent,3,3),
           gene["chr"], 
           pos = (as.numeric(gene["start"]) + samp +2),
           context = pent)
      }))
      # make sure that none of the positions are actually overlapping with existing positions
      temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
      toRm = which(temp %in% posToFilter | duplicated(temp))
      while (length(toRm) > 0) {
         TNs_sample = sample(1:nrow(exonSeqs),
                             size = length(toRm),
                             prob = weights, replace = T)
         new = t(sapply(TNs_sample, function(i){
            gene = unlist(exonSeqs[i,])
            matches = gregexpr(pattern = pent, text = gene["sequence"])[[1]]
            if(length(matches) == 1){
               samp = matches[1]
            } else{
               samp = sample(matches,size = 1)
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
   TNs = do.call(rbind,TNs)
   TNs = data.frame(chr=TNs[,"chr"], pos = as.integer(TNs[,"pos"]),
                    ref=TNs[,"ref"], alt = NA, context=TNs[,"context"], 
                    stringsAsFactors = F)
   
   save(TNs,file = paste0("data/rdata/", tissue, "/TNs.RData"))
   # write.table(TNs, file=paste0("data/procData/", tissue, "/TNs.tsv"),
   #             col.names = F, row.names = F, sep = "\t", quote = F)
   
   
   # combine TPs and TNs
   exomemuts = cbind(submuts[,c("Chromosome", "Start_Position", 
                                "Reference_Allele","Tumor_Seq_Allele2")], 
                     pents)
   colnames(exomemuts) = colnames(TNs)
   TNs$mutated = 0
   exomemuts$mutated = 1
   Muts = rbind(TNs,exomemuts)
   Muts = Muts[order(Muts[,1], Muts[,2]),]
   
   
   # save
   ids = do.call(paste,c(Muts, sep = "_"))
   MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                   ids)
   save(Muts, file = paste0("data/rdata/", tissue, "/Muts.RData"))
   dir.create(paste0("data/procData/", tissue),showWarnings=F)
   write.table(MutsBed, file = paste0("data/procData/", tissue, "/Muts_noUTR.bed"),
               col.names = F, row.names = F, sep = "\t", quote = F)
   # MutsBedwithIDs = cbind(MutsBed, 1:nrow(MutsBed))
   # write.table(MutsBedwithIDs, file = paste0("data/procData/", tissue, "/MutsWithIDs.bed"),
   #             col.names = F, row.names = F, sep = "\t", quote = F)
}
#####

# get number of mutations per sample per tissue #####
nPerTissue = sapply(names(tissue2Cancer), function(tissue){
   print(tissue)
   # subset to tissue
   submuts = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]
   table(submuts$Tumor_Sample_Barcode)
})
nPerTissue = nPerTissue[sort(names(nPerTissue))]
save(nPerTissue, file="data/rdata/nPerTissue.RData")
####