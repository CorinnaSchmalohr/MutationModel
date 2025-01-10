# preparation #####
library(stringr)
library(GenomicRanges)
tissue2Cancer = list("lung" = "Lung adenocarcinoma",
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
# load sequences for TN generation
exonSeqs = read.table("data/processedData/codExons_noUTR_filtered.withsequences.bed")
# in this file, covered regions are extended by two bases in 
# both direction, in order to be able to find the 5mers.
colnames(exonSeqs) = c("chr", "start", "end", "id","sequence")
exonSeqs = exonSeqs[exonSeqs$chr %in% chrs,]

#####



# for each tissue #####
load("data/MutTables/exomeTrainData/muts.RData")
set.seed(235)
for(tissue in names(tissue2Cancer)){
   print(tissue)
  submuts = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]
  
  # exclude positions that were mutated more than once, 
  # to exclude possible selection
  positions = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
  multiple = duplicated(positions) | duplicated(positions, fromLast=T)
  submuts = submuts[!multiple,]
  
   # loa      d intermediate files from previous script
   load(paste0("data/rdata/", tissue, "/submuts.RData")) #submuts
   posToFilter = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
   context = submuts$CONTEXT
   contextL = unique(str_length(context))
   pents = substr(context,start=(contextL-1)/2-1, stop=(contextL-1)/2+3)
   load(paste0("data/rdata/", tissue, "/TNs.RData")) #TNs

   
   # combine TPs and TNs
   exomemuts = cbind(submuts[,c("Chromosome", "Start_Position", 
                                "Reference_Allele","Tumor_Seq_Allele2", 
                                "Tumor_Sample_Barcode")], 
                     pents)
   exomemuts$mutated = 1
   colnames(exomemuts) = c("chr", "pos", "ref", "alt", "TumorID", "context", "mutated")
   TNs$mutated = 0
   TNs$TumorID = NA
   exomemuts = exomemuts[,colnames(TNs)]
   
   MutsWIDs = rbind(TNs,exomemuts)
   MutsWIDs = MutsWIDs[order(MutsWIDs[,1], MutsWIDs[,2]),]
   
   
   # save
   save(MutsWIDs, file = paste0("data/rdata/", tissue, "/MutsWithTumorIDs.RData"))
}
#####
