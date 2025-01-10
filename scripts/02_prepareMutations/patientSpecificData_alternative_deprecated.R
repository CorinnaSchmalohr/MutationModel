# preparation #####
library(stringr)
library(GenomicRanges)
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


# load sequences for TN generation
exonSeqs = read.table("data/procData/codExons.withsequences.bed")
# in this file, covered regions are extended by two bases in
# both direction, in order to be able to find the 5mers.
colnames(exonSeqs) = c("chr", "start", "end", "id","sequence")
exonSeqs = exonSeqs[exonSeqs$chr %in% chrs,]
#####


# for each tissue #####
set.seed(235)
for(tissue in names(tissue2Cancer)){
   print(tissue)
   
   # load intermediate files from previous script
   load(paste0("data/rdata/", tissue, "/submuts.RData")) #submuts
   samps = table(submuts$Tumor_Sample_Barcode)
   samps = names(sort(samps, decreasing=T)[1:20])
   
   # get TNs #####
   sapply(samps, function(samp){
      TPs = submuts[submuts$Tumor_Sample_Barcode != samp,]
      posToFilter = paste(TPs$Chromosome, TPs$Start_Position, sep = "_")
      context = TPs$CONTEXT
      contextL = unique(str_length(context))
      pents = substr(context,start=(contextL-1)/2-1, stop=(contextL-1)/2+3)
      TNs = sapply(unique(pents), function(pent){
         cat(pent, ' ')
         nmatch = str_count(string = exonSeqs$sequence, pattern = pent)
         weights = nmatch/sum(nmatch)
         toRm = rep(1,sum(pents == pent)) 
         pos = NULL
         while(length(toRm > 0)){
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
               c(gene["chr"], 
                 pos = (as.numeric(gene["start"]) + samp +2),
                 ref = substr(pent,3,3),
                 context = pent)
            }))
            pos = rbind(pos,new)
            temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
            #  make sure that none of the positions are actually overlapping with existing positions:
            toRm = which(temp %in% posToFilter | duplicated(temp))
            if(length(toRm > 0)){
               pos = pos[-toRm,]
            }
            
         }
      })
      TNs = do.call(rbind,TNs)
      TNs = data.frame(chr=TNs[,"chr"], pos = as.integer(TNs[,"pos"]),
                       ref=TNs[,"ref"], alt = NA, context=TNs[,"context"], 
                       stringsAsFactors = F)
      # load(paste0("data/rdata/", tissue, "/TNs.RData")) #TNs
      
      # combine TPs and TNs
      TPs = cbind(TPs[,c("Chromosome", "Start_Position", 
                         "Reference_Allele","Tumor_Seq_Allele2", 
                         "Tumor_Sample_Barcode")], 
                  pents)
      TPs$mutated = 1
      colnames(TPs) = c("chr", "pos", "ref", "alt", "TumorID", "context", "mutated")
      TNs$mutated = 0
      TNs$TumorID = NA
      exomemuts = exomemuts[,colnames(TNs)]
      
      MutsWIDs = rbind(TNs,exomemuts)
      MutsWIDs = MutsWIDs[order(MutsWIDs[,1], MutsWIDs[,2]),]
      
      # save
      save(MutsWIDs, file = paste0("data/rdata/", tissue, "/MutsWithTumorIDs.RData"))
   }
   #####
   