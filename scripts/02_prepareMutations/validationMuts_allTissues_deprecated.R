### Packages ######
library(readxl)
library(tibble)
library(GenomicRanges)
library(data.table)
library(stringr)
library(dplyr)
library(Biostrings)

###############


# read in reference gtf #####
#load("data/procData/codExons.RData")
load("data/processedData/codExons_noUTR_filtered.RData")
gr_exon = GRanges(seqnames=codExons_filtered$chr, 
                  ranges = IRanges(start = codExons_filtered$start, end = codExons_filtered$end))
# and the respective sequences. careful, ranges were extended by 2 on either
# end to get complete context
exonSeqs = read.table("data/processedData/codExons_noUTR_filtered.withsequences.bed")
colnames(exonSeqs) = c("chr", "start", "stop", "name", "sequence")

exonSeqs = exonSeqs[!exonSeqs$chr %in% c("chrM", "chrY", "chrX"),]
######



# function to process a table of mutations to filter and create TN positions #####
load("data/processedData/pentamers.RData")

processMutations = function(dat, nsamples = NA, grExon = gr_exon, pent2context=pent2context){
  # get table of TPs
  # remove mutations that occured more than once.
  posToFilterTotal = cbind(dat$chr,
                           paste(dat$chr, dat$pos, sep = "_"))
  if(!is.na(nsamples)){
    dat = dat[nsamples<2,]
  }
  dupl = duplicated(paste(dat$chr, dat$pos, sep = "_")) | 
    duplicated(paste(dat$chr, dat$pos, sep = "_"), fromLast = T)
  dat = dat[!dupl,]
  # only chrs 1-22
  dat = dat[dat$chr %in% paste0("chr", 1:22),]
  # create a GR object
  gr_data <- makeGRangesFromDataFrame(dat, 
                                      start.field = "pos", 
                                      end.field = "pos", 
                                      seqnames.field = "chr")
  # Find overlap between data and exon granger object
  ov = !is.na(findOverlaps(query=gr_data, subject=grExon, select = "first"))
  print(table(ov)) 
  TPs = dat[ov,]
  # get sequence context for positions
  TPs_bed = data.frame(TPs$chr, TPs$pos-3, TPs$pos+2)
  write.table(TPs_bed, file="temp/TPs_bed.bed", quote=F, col.names=F, row.names=F, sep="\t")
  cmd = paste0("/data/public/cschmalo/anaconda3/envs/MutModel/bin/bedtools getfasta -fi data/rawdata/GRCh37.primary_assembly.genome.fa ", 
               "-bed temp/TPs_bed.bed ",
               "-name -tab | cut -f 2")
  
  context = system(cmd, intern = T)
  TPs$context = context
  if(!all.equal(substring(context, 3,3), TPs$ref)){
    error("context is not matching reference base. check script")
  }
  print("getting TNs")
  TNs = lapply(unique(TPs$chr), function(chr){ 
    cat(chr, ' ')
    #iterate through chromosomes so that equal number of TNs and TPs per chr
    # subset data to chromosome
    TPs_chr = TPs[TPs$chr == chr,]
    exonSeqs_chr = exonSeqs[exonSeqs$chr == chr,]
    pents = TPs_chr$context
    posToFilter = posToFilterTotal[posToFilterTotal[,1] == chr,]
    
    TNs_chr = lapply(unique(pent2context[pents]), function(pent){
      # cat(pent, ' ')
      pent2search = fivemers[pent,]
      nmatch1 = str_count(string = exonSeqs_chr$sequence, pattern = pent2search[1])
      nmatch2 = str_count(string = exonSeqs_chr$sequence, pattern = pent2search[2])
      nmatch = nmatch1 + nmatch2
      weights = nmatch/sum(nmatch)
      toRm = rep(1,sum(pents %in% pent2search)) 
      pos = NULL
      while(length(toRm) > 0){
        TNs_sample = sample(1:nrow(exonSeqs_chr),
                            size = length(toRm),
                            prob = weights, replace = T)
        new = t(sapply(TNs_sample, function(i){
          # cat(i,' ')
          gene = unlist(exonSeqs_chr[i,])
          matches1 = cbind(gregexpr(pattern = pent2search[1], text = gene["sequence"])[[1]],pent2search[1])
          matches2 = cbind(gregexpr(pattern = pent2search[2], text = gene["sequence"])[[1]],pent2search[2])
          matches = rbind(matches1, matches2)
          matches = matches[matches[,1] != "-1",, drop = F]
          if(nrow(matches) == 1){
            sampMatch = matches[1,]
          } else{
            sampMatch = matches[sample(1:nrow(matches),size = 1),]
          }
          samp = as.integer(sampMatch[1])
          c(gene["chr"], 
            pos = (as.numeric(gene["start"]) + samp +2),
            ref = unname(substr(gene["sequence"],start = samp+2, stop = samp+2)),
            context = sampMatch[2])
        }))
        pos = rbind(pos,new)
        temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
        #  make sure that none of the positions are actually overlapping with existing positions:
        toRm = which(temp %in% posToFilter | duplicated(temp))
        if(length(toRm) > 0){
          pos = pos[-toRm,]
        }
      }
      return(pos)
    })
    TNs_chr = do.call(rbind,TNs_chr)
    
    #old
    # TNs_chr = lapply(unique(pents), function(pent){
    #   # cat(pent, ' ')
    #   nmatch = str_count(string = exonSeqs_chr$sequence, pattern = pent)
    #   weights = nmatch/sum(nmatch)
    #   toRm = rep(1,sum(pents == pent)) 
    #   pos = NULL
    #   while(length(toRm > 0)){
    #     TNs_sample = sample(1:nrow(exonSeqs_chr),
    #                         size = length(toRm),
    #                         prob = weights, replace = T)
    #     new = t(sapply(TNs_sample, function(i){
    #       gene = unlist(exonSeqs_chr[i,])
    #       matches = gregexpr(pattern = pent, text = gene["sequence"])[[1]]
    #       if(length(matches) == 1){
    #         samp = matches[1]
    #       } else{
    #         samp = sample(matches,size = 1)
    #       }
    #       c(gene["chr"], 
    #         pos = (as.numeric(gene["start"]) + samp +2),
    #         ref = substr(pent,3,3),
    #         context = pent)
    #     }))
    #     pos = rbind(pos,new)
    #     temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
    #     #  make sure that none of the positions are actually overlapping with existing positions:
    #     toRm = which(temp %in% posToFilter | duplicated(temp))
    #     if(length(toRm > 0)){
    #       pos = pos[-toRm,]
    #     }
    #     
    #   }
    #   return(pos)
    # })
    # TNs_chr = do.call(rbind,TNs_chr)
  })
  cat('\n')
  TNs = do.call(rbind,TNs)
  TNs = data.frame(chr=TNs[,"chr"], pos = as.integer(TNs[,"pos"]),
                   ref=TNs[,"ref"],  alt = NA, context=TNs[,"context"], 
                   stringsAsFactors = F)
  # Combine TNs and TPs
  TPs$mutated <- 1
  TNs$mutated <- 0
  Muts <- rbind(TPs, TNs)
  Muts = Muts[order(Muts[,1], Muts[,2]),]
  return(Muts)
}
#####



# SomaMutDB #####
tissues = c("brain","breast", "colon","esophagus", "kidney", "liver", "lung","prostate", "skin")
colname = c("chr","pos","ref","alt","tissue","tissue_detail","sex","age","gene_id","gbiotype","symbol","strand",
            "feature","distance_upstream_downstream","sift","polyphen","cadd_phred","regulatory_id","regulatory_type")
nsamples = NA
dir.create("data/MutTables/healthyTissues/", showWarnings = F)
sapply(tissues, function(tissue){
  print(tissue)
  rawdata = read.table(paste0("data/rawdata/validation/SomaMutDB/",tissue,".tsv"), sep = "\t")
  colnames(rawdata) = colname
  
  # Remove indels 
  rawdata <- rawdata[nchar(rawdata$ref) == 1 & nchar(rawdata$alt) == 1,]
  
  data_adj <- data.frame(chr = paste0("chr", rawdata$chr), 
                         pos = rawdata$pos,
                         ref = rawdata$ref,
                         alt = rawdata$alt)
  Muts = processMutations(dat=data_adj, nsamples = nsamples)
  save(Muts, file = paste0("data/MutTables/healthyTissues/", tissue, "_SomaMutDB.RData"))
  
  ids = do.call(paste,c(Muts, sep = "_"))
  MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                  ids)
  write.table(MutsBed, file = paste0("data/MutTables/healthyTissues/", tissue, "_SomaMutDB.bed"),
              col.names = F, row.names = F, sep = "\t", quote = F)
})


######

