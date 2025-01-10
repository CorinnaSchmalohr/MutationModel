# this file extracts all QTL from GTEx for each tissue and saves 
# them as a GenomicRanges object with nominal p-value and slope as 
# genomic feature scores.
library(GenomicRanges)
library(rtracklayer)
if(!dir.exists("data/predictors/GTEx_eQTL")){
   dir.create("data/predictors/GTEx_eQTL")
}
files = list.files("data/rawdata/GTEx_Analysis_v8_eQTL/", 
           pattern="signif_variant_gene_pairs.txt.gz")
tissues = sapply(files, function(x){
   strsplit(x, "[.]")[[1]][1]
   }, USE.NAMES=F)
# v8 GTEx is in hg38 --> liftover to hg19
ch = import.chain("data/rawdata/genome/hg38ToHg19.over.chain")
eqtl = lapply(tissues, function(tissue){
   print(tissue)
   eqtl_tab = read.table(paste0("data/rawdata/GTEx_Analysis_v8_eQTL/",
                                tissue,
                                ".v8.signif_variant_gene_pairs.txt.gz"),
                         sep = "\t", as.is = T, header = T)
   pos = do.call(rbind,strsplit(eqtl_tab[,1], split="_"))
   varl = nchar(pos[,3]) # I decided against subsetting to snps only, contains indels as well
   gr = GRanges(seqnames=pos[,1], 
                ranges=IRanges(start=as.integer(pos[,2]), 
                               end=as.integer(pos[,2])+varl-1),
                val = -log(eqtl_tab$pval_nominal))
   genome(gr) = "hg38"
   gr = liftOver(gr, ch)
   gr = unlist(gr)
   genome(gr) = "hg19"
   gr = sort(gr)
   
   # combine multiple eQTLs for the same SNPs by taking mean
   newGR <- reduce(gr, ignore.strand = TRUE, with.revmap = TRUE)
   vals = gr$val
   toCombine = as.list(mcols(newGR)$revmap)
   newVals = sapply(toCombine, function(x){
     mean(vals[x])
   })
   newGR$val = newVals
   gr = newGR
   gr$revmap = NULL
   rm(newGR, vals, toCombine, newVals)
   # save
   gr = sort(gr)
   save(gr, file=paste0("data/predictors/GTEx_eQTL/", tissue, ".RData"))
   gr = GRanges(seqnames=pos[,1], 
                ranges=IRanges(start=as.integer(pos[,2]), 
                               end=as.integer(pos[,2])+varl-1),
                val = abs(eqtl_tab$slope))
   genome(gr) = "hg38"
   gr = liftOver(gr, ch)
   gr = unlist(gr)
   genome(gr) = "hg19"
   gr = sort(gr)
   newGR <- reduce(gr, ignore.strand = TRUE, with.revmap = TRUE)
   vals = gr$val
   toCombine = as.list(mcols(newGR)$revmap)
   newVals = sapply(toCombine, function(x){
     mean(vals[x])
   })
   newGR$val = newVals
   gr = newGR
   rm(newGR, vals, toCombine, newVals)
   gr$revmap = NULL
   gr = sort(gr)
   save(gr, file=paste0("data/predictors/GTEx_eQTL/", tissue, "_slope.RData"))
})
# create a gr combining all brain p-vals
brainPvals = lapply(tissues[grep(tissues, pattern="Brain")], function(tissue){
   pv = get(load(paste0("data/predictors/GTEx_eQTL/", tissue, ".RData")))
   return(pv)
})
gr = do.call(c,brainPvals)
newGR <- reduce(gr, ignore.strand = TRUE, with.revmap = TRUE)
vals = gr$val
toCombine = as.list(mcols(newGR)$revmap)
newVals = sapply(toCombine, function(x){
  mean(vals[x])
})
newGR$val = newVals
gr = newGR
gr$revmap = NULL
rm(newGR, vals, toCombine, newVals)
# save
gr = sort(gr)
save(gr, file=paste0("data/predictors/GTEx_eQTL/brain_combined.RData"))
# and the same for the slopes
brainSlopes = lapply(tissues[grep(tissues, pattern="Brain")], function(tissue){
   slope = get(load(paste0("data/predictors/GTEx_eQTL/", tissue, "_slope.RData")))
   return(slope)
})
gr = do.call(c,brainSlopes)
newGR <- reduce(gr, ignore.strand = TRUE, with.revmap = TRUE)
vals = gr$val
toCombine = as.list(mcols(newGR)$revmap)
newVals = sapply(toCombine, function(x){
  mean(vals[x])
})
newGR$val = newVals
gr = newGR
gr$revmap = NULL
rm(newGR, vals, toCombine, newVals)
save(gr, file=paste0("data/predictors/GTEx_eQTL/brain_slopes_combined.RData"))




