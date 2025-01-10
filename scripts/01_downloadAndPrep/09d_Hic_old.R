dir.create("data/predictors/HiC_Encode")
library(rtracklayer)
library(rhdf5)
library(strawr) # install_github("aidenlab/straw/R")
library(Matrix)
library(GenomicRanges)


# for each tissue
tissues = list.files("data/rawdata/HiC_Encode/")
dumpVar = sapply(tissues, function(tissue){
   print(tissue)
   # load meta file
   meta = read.table(paste0("data/rawdata/HiC_Encode/", tissue, "/subMeta.tsv"),
                     sep = "\t", header = T)
   intMeta  = meta[grep("chromatin interactions", meta$Output.type),]
   # for each experiment
   GRs = apply(intMeta, function(experiment){
      ## get list of interaction matrix files
      sub = intMeta[intMeta$Experiment.accession == experiment,]
      print(paste0(experiment, ": ", nrow(sub), " files"))
      ## make GR for each interaction file
      GRs = apply(sub, 1, function(row){
         cat(row["File.accession"], ' ')
         ### load interaction matrix
         format = row["File.format"]
         ### either h5 with rhdf5
         if(format %in% c("h5", "hdf5")){
            # prepare meta info for h5 file
            filePath = paste0("data/rawdata/HiC_Encode/", tissue, "/",
                              row["File.accession"],".h5")
            # load available chromosomes
            chrs = h5read(filePath, "/chrs")
            # get indices of chromosomes
            chrbinRanges = t(h5read(filePath, "/chr_bin_range"))+1
            dimnames(chrbinRanges) = list(chrs, c("first", "last"))
            # get bin borders and balance values
            bins = data.frame(t(h5read(filePath, "/bin_positions")))
            colnames(bins) = c("chr", "start", "end")
            bins$chr = chrs[bins$chr+1] # is 0-indexed.
            # subset to autosomes
            chrs = chrs[grep("[XYM]", chrs, invert=T)]
            
            # iterate through chromosomes and load data
            GRperChr = lapply(chrs, function(chrom){
               subBin = bins[bins$chr == chrom,]
               chrIndex = chrbinRanges[chrom, "first"]:chrbinRanges[chrom, "last"]
               ints = h5read(file = filePath,
                             name = "/interactions", 
                             index = list(chrIndex, chrIndex))
               gr = GRanges(seqnames = subBin$chr, 
                            ranges = IRanges(start=subBin$start, end=subBin$end),
                            score = rowMeans(ints, na.rm = T))
               gr = gr[!is.na(score(gr))]
               h5closeAll()
               return(gr)
            })
            GR = suppressWarnings(do.call(c,GRperChr))
         }
         
         ### or hic file with straw
         else if(format == "hic"){
            GRperChr = lapply(as.character(1:22), function(chrom){
               hic = straw(matrix = "oe", # observed, oe
                           norm = "KR", 
                           fname = paste0("data/rawdata/HiC_Encode/", 
                                          tissue, "/", 
                                          row["File.accession"], ".hic"), 
                           chr1loc=chrom, chr2loc=chrom, 
                           unit="BP", binsize=25000)
               bins = sort(unique(c(hic$x, hic$y)))
               bins2Index = setNames(seq_len(length(bins)), bins)
               mat = sparseMatrix(i = bins2Index[as.character(hic$x)], 
                                  j=bins2Index[as.character(hic$y)], 
                                  x=as.numeric(hic$counts), 
                                  symmetric = T)
               # calculate mean over columns. cannot simply use rowMeans
               # because 0-entries of sparse matrix are used too.
               vals = rowSums(mat, na.rm = T) / rowSums(abs(sign(mat)), na.rm = T)
               gr = GRanges(seqnames=paste0("chr", chrom),
                            ranges=IRanges(start=bins,width=25000),
                            score = vals)
               gr = gr[!is.na(score(gr))]
               return(gr)
            })
            GR = suppressWarnings(do.call(c,GRperChr))
         } else{
            print("unknown format")
         }
         return(GR)
      })
      ## save GR-List
      GRL = GRangesList(GRs)
      names(GRL) = sub$File.accession
      save(GRs, file = paste0("data/predictors/HiC_Encode/", tissue,
                              "_ints_", experiment,".RData"))
   })
   ## get compartment files
   compMeta = meta[meta$Output.type == "genome compartments",]
   ## for each compartment file
   dumpVar2 = sapply(unique(intMeta$Experiment.accession), function(experiment){
      
   })
   temp = import.bw(con = paste0("data/rawdata/HiC_Encode/", tissue, "/",
                                 x, ".bigWig"))
   start(temp) = start(temp)-1
   sd(score(temp), na.rm = T)
   ### load bw with rtracklayer
   ### scale
   ### generate and return GR
   ## combine to GR-List
   ## save GRL (for each experiment, name: organ_tissue/cellline_comp.RData)
})
})



# read in bigWigs to compare comparment bigWigs between different tissues
testGR = GRanges(seqnames=rep("chr1", 1000000), 
                 ranges=IRanges(start = seq(1,100000000, by = 100), width=1))
comp = sapply(meta$File.accession, function(x){
   cat(x, ' ')
   temp = import.bw(con = paste0("data/rawdata/HiC_Encode/", tissue, "/",
                                 x, ".bigWig"))
   start(temp) = start(temp)-1
   sd(score(temp), na.rm = T)
   
   ov = as.data.frame(findOverlaps(query=testGR, subject=temp))
   map = split(score(temp)[ov$subjectHits], ov$queryHits)
   res = sapply(map, mean, na.rm = T)
   res = res[as.character(1:length(testGR))]
   return(res)
   # return(temp)
   # return(temp[seqnames(temp) == "chr1"])
})
comp2 = scale(comp)
meanSignal = rowMeans(comp2, na.rm=T)
plot(start(testGR),meanSignal, type = "l", 
     ylim = c(min(meanSignal, na.rm = T)*2, max(meanSignal, na.rm = T)*2),
     xlab = "genomic position chr1", ylab = "HiC compartment PC1")
sapply(1:ncol(comp2) ,function(i){
   lines(start(testGR),comp2[,i], col = rainbow(ncol(comp))[i])
})
lines(start(testGR),meanSignal, lwd = 2)
abline(h = 0)
legend("bottomright", col = c(rainbow(ncol(comp)), "black"), 
       legend=c(colnames(comp), "mean"),
       lwd = c(rep(1,ncol(comp)),2), cex = 0.5, ncol = 2)


abline(h = 0, lty = 2)
multi = read.table("multi.bg")
multi = multi[multi$V1 == "chr1",]
multi2 = cbind(multi[,1:3],scale(multi[,-(1:3)]))
plot(multi2$V2,rowMeans(multi2[,-(1:3)], na.rm=T), type = "l", xlim = c(1,100000000))
sapply(4:ncol(multi2),function(i){
   lines(multi2$V2, multi2[,i], col = i)
})
#####
