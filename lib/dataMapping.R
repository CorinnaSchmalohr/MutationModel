# Sys.setenv(PATH = "/data/anaconda3/bin/") # When error sh: 1: bedtools: not found apperas (hard-coding path is not optimal but solved it for now)


mapPredictors = function(x, posFile,nthreads = 1, chrL = "data/processedData/chrLengths.RData"){
   require(GenomicRanges)
   require(vcfR)
   require(stringr)
   posBed = read.table(posFile, sep = "\t", stringsAsFactors=F)
   
   posGR = GRanges(seqnames = posBed$V1, 
                   ranges = IRanges(start = posBed$V2+1, width = 1))
   fileCheck = sapply(x[,10], file.exists)
   if(any(!fileCheck)){
      stop(paste0(c("files not found: ", x[,10][!fileCheck]), 
                   collapse=" "))
   }
   res = apply(x,1,function(y){
      cat(y[3]," ")
      sFile = y[10]
      type = y[4]
      range = as.integer(y[5])
      measure = y[6]
      colName = y[7]
      if(type == "BW"){
         vals = getValFromBW(posFile=posFile, bwFile=sFile,chrL = chrL,
                             range = range, measure=measure)
      } else if (type == "GR"){
         if(measure %in% c("ifany",  "nHits", "distance")){
            vals = getOverlapFromGR(posGR=posGR, grFile = sFile, chrL = chrL,
                                    range = range, measure=measure)
         } else{
            vals = getValFromGR(posGR = posGR, grFile = sFile, chrL = chrL,
                                range = range, colName=colName, measure = measure)
         }
      } else if (type == "GRL"){
         vals = getValFromGRL(posGR = posGR, grlFile = sFile, chrL = chrL,
                              range = range, colName=colName, measure = measure)
      } else if (type == "GC"){
         vals = getGCcontent(posFile = posFile, 
                             referenceFA = sFile,
                             range = range)
      } else if (type == "SC"){
         # careful: sequence context will return a multi-column result.
         vals = getSeqcontext(posFile = posFile, 
                              referenceFA = sFile,
                              range = range)
      } else if (type == "ST"){
         vals = SIFT(posBed = posBed, 
                     posFile= posFile, 
                     SIFTpath = sFile)
      }else {
         stop(paste0("unknown value in 'type' column. current row:\n",
                     paste0(y, collapse=" ")))
      }
      if(type != "BW"){
         cat("\n")
      }
      return(vals)
   }, simplify = F)
   cat("\n")
   res = as.data.frame(res)
   colnames(res) = x$abbreviation
   res = sapply(colnames(res), function(i){
     vals = res[[i]]
     method = x$transform[x$abbreviation == i]
     if(is.na(method)){
       vals=vals
     } else  if(method == "sqrt"){
       vals=sqrt(vals)
     } else if(method == "log"){
       vals = log(vals+0.000001)
     } else if(method == "missing0"){
       vals[is.na(vals)] = 0
     } 
     vals[is.na(vals)] = mean(vals, na.rm = T)
     return(vals)
   })
   res = as.data.frame(res)
   return(res)
}

#' SIFTpath = "data/rawdata/SIFT_GRCh37.74/GRCh37.74"
#' posFile = "data/procData/breast/Muts.bed"
#' posBed= read.table(posFile, sep = "\t", stringsAsFactors=F)
SIFT = function(posBed, posFile, SIFTpath){
   # create VCF
   Muts = t(sapply(strsplit(posBed$V4, split="_"), function(x){
      x[x!=""]
   }))
   Muts = as.data.frame(Muts)
   Muts[Muts == "NA"] = NA
   colnames(Muts) = c("chr", "pos", "ref", "alt", "context", "mutated")
   origOrder = paste0(substr(Muts$chr, 4,10), "_",Muts$pos)
   
   bases = c("A", "T", "C", "G")
   vcffile = lapply(1:nrow(Muts), function(i){ 
      CHROM = as.numeric(strsplit(x = Muts$chr[i], split = "chr")[[1]][2])
      POS = as.integer(Muts$pos[i])
      REF = Muts$ref[i]
      ID = paste0(CHROM,"_",POS)
      QUAL = 100
      FILTER = "PASS"
      INFO = paste0("ST=+")
      
      alleles = bases[-which(bases == REF)] # all alternativ alleles 
      chr3 = rep(x = CHROM, 3)
      pos3 = rep(x = POS, 3)
      ref3 = rep(x = REF, 3)
      id3 = rep(x = ID, 3)
      qual3 = rep(x = QUAL, 3)
      filter3 = rep(x = FILTER, 3)
      info3 = rep(x = INFO, 3)
      return(cbind(CHROM = chr3, POS = pos3, ID = id3, REF = ref3, 
                   ALT = alleles, QUAL = qual3, FILTER = filter3, INFO = info3))
     
   })
   vcffile = as.data.frame(do.call(rbind, vcffile))
   # Sort by chrom and position
   vcffile = vcffile[order(as.numeric(vcffile$CHROM), 
                           as.numeric(vcffile$POS)),] 
   # Write Meta-info lines
   vcfPath = tempfile(tmpdir = "temp")
   fileConn<-file(paste0(vcfPath, ".vcf"))
   writeLines(c("##fileformat=VCFv4.0", 
                paste0("##source=",posFile), 
                "##reference=GRCh37", 
                "##INFO=<ID=ST,Number=1,Type=String,Description=Strand>"), fileConn)
   close(fileConn)
   # Write VCF file columns 
   write.table(vcffile, file = paste0(vcfPath, ".vcf"), 
               append = T, sep = "\t", quote = F,
               col.names = F, row.names = F)
   # run SIFT
   cmd = paste0("java -jar lib/SIFT4G_Annotator.jar -c ",
                " -i ", vcfPath, ".vcf",
                " -d ", SIFTpath,
                " -r ", vcfPath, "_SIFToutput")
   system(command = cmd, ignore.stdout = T) 
   
   # Read results
   res = read.vcfR(file = list.files(path = paste0(vcfPath, "_SIFToutput/"),
                                     pattern = ".vcf", full.names=T), verbose=F)
   unlink(x = c(paste0(vcfPath, ".vcf"), paste0(vcfPath, "_SIFToutput")),
          recursive = T)
   res = as.data.frame(res@fix)
   
   # Turn Variant type output in a score 
   res$VarType = do.call(rbind,strsplit(res$V8, split = "|", fixed = T))[,6]
   res$VarScores = sapply(res$VarType, function(x){
      if(x == "NONCODING" | is.na(x)){ # If NA or Noncoding
         return(0)
      } else if(x == "SYNONYMOUS"){ 
         return(0)
      } else if(x == "NONSYNONYMOUS"){
         return(0.5)
      } else { # If START-LOST, STOP-GAIN or STOP-LOSS
         return(1)
      }
   })
   
   scoreByID = sapply(split(res$VarScores, res$V3), mean, na.rm = T)
   output = scoreByID[origOrder]
   return(output)
}

# Extract value from BigWig ####
#' @title Extract value from BigWig
#'
#' @description
#' A function to extract values from a bwFile for a given file of
#' genomic positions.
#'
#' @param posFile a character string indicating a file path to a bed file. 
#' The bed file should only contain ranges of length one (i.e. single basepairs).
#' This function might also work for other beds, but was not tested for that.
#' @param bwFile a character string indicating a file path to a bigWig file.
#' @param range NULL or a positive integer, representing the range around
#' the positions in posFile that should be considered. Passed on as
#' \code{-sampleAroundCenter} argument to bigWigAverageOverBed.
#' @param measure one of 'covered', 'sum', 'mean0', or 'mean' (default). 
#' covered: number of bases within posFile covered by bigWig.
#' sum: sum of values over all bases covered.
#' mean0: average over bases with non-covered bases counting as zeroes.
#' mean: average over just covered bases.
#' @return A vector of values. One value for each range indicated in the 
#' bed \code{posFile}, containing read-outs from the bigWig file. Ergo, length
#' of the returned vector corresponds to the number of rows in the bed file.
#' 
#' @export
#'
#' @examples
# load("data/rdata/breast/Muts.RData")
# vals = getValFromBW(posFile = "data/procData/breast/MutsWithIDs.bed",
#                     bwFile = "data/predictors/DNAaccessibility/ovary.bigWig", 
#                     measure="mean0", range=50)
getValFromBW = function(posFile, bwFile, range = NA, measure, chrL){
   measure = match.arg(measure, choices=c('covered', 'sum', 'mean0', 'mean'))
   stopifnot("argument posFile must be provided" = !missing(posFile), 
             "argument bwFile must be provided" = !missing(bwFile),
             "argument posFile must be a character string indicating a file path" = length(posFile) == 1 & is.character(posFile),
             "argument bwFile must be a character string indicating a file path" = length(bwFile) == 1 & is.character(bwFile),
             "file in posFile does not exist" = file.exists(posFile),
             "file in bwFile does not exist" = file.exists(bwFile),
             "argument range has to be a non-negative number or NA" = ifelse(is.na(range), yes=T, no=range >= 0))
   # check if bwFile file type is correct
   pos <- regexpr("\\.([[:alnum:]]+)$", bwFile)
   ending = ifelse(pos > -1L, substring(bwFile, pos + 1L), "")
   if(!ending %in% c("bw", "bigWig", "bigwig")){
      warning(paste0("are you certain that bwFile points to a bigWig? ",
                     "File ending is '.", ending, "'. Continuing anyways."))
   } 
   # check if posFile file type is correct
   pos <- regexpr("\\.([[:alnum:]]+)$", posFile)
   ending = ifelse(pos > -1L, substring(posFile, pos + 1L), "")
   if(!ending %in% c("bed")){
      warning(paste0("are you certain that posFile points to a bed? ",
                     "File ending is '.", ending, "'. Continuing anyways."))
   } 
   load(chrL)
   if(!is.na(range)){
      temp = read.table(posFile)
      temp$V2 = as.integer(pmax(temp$V2-(range),1))
      temp$V3 = as.integer(pmin(temp$V3+(range), chrLengths[temp$V1]))
      tmpfile = tempfile(tmpdir = "temp", fileext = ".bed")
      write.table(temp, file = tmpfile, quote=F, row.names=F, col.names=F)
      cmd = paste0("lib/bigWigAverageOverBed ",
                   bwFile, " ",
                   tmpfile, " stdout")
      
   } else{
      cmd = paste0("lib/bigWigAverageOverBed ",
                   bwFile, " ",
                   posFile, " ", "stdout")
   }
   
   res = read.table(text = system(cmd, intern = T), stringsAsFactors=F)
   colnames(res) = c("name", "size", "covered", "sum", "mean0", "mean")
   if(!is.na(range)){
     file.remove(tmpfile)
   }
   # The output columns are:
   # name - name field from bed, which should be unique
   # size - size of bed (sum of exon sizes
   # covered - # bases within exons covered by bigWig
   # sum - sum of values over all bases covered
   # mean0 - average over bases with non-covered bases counting as zeroes
   # mean - average over just covered bases
   # option -sampleAroundCenter=100 --> range around position
   if(measure == "mean0"){
     return(res[,measure])
   } else if(measure == "mean"){
     return(res[,"sum"]/res[,"covered"]) # so that non-covered bases are actually NA
   } else{
     stop("non-valid measure type")
   }
   
}
#####


# extract value from GR #####
#' @title Extract value from GenomicRanges object
#'
#' @description
#' A function to extract values from a bwFile for a given file of
#' genomic positions.
#'
#' @param posGR a GenomicRanges object. This object should only contain 
#' ranges of length one (i.e. single basepairs).
#' @param grFile a character string indicating a file path to an .RData file 
#' containing a  GenomicRanges object.
#' @param range NA or a positive integer, representing the range around
#' the positions in posFile that should be considered. The ranges in posGR 
#' are extended in each direction by \code{range} before computing overlaps.
#' @param colName a character indicating the column name of the metadata column
#' of the GR object in grFile that should be returned. If NA or not among the
#' column names of the metadata, the first  metadata column is returned, with
#' a warning for the latter case.
#' @param measure either 'mean' or 'mean0'.
#' mean0: average over bases with non-covered bases counting as zeroes.
#' mean: average over just covered bases.
#' @return A vector of values. One value for each range indicated in the 
#' bed \code{posGR}, containing read-outs from the grFile ranges. Ergo, length
#' of the returned vector corresponds to the number of ranges in posGR.
#' 
#' @export
#'
#' @examples
# posGR = get(load("data/rdata/breast/Muts.RData"))
# posGR = GRanges(seqnames = posGR$chr, ranges = IRanges(start = posGR$pos, width = 1))
# temp = getValFromGR(grFile = "data/predictors/cancerExpression/BRCA.RData",
#                     posGR = posGR,
#                     range = NA,
#                     colName="score",
#                     measure = "mean0")
getValFromGR = function(posGR, grFile, range = NA, chrL,
                        colName="score", measure){
   libGR = require(GenomicRanges)
   stopifnot("GenomicRanges package must be installed" = libGR,
             "argument posGR must be provided" = !missing(posGR), 
             "argument grFile must be provided" = !missing(grFile),
             "argument posGR must be a GenomicRanges object" = is(posGR,"GRanges"),
             "posGR includes ranges longer than 1. Are these really mutation positions?" = max(width(posGR)) == 1,
             "argument grFile must be a character string indicating a file path" = length(grFile) == 1 & is.character(grFile),
             "file in grFile does not exist" = file.exists(grFile),
             "argument range has to be a non-negative number or NA" = ifelse(is.na(range), yes=T, no=range >= 0),
             "measure has to be 'mean' or 'mean0'" = measure %in% c("mean", "mean0"))
   gr = get(load(grFile))
   stopifnot("file in grFile is not a GR object"= is(posGR,"GRanges"))
   # check that colName is in colnames
   scoresAvail= colnames(mcols(gr))
   if(!colName %in% scoresAvail){
      warning(paste0("'", colName, "' is not a metacolumn in grFile. ",
                     "Using first metacolumn '",scoresAvail[1],  "' instead."))
      colName=scoresAvail[1]
   }
   if(measure == "mean"){
      default = NA
   }else if(measure == "mean0"){
      default = 0
   }
   if(is.na(range)){
      # easy case: posGR has only 1bp positions
      ov = findOverlaps(query=posGR, subject=gr, select = "all")
      if(length(ov) == 0){
         return(rep(default, length(posGR)))
      }
      map = split(mcols(gr)[,colName][subjectHits(ov)], queryHits(ov))
      res = sapply(map, mean, na.rm = T)
      res = res[as.character(1:length(posGR))]
      res[is.na(res)] = default
      return(unname(res))
      # names(res)
   } else{
      load(chrL)
      # hard case: posGR is extended by range
      start(posGR) = pmax(start(posGR)-range,1)
      end(posGR) = pmin(end(posGR) + range , chrLengths[as.character(seqnames(posGR))])
      
      # because of memory issues, process bit by bit (max 50k pos)
      splitVar = gl(ceiling(length(posGR)/50000), k= 50000, length = length(posGR))
      posGRsplit = split(posGR, splitVar)
      res = lapply(posGRsplit, function(subPosGR){
         ov = findOverlaps(query=subPosGR, subject=gr, select = "all")
         valMap = split(mcols(gr)[,colName][subjectHits(ov)], queryHits(ov))
         overlaps = pintersect(gr[subjectHits(ov)], subPosGR[queryHits(ov)])
         if(length(overlaps) == 0){
            return(rep(default, length(subPosGR)))
         }
         widthMap = split(width(overlaps), queryHits(ov))
         res = sapply(seq_len(length(valMap)), function(i){
            v = valMap[[i]]
            if(length(v)==1 & is.na(default)) {
               return(v)
            }else{
               l = widthMap[[i]]
               maxW = max(sum(l),width(subPosGR)[i])
               v=c(v,default)
               l=c(l,maxW-sum(l))
               return(sum(v*l/maxW, na.rm = T))
            }
         })
         names(res) = names(valMap)
         res = res[as.character(1:length(subPosGR))]
         res[is.na(res)] = default
         return(unname(res))
      })
      res = do.call(c, res)
      # return(res)
   }
}
#####


# check if overlap with GR #####
#' @title Extract value from GenomicRanges object
#'
#' @description
#' A function to extract values from a genomicRanges file for a given 
#' GenomicRanges object with single-basepair genomic positions.
#'
#' @param posGR a GenomicRanges object. This object should only contain 
#' ranges of length one (i.e. single basepairs).
#' @param grFile a character string indicating a file path to an .RData file 
#' containing a  GenomicRanges object.
#' @param range NA or a positive integer, representing the range around
#' the positions in posFile that should be considered. The ranges in posGR 
#' are extended in each direction by \code{range} before computing overlaps.
#' @param measure either "ifany",  "nHits", or "distance". 
#' @return A vector of values.
#' One value for each range indicated in \code{posGR}, containing read-outs 
#' from the grFile ranges, depending on the variable measure. When 
#' \code{measure = "ifany"}, the return value will be a binary indicator 
#' whether each position in posGR overlapped grFile's ranges at all. When 
#' \code{measure = "nHits"}, the return value will be a count of such overlaps.
#' 
#' @export
#'
#' @examples
# posGR = get(load("data/rdata/breast/Muts.RData"))
# posGR = GRanges(seqnames = posGR$chr, ranges = IRanges(start = posGR$pos, width = 1))
# grFile = "data/predictors/telomeres.RData"
# temp = getOverlapFromGR(grFile = "data/predictors/cancerExpression/BRCA.RData",
#                         posGR = posGR,
#                         range = NULL,
#                         measure="ifany")
getOverlapFromGR = function(posGR, grFile, range = NA, chrL,
                            measure=c("ifany", "nHits", "distance")){
   libGR = require(GenomicRanges)
   stopifnot("GenomicRanges package must be installed" = libGR,
             "argument posGR must be provided" = !missing(posGR), 
             "argument grFile must be provided" = !missing(grFile),
             "argument posGR must be a GenomicRanges object" = is(posGR,"GRanges"),
             "posGR includes ranges longer than 1. Are these really mutation positions?" = max(width(posGR)) == 1,
             "argument gr must be a character string indicating a file path" = length(grFile) == 1 & is.character(grFile),
             "file in grFile does not exist" = file.exists(grFile),
             "argument range has to be a non-negative number or NULL" = ifelse(is.na(range), yes=T, no=range >= 0),
             "measure has to either 'ifany' or 'nHits'" = measure %in% c("ifany", "nHits", "distance"))
   gr = get(load(grFile))
   stopifnot("file in grFile is not a GR object"= is(posGR,"GRanges"))
   
   if(!is.na(range)){
      load(chrL)
      start(posGR) = pmax(start(posGR)-range,1)
      end(posGR) = pmin(end(posGR) + range , chrLengths[as.character(seqnames(posGR))])
   }
   if(measure =="ifany"){
      res = as.integer(!is.na(findOverlaps(query=posGR, subject=gr, select = "first")))
   } else if(measure == "nHits"){
      res = countOverlaps(query=posGR, subject=gr)
   } else{
      res = distanceToNearest(posGR, gr)
      res = elementMetadata(res)@listData$distance
   }
   return(res)
} 
#####

# extract value from GRL #####
#' @title Extract value from GenomicRangesList object
#'
#' @description
#' A function to extract values from a genomicRanges file for a given 
#' GenomicRanges object with single-basepair genomic positions.
#'
#' @param posGR a GenomicRanges object. This object should only contain 
#' ranges of length one (i.e. single basepairs).
#' @param grlFile a character string indicating a file path to an .RData file 
#' containing a  GenomicRangesList object.
#' @param range NA or a positive integer, representing the range around
#' the positions in posFile that should be considered. The ranges in posGR 
#' are extended in each direction by \code{range} before computing overlaps.
#' @param colName a character indicating the column name of the metadata column
#' of the GR object in grFile that should be returned. If NA or not among the
#' column names of the metadata, the first  metadata column is returned, with
#' a warning for the latter case.
#' @param measure either "mean" or "mean0". 
#' mean0: average over bases with non-covered bases counting as zeroes.
#' mean: average over just covered bases.
#' @return A vector of values.
#' One value for each range indicated in posGR. The values are extracted
#' from each entry of the GR-list in grlFile and averaged.
#' 
#' @export
#'
#' @examples
#' grlFile = "data/predictors/HiC_Encode/brain_comp.RData"
#' library(GenomicRanges)
#' posGR = get(load("data/rdata/breast/Muts.RData"))
#' posGR = GRanges(seqnames = posGR$chr, ranges = IRanges(start = posGR$pos, width = 1))
#' temp = getValFromGRL(posGR=posGR, grlFile=grlFile)
getValFromGRL = function(posGR, grlFile, range = NA, chrL,
                         colName="score", measure){
   libGR = require(GenomicRanges)
   stopifnot("package GenomicRanges cannot be loaded" = libGR,
             "argument posGR must be provided" = !missing(posGR), 
             "argument grlFile must be provided" = !missing(grlFile),
             "argument posGR must be a GenomicRanges object" = is(posGR,"GRanges"),
             "posGR includes ranges longer than 1. Are these really mutation positions?" = max(width(posGR)) == 1,
             "argument grlFile must be a character string indicating a file path" = length(grlFile) == 1 & is.character(grlFile),
             "file in grlFile does not exist" = file.exists(grlFile),
             "argument range has to be a non-negative number or NA" = ifelse(is.na(range), yes=T, no=range >= 0),
             "measure has to be either 'mean' or 'mean0'" = measure %in% c("mean", "mean0"))
   grl = get(load(grlFile))
   stopifnot("file in grlFile is not a GR object"= is(posGR,"GRanges"))
   # check that colName is in colnames
   scoresAvail= unique(sapply(grl, function(x){colnames(mcols(x))}))
   if(!colName %in% scoresAvail){
      warning(paste0("'", colName, "' is not a metacolumn in grFile. ",
                     "Using first metacolumn '",scoresAvail[1],  "' instead."))
      colName=scoresAvail[1]
   }
   if(measure == "mean"){
      default = NA
   }else if(measure == "mean0"){
      default = 0
   }
   if(is.na(range)){
      # easy case: posGR has only 1bp positions
      resPerGR = sapply(grl, function(gr){
         ov = findOverlaps(query=posGR, subject=gr, select = "all")
         if(length(ov) == 0){
            return(rep(default, length(posGR)))
         }
         map = split(mcols(gr)[,colName][subjectHits(ov)], queryHits(ov))
         res = sapply(map, mean, na.rm = T)
         res = res[as.character(1:length(posGR))]
         res[is.na(res)] = default
         return(unname(res))
      })
   } else{
      # hard case: posGR is extended by range
      load(chrL)
      start(posGR) = pmax(start(posGR)-range,1)
      end(posGR) = pmin(end(posGR) + range , chrLengths[as.character(seqnames(posGR))])
      resPerGR = sapply(grl, function(gr){
         # because of memory issues, process chromosome by chromosome
         res = lapply(levels(seqnames(posGR)), function(cr){
            subPosGR = posGR[seqnames(posGR) == cr,]
            ov = findOverlaps(query=subPosGR, subject=gr, select = "all") ; gc()
            valMap = split(mcols(gr)[,colName][subjectHits(ov)], queryHits(ov))
            overlaps <- pintersect(gr[subjectHits(ov)], subPosGR[queryHits(ov)])
            if(length(overlaps) == 0){
               return(rep(default, length(subPosGR)))
            }
            widthMap = split(width(overlaps), queryHits(ov))
            res = sapply(seq_len(length(valMap)), function(i){
               v = valMap[[i]]
               if(length(v)==1 & is.na(default)) {
                  return(v)
               }else{
                  l = widthMap[[i]]
                  maxW = max(sum(l),width(subPosGR)[i])
                  v=c(v,default)
                  l=c(l,maxW-sum(l))
                  return(sum(v*l/maxW, na.rm = T))
               }
            })
            names(res) = names(valMap)
            res = res[as.character(1:length(subPosGR))]
            res[is.na(res)] = default
            names(res) = which(seqnames(posGR) == cr)
            return(res)
         })
         res = do.call(c, res)
         res = res[as.character(1:length(posGR))]
         res[is.na(res)] = default
         return(unname(res))
      })
   }
   res = rowMeans(resPerGR, na.rm= T)
   res[is.na(res)] = default
   return(res)
}

#####

# get GCcontent #####
getGCcontent = function(posFile, 
                        referenceFA,
                        range = 100){
   referenceFAI = paste0(referenceFA, ".fai")
   stopifnot("refenceFA was not found" = file.exists(referenceFA),
             "no index found for referenceFA" = file.exists(referenceFAI),
             "posFile was not found" = file.exists(posFile),
             "range has to be a non-negative integer" = range >= 0)
   cmd = paste("bedtools slop",
               "-i",posFile,
               "-g", referenceFAI,  
               "-b", range, 
               "| bedtools nuc -fi ", referenceFA, "-bed - ")
   # res = read.table(system(cmd, intern = T))
   res = read.table(text = system(cmd, intern = T), header = T, 
                    comment.char="", stringsAsFactors=F)
   return(res[,"X6_pct_gc"])
}
#####

# get sequence context #####
#posFile = "data/procData/breast/MutsWithIDs.bed"
# referenceFA = "data/rawdata/GRCh37.primary_assembly.genome.fa"
# range = 1
getSeqcontext = function(posFile, 
                         referenceFA,
                         range = 1){
   referenceFAI = paste0(referenceFA, ".fai")
   stopifnot("refenceFA was not found" = file.exists(referenceFA),
             "no index found for referenceFA" = file.exists(referenceFAI),
             "posFile was not found" = file.exists(posFile),
             "range has to be a non-negative integer" = range >= 0)
   #cmd = paste("bedtools slop",
  #             "-i",posFile,
  #             "-g", referenceFAI,  
  #             "-b", range,
  #             "| bedtools getfasta -fi", referenceFA,
  #             "-bed - -tab | cut -f 2")  #### !!!! ins't working atm this way as I needed to change the PATH variable to run bedtools !!! 
   cmd = paste("bedtools slop",
               "-i",posFile,
               "-g", referenceFAI,  
               "-b", range,
               "| bedtools getfasta -fi", referenceFA,
               "-bed - -tab")
   res = system(cmd, intern = T)
   res = substr(res, nchar(res)-2, nchar(res)) # to compensate cut -f 2
   res2 = do.call(rbind,strsplit(res, split = ""))
   res3 = data.frame(res2, stringsAsFactors=T)
   if(range > 1){
      colnames(res3) = c(paste0("preceding_" , range:1,"_base"),
                         "ref", 
                         paste0("following_", 1:range, "_base"))
   } else{
      colnames(res3) = c("precedingBase", "ref",
                         "followingBase")
   }
   
   return(res3)
}
#####

