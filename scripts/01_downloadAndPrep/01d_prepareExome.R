library(GenomicRanges)
library(rtracklayer)


# read in reference gtf and basic filtering #####
# gtf <- import("data/rawdata/genome/gencode.v43lift37.basic.annotation.gtf.gz") # no proper handling of tags, etc.
gtf = read.table("data/rawdata/genome/gencode.v43lift37.basic.annotation.gtf.gz", sep="\t")
colnames(gtf)[1:9] = c("chr", "source", "type", "start", "end", 
                       "score", "strand", "frame", "info") ; nrow(gtf)
gtf = gtf[gtf$chr %in% paste0("chr",1:22),] ; nrow(gtf)
gtf = gtf[grep("tag basic", gtf$info),] ; nrow(gtf)
gtf = gtf[grep("artifact", gtf$info, invert = T),] ; nrow(gtf)
# gtf = gtf[grep("CCDS", gtf$info),] ; nrow(gtf) # too stringent
gtf = gtf[grep("transcript_support_level [45]", gtf$info, invert = T),] ; nrow(gtf) # low transcript support level, see here: https://www.ensembl.org/info/genome/genebuild/transcript_quality_tags.html
#####


# split genome into three regions:
# translated (coding regions), 
# transcribed (noncoding regions, e.g. introns, UTRs, lncRNAs, ...), 
# untranscribed (the rest) ######
transcripts = makeGRangesFromDataFrame(gtf[gtf$type == "transcript",] )
coding = makeGRangesFromDataFrame(gtf[gtf$type == "CDS",])
nonCoding = setdiff(x = transcripts, 
                    y = coding, 
                    ignore.strand=TRUE) 
# findOverlaps(nonCoding, coding) # test
# transcripts = reduce(transcripts) # merge adjacent regions
save(transcripts, file = "data/predictors/transcripts.RData")
# coding = reduce(coding) # merge adjacent regions
save(coding, file = "data/predictors/coding.RData")
# nonCoding = reduce(nonCoding) # merge adjacent regions
save(nonCoding, file = "data/predictors/nonCoding.RData")
#####



# filter gtf for mappability etc. #####
load("data/predictors/coding.RData")
filterSteps = c(NULL,codExons_noUTR = sum(width(coding)))
coding = reduce(coding)
filterSteps = c(filterSteps,reduce = sum(width(coding)))

# load("data/predictors/MappabilityAlign24mer_regionsToKeep.RData")
# codExons_filtered = intersect(x = coding, y = gr)
# filterSteps = c(filterSteps,no24mers = sum(width(codExons_filtered)))
# 
# load("data/predictors/MappabilityAlign40mer_regionsToKeep.RData")
# codExons_filtered = intersect(x = coding, y = gr)
# filterSteps = c(filterSteps,no40mers = sum(width(codExons_filtered)))

# load("data/predictors/MappabilityAlign100mer_regionsToKeep.RData")
# codExons_filtered = intersect(x = coding, y = gr)
# filterSteps = c(filterSteps,no100mers = sum(width(codExons_filtered)))

load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData") # remove anything ==1
codExons_filtered = setdiff(x = coding, y = gr)
filterSteps = c(filterSteps, consensus = sum(width(codExons_filtered)))

# load("data/predictors/hg19.trf.RData")  # we are actually interested in tandem repeats
# codExons_filtered = setdiff(x = codExons_filtered, y = gr)
# filterSteps = c(filterSteps, trf = sum(width(codExons_filtered)))
#
# load("data/predictors/hg19.repeatMasker.UCSC.RData")  # remove anything == 1
# codExons_filtered2 = setdiff(x = codExons_filtered, y = gr)
# filterSteps = c(filterSteps,repeatMasker = sum(width(codExons_filtered)))

par(mar = c(3,10,1,1))
barplot(filterSteps, horiz = T, las = 1)
######


# save result #####
# GR to gtf file 
export(codExons_filtered, "data/processedData/codExons_noUTR_filtered.gtf")

# GR back to data frame 
codExons_filtered = as.data.frame(DataFrame(coding))
colnames(codExons_filtered) = c("chr", "start", "end", "width", "strand")
save(codExons_filtered,file = "data/processedData/codExons_noUTR_filtered.RData")
#####


