# BigWigFile (in rtracklayer)
# https://rdrr.io/bioc/rtracklayer/man/BigWigFile.html
# this is way too slow. #####
library(rtracklayer)
load("data/rdata/luad/Muts.RData")
MutsR = GRanges(seqnames=Muts$chr, ranges=IRanges(Muts$pos, width=1))
subMutsR = MutsR[1:10000]
sel = BigWigSelection(ranges = subMutsR)

# prepare bigWig package
library(bigWig)
bw=load.bigWig("data/rawdata/DNAbinding/breast/ENCFF143AFG.bigWig")
MutsBed = data.frame("chr" = Muts$chr, start = Muts$pos, end = Muts$pos+1)
subMuts = MutsBed[1:10000,]


system.time({bigWig_res = bed.region.bpQuery.bigWig(bw, subMuts, op = "avg")})
system.time({rtrayklayer_res = import.bw(con = "data/rawdata/DNAbinding/breast/ENCFF143AFG.bigWig",
                                         selection = sel, 
                                         as = "NumericList")})

# using commandline lib/bigWigAverageOverBed is the fastest. 

unload.bigWig(bw)
remove(bw)
#####




cat data/procData/luad/MutsWithIDs.bed | head -n 5000 > temp.bed
time lib/bigWigAverageOverBed -bedOut=temp.out.bed \
data/rawdata/DNAbinding/breast/ENCFF143AFG.bigWig \
data/procData/luad/MutsWithIDs.bed temp.out.tab


wiggletools mean data/rawdata/DNAbinding/breast/ENCFF143AFG.bigWig  \
data/rawdata/DNAbinding/breast/ENCFF529OVO.bigWig > test.wig
# data/rawdata/DNAbinding/breast/ENCFF728WWJ.bigWig
f = list.files("data/rawdata/DNAbinding/breast/")
cmd = paste0("time wiggletools mean ",
             paste(list.files("data/rawdata/DNAbinding/breast",
                              pattern=".bigWig",full.names=T), collapse=" "),
             " > test.wig")

multiBigwigSummary bins \
   --numberOfProcessors 6 \
   --region chr2 \
   --outFileName test_multiBigwigSummary.npz \
   --outRawCounts test_multiBigwigSummary.tab \
   --bin
   --bwfiles data/rawdata/DNAbinding/luad/ENCFF011TOU.bigWig  \
   data/rawdata/DNAbinding/luad/ENCFF017IRQ.bigWig


