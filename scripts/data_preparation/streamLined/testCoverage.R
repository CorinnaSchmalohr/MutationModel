library(IRanges)
chrs = sort(paste0("chr", c(1:22)))
chrLengths = read.table("data/rawdata/GRCh37.p11.genome.fa.fai")
chrLengths = setNames(chrLengths[,2],chrLengths[,1])
chrLengths = chrLengths[chrs]

tissues = c("luad", "breast", "skin", "colon", "ovary", "prostate", "kidney")
t2c = c(luad = "luad", breast="brca", skin = "skcm", colon = "coadread", 
        ovary = "ov", kidney = "kirc" , prostate = "prad")
# load regions where >50% are covered
covs = lapply(tissues, function(tissue){
   temp = read.table(paste0("data/procData/", tissue, "/covered_regions.bed"))
   return(temp[temp$V1 == "chr1",])
})
#load individual coverage files
covsDetail = sapply(tissues, function(tissue){
   files = list.files(paste0("data/rawdata/coverage/", 
                             
                             t2c[tissue],"/beds/"),pattern=".bed")
   covsTissue = sapply(files, function(f){
      temp = read.table(paste0("data/rawdata/coverage/", 
                               t2c[tissue],"/beds/",f), stringsAsFactors=F)
      return(temp[temp$V1 == 1,2:3])
   }, simplify=F)
   return(covsTissue)
}, simplify=F)
save(covsDetail, file="data/temp.RData")

load("data/streamLined/gtf_exons.RData") #gtf_exons
gtf_exons = gtf_exons[gtf_exons$chr == "chr1",]
gtf_exons = data.frame(lapply(gtf_exons, function(x){
   if(is.factor(x)){return(as.character(x))} else{return(x)}
}), stringsAsFactors=F)
info = t(sapply(gtf_exons$info, function(x){
   x = strsplit(as.character(x),split = ";", fixed = T)[[1]]
   id = strsplit(x[grep("gene_id",x)], " ")[[1]][2]
   transcript_id =  strsplit(x[grep("transcript_id",x)], " ")[[1]][3]
   type = strsplit(x[grep("gene_type",x)], " ")[[1]][3]
   name =  strsplit(x[grep("gene_name",x)], " ")[[1]][3]
   transcript_type = strsplit(x[grep("transcript_type",x)], " ")[[1]][3]
   transcript_status = strsplit(x[grep("transcript_status",x)], " ")[[1]][3]
   return(c(gene_id = id,transcript_id = transcript_id,type = type,name = name,
            transcript_type = transcript_type, transcript_status = transcript_status))
}, USE.NAMES=F))
gtf_exons = cbind(gtf_exons, info)
gtf_exons$info = NULL
codExons = gtf_exons[gtf_exons$type == "protein_coding" & 
                        gtf_exons$transcript_type == "protein_coding" & 
                        gtf_exons$transcript_status == "KNOWN",] #
illumina = read.table("data/rawdata/coverage/nexterarapidcapture_exome_targetedregions.bed",
                      sep = "\t", skip =1)
illumina = illumina[illumina$V1 == "chr1",]

# visualize coverage
tcols = setNames(rainbow(7),tissues)
nFiles = 30
n = nFiles*length(covsDetail)
png("fig/coverageComparison.png",
    width=2000, height=1300, pointsize=80)
par(mar = c(3,2,0.5,0.5))
# n = sum(sapply(covsDetail, length))
plot(NA, xlim = c(0, chrLengths["chr1"]),ylim = c(1,n+1), yaxt = "n", 
     xlab = "bp", ylab = "", mgp=c(1.5,0.5,0))
i = 1
for(tissue in tissues){
   for(f in 1:nFiles){ #length(covsDetail[[tissue]]
      segments(y0=i,y1=i, col=tcols[tissue],
               x0=covsDetail[[tissue]][[f]][,1], 
               x1=covsDetail[[tissue]][[f]][,2],
               lty=1, lwd=2,  lend = 2)
      i = i+1
   }
}
segments(y0 = n+20, y1 = n+20, x0 = codExons$start, x1 = codExons$end, 
         col = 1,lty=1, lwd=2, lend = 2)
segments(y0 = n+40, y1 = n+40, x0 = illumina$V2, x1 = illumina$V3, 
         col = "grey",lty=1, lwd=2, lend = 2)
dev.off()
