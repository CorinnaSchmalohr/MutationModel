library(IRanges)

# create Folders #####
dir.create("data/streamLined")
#####


# prepare chromosome lengths #####
chrs = sort(paste0("chr", c(1:22)))
chrLengths = read.table("data/rawdata/GRCh37.p11.genome.fa.fai")
chrLengths = setNames(chrLengths[,2],chrLengths[,1])
chrLengths = chrLengths[chrs]
#####

# create Bins #####
sizes = c("bins1Mb" = 1000000, "bins100kb" = 100000, 
          "bins10kb" = 10000)# "bins1kb" = 1000, "bins100bp" = 100
for(size in names(sizes)){
   print(size)
   bins = do.call(rbind,lapply(names(chrLengths), function(chr){
      start = as.integer(seq(1,chrLengths[chr],
                             by=sizes[size]))
      end = as.integer(c(start[2:length(start)]-1,
                         chrLengths[chr]))
      bins = data.frame(chr = chr, 
                        start = start, 
                        end = end)
   }))
   save(bins, file = paste0("data/rdata/", size, "_allChrs.RData"))
   
   # save bed file
   Bins = bins
   Bins$start = as.integer(Bins$start - 1)
   Bins$id = 1:nrow(Bins)
   write.table(Bins, file=paste0("data/procData/", size, "_allChrs.bed"), 
               quote=F, sep="\t", col.names=F, row.names=F)
}
#####

# prepare exome #####

# this script takes an annotation gtf and filters it for 
# protein-coding genes, reformats it and adds the sequences of each gene
system("scripts/streamLined/lib/prepareExomeGTF.sh")

# exomes with sequences to sample from for TNs
ref_coding = read.table("data/streamLined/gencode.v19.annotation.codingonly.withsequences.bed", as.is = T)
colnames(ref_coding) = c("name", "chr", "start", "stop","sequence")
ref_coding = ref_coding[!ref_coding$chr %in% c("chrM", "chrY", "chrX"),]
save(ref_coding, file = "data/streamLined/ref_coding.RData")

# read in reference gtf
gtf = read.table("data/rawdata/gencode.v19.annotation.gtf.gz", 
                 sep = "\t", as.is = T)
colnames(gtf)[1:9] = c("chr", "source", "feature", "start", "end", 
                       "score", "strand", "frame", "info")
# get only genes
gtf_genes = gtf[gtf$feature == "gene",]
splt = strsplit(as.character(gtf_genes$info),split = ";", fixed = T)
info = t(sapply(splt, function(x) {
   id = strsplit(x[grep("gene_id",x)], " ")[[1]][2]
   transcit_id =  strsplit(x[grep("transcript_id",x)], " ")[[1]][3]
   type = strsplit(x[grep("gene_type",x)], " ")[[1]][3]
   name =  strsplit(x[grep("gene_name",x)], " ")[[1]][3]
   transcript_type = strsplit(x[grep("transcript_type",x)], " ")[[1]][3]
   transcript_status = strsplit(x[grep("transcript_status",x)], " ")[[1]][3]
   return(c(gene_id = id,transcript_id = transcit_id,type = type,name = name, 
            transcript_type = transcript_type, transcript_status = transcript_status))
}))
gtf_genes = cbind(gtf_genes,info,stringsAsFactors = F)
gtf_genes = gtf_genes[! gtf_genes$chr %in% c("chrM", "chrY", "chrX"),]
save(gtf_genes,file = "data/streamLined/gtf_genes.RData")

gtf_exons = gtf[gtf$feature == "exon",]
save(gtf_exons,file = "data/streamLined/gtf_exons.RData")
#####




