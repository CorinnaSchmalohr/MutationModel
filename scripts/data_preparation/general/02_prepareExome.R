# read in reference exome annotation bed/fasta file with sequences #####
ref_coding = read.table("data/procData/gencode.v19.annotation.codingonly.withsequences.bed", as.is = T)
colnames(ref_coding) = c("name", "chr", "start", "stop","sequence")
ref_coding = ref_coding[!ref_coding$chr %in% c("chrM", "chrY", "chrX"),]
save(ref_coding, file = "data/rdata/ref_coding.RData")
#####

# read in reference gtf
gtf = read.table("data/rawdata/gencode.v19.annotation.gtf.gz", sep = "\t", as.is = T)
colnames(gtf)[1:9] = c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "info")
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
save(gtf_genes,file = "data/rdata/gtf_genes.RData")

gtf_exons = gtf[gtf$feature == "exon",]
save(gtf_exons,file = "data/rdata/gtf_exons.RData")

#####