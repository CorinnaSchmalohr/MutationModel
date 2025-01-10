library(GenomicRanges)
library(biomaRt)

# gtex #####
dir.create("data/predictors/GTEx_expression/")
# get gene positions for each gene
gtf = read.table(file = "data/rawdata/gencode.v19.genes.v7.patched_contigs.gtf", 
                 sep = "\t")
gtf = gtf[gtf$V3 == "gene",]
#gtfs are 1-based
gtf$geneID=sapply(gtf$V9, function(x){
   x = strsplit(x,split="; ", fixed=T)[[1]]
   x = x[grep("gene_id",x)]
   x = strsplit(x," ")[[1]][2]
   return(x)
})
genePos = gtf[,c(1,4,5)]
rownames(genePos) =gtf$geneID

# load gtex data and gene positions from gtex
gtex = read.table(file = "data/rawdata/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz",
                  header = T, skip = 2, sep = "\t", as.is = T)
gtexPositions = genePos[gtex$gene_id,]
row.names(gtex) = do.call(rbind,strsplit(x = gtex$gene_id,
                                         split = ".", fixed = T))[,1]
# create genomicRanges object for each tissue and save
dumpVar = sapply(colnames(gtex)[-(1:2)], function(tissue){
   tissueName = gsub(pattern="[.]", replacement="", x=tissue)
   print(tissueName)
   gr = GRanges(seqnames=paste0("chr", gtexPositions$V1), 
                ranges=IRanges(start=gtexPositions$V4, end=gtexPositions$V5),
                val = log(gtex[,tissue]+1))
   gr = sort(gr)
   save(gr, file = paste0("data/predictors/GTEx_expression/",
                          tissueName, ".RData"))
})
# create GR with all brain tissues averaged
brainExpr = gtex[,grep(colnames(gtex),pattern="Brain")]
gr = GRanges(seqnames=paste0("chr", gtexPositions$V1), 
             ranges=IRanges(start=gtexPositions$V4, end=gtexPositions$V5),
             val = rowMeans(log(brainExpr+1)))
gr = sort(gr)
save(gr, file = paste0("data/predictors/GTEx_expression/",
                       "brain_combined", ".RData"))
#####



# cancer Expression #####
dir.create("data/predictors/cancerExpression/")
# load annotated positions for entrez IDs
ensembl = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",  GRCh = 37)
IDmap <- getBM(attributes = c("entrezgene_id", 
                              "chromosome_name", "start_position", "end_position"),
               mart = ensembl)
IDmap$entrezgene_id = as.character(IDmap$entrezgene_id)
IDmap =IDmap[!is.na(IDmap$entrezgene_id),]
IDmap =IDmap[IDmap$chromosome_name %in% as.character(1:22),]
IDmap$chromosome_name = paste0("chr", IDmap$chromosome_name)
IDmap = IDmap[order(IDmap$chromosome_name, IDmap$start_position),]

# load and process cancer expression from cBioportal (cancer expression) 
filepaths = list.dirs("data/rawdata/cancer_expr", recursive=F, full.names=F)
tissueNames = do.call(rbind,strsplit(filepaths, split="_"))[,1]
dumpVar = sapply(tissueNames, function(tissue){
   print(tissue)
   expr = read.table(paste0("data/rawdata/cancer_expr/", tissue,
                            "_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt"),
                     header = T, sep = "\t", as.is = T)
   meanExpr = setNames(object=apply(expr[,-(1:2)],1,median), expr[,2])
   meanExpr = log2(meanExpr + 1)
   
   gr = GRanges(seqnames=IDmap$chromosome_name, 
                ranges=IRanges(start=IDmap$start_position, 
                               end=IDmap$end_position),
                val = unname(meanExpr[IDmap$entrezgene_id]))
   save(gr, file = paste0("data/predictors/cancerExpression/",
                          tissue, ".RData"))
})

#####

# PanCan Expression normal #####
dir.create("data/predictors/TCGAnormals_expression/")
# load annotated positions for entrez IDs
ensembl = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",  GRCh = 37)
IDmap <- getBM(attributes = c("entrezgene_id", 
                              "chromosome_name", "start_position", "end_position"),
               mart = ensembl)
IDmap$entrezgene_id = as.character(IDmap$entrezgene_id)
IDmap =IDmap[!is.na(IDmap$entrezgene_id),]
IDmap =IDmap[IDmap$chromosome_name %in% as.character(1:22),]
IDmap$chromosome_name = paste0("chr", IDmap$chromosome_name)
IDmap = IDmap[order(IDmap$chromosome_name, IDmap$start_position),]

# load and process normal expression from cBioportal (cancer expression) 
filepaths = list.dirs("data/rawdata/cancer_expr", recursive=F, full.names=F)
tissueNames = do.call(rbind,strsplit(filepaths, split="_"))[,1]
dumpVar = sapply(tissueNames, function(tissue){
   print(tissue)
   if(file.exists(paste0("data/rawdata/cancer_expr/", tissue,
                         "_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_mRNA_median_normals.txt"))){
      expr = read.table(paste0("data/rawdata/cancer_expr/", tissue,
                               "_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_mRNA_median_normals.txt"),
                        header = T, sep = "\t", as.is = T)
      meanExpr = setNames(object=apply(expr[,-(1:2)],1,median), expr[,2])
      meanExpr = log2(meanExpr + 1)
      
      gr = GRanges(seqnames=IDmap$chromosome_name, 
                   ranges=IRanges(start=IDmap$start_position, end=IDmap$end_position),
                   val = unname(meanExpr[IDmap$entrezgene_id]))
      save(gr, file = paste0("data/predictors/TCGAnormals_expression/",
                             tissue, ".RData"))
   }
})

#####

# cancer Expression old  #####
# dir.create("data/predictors/cancerExpression/")
# # load gene positions 
# gtf22 = read.table("data/rawdata/gencode.v22.annotation.gtf.gz",
#                    sep = "\t", as.is = T)
# gtf22 = gtf22[gtf22$V3=="gene",]
# geneIDs = sapply(gtf22$V9, function(x){
#    x = strsplit(x,split="; ", fixed=T)[[1]]
#    x = x[grep("gene_id",x, fixed="T")]
#    x = strsplit(x," ")[[1]][2]
#    return(x)
# })
# rownames(gtf22) = geneIDs
# 
# 
# # load and process cancer expression from cBioportal (cancer expression) #
# files = list.files("data/rawdata/cancer_expr", 
#                   recursive = F, full.names=F, pattern="_averaged.tsv")
# tissueNames = do.call(rbind,strsplit(files, split="_"))[,1]
# dumpVar = sapply(tissueNames, function(tissue){
#    print(tissue)
#    expr = read.table(paste0("data/rawdata/cancer_expr/",tissue,
#                             "_averaged.tsv"), 
#                      header = T, as.is = T, sep = "\t")
#    # data is mean of log1
#    if(any(! expr$Ensembl_ID %in% rownames(gtf22))){
#       print("caution, Ensembl IDs were missing! Check script.")
#    }
#    pos = gtf22[expr$Ensembl_ID,c(1,4,5)]
#    pos = pos[order(pos$V1, pos$V4, pos$V5),]
#    gr = GRanges(seqnames=pos$V1, 
#                 ranges=IRanges(start=pos$V4, end=pos$V5), 
#                 val = expr$meanExpression)
#    gr = sort(gr)
#    save(gr, file = paste0("data/predictors/cancerExpression/",
#                           tissue, ".RData"))
# })
# # combine COAD and READ 
# coadread = sapply(c("COAD", "READ"), function(tissue){
#    load(paste0("data/predictors/cancerExpression/",
#                tissue, ".RData"))
#    return(gr)
# })
# gr = coadread[[1]]
# meanScores = rowMeans(cbind(coadread[[1]]$val,coadread[[2]]$val))
# gr$val = meanScores
# save(gr, file = paste0("data/predictors/cancerExpression/",
#                        "COADREAD", ".RData"))
#####