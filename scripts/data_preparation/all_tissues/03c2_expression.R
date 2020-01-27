# load gtex data and extract lung expression  #####
tissue2GTEx = list("luad"   = "Lung",
                   "breast" = "Breast...Mammary.Tissue",
                   "skin"   = c("Skin...Sun.Exposed..Lower.leg.",
                                "Skin...Not.Sun.Exposed..Suprapubic."),
                   "colon"  = c("Colon...Sigmoid",
                                "Colon...Transverse"),
                   "ovary"  = "Ovary",
                   "kidney" = "Kidney...Cortex",
                   "prostate" = "Prostate")
gtex = read.table(file = "data/rawdata/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz",
                  header = T, skip = 2, sep = "\t", as.is = T)
row.names(gtex) = do.call(rbind,strsplit(x = gtex$gene_id,
                                         split = ".", fixed = T))[,1]
for(tissue in names(tissue2GTEx)){
   load(paste0("data/rdata/", tissue, "/Muts.RData"))
   name = tissue2GTEx[[tissue]]
   healthyExpr = gtex[Muts$baseGeneID,name]
   healthyExpr = log2(healthyExpr + 1)
   if(!is.null(ncol(healthyExpr))){
      healthyExpr = rowMeans(healthyExpr, na.rm = T)
   }
   save(healthyExpr,file = paste0("data/rdata/", tissue, "/healthyExpr.RData"))
}
rm(gtex, healthyExpr)
#####




# load and process cancer expression from cBioportal #####
tissue2Cancerexpr = c("luad"   = "luad",
                      "breast" = "brca",
                      "skin"   = "skcm",
                      "colon"  = "coadread",
                      "ovary"  = "ov",
                      "kidney" = "kirc",
                      "prostate" = "prad")
library(biomaRt)
ensembl = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror="useast")
# ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
IDmap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
               mart = ensembl)
IDmap$entrezgene_id = as.character(IDmap$entrezgene_id)
for(tissue in names(tissue2Cancerexpr)){
   print(tissue)
   load(paste0("data/rdata/", tissue, "/Muts.RData"))
   expr = read.table(paste0("data/rawdata/cancer_expr/", tissue2Cancerexpr[tissue],
                            "_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt"),
                     header = T, sep = "\t", as.is = T)
   expr$Entrez_Gene_Id = as.character(expr$Entrez_Gene_Id)
   exprs = expr[,3:ncol(expr)]
   exprs = apply(exprs,1,median)
   exprs = log2(exprs + 1)
   exprs = data.frame("Entrez_Gene_Id" = expr$Entrez_Gene_Id, 
                      "exprs" = exprs, stringsAsFactors = F)
   cancerExpr = sapply(Muts$baseGeneID,function(x){
      entrezid = IDmap[IDmap$ensembl_gene_id == x,2]
      res = exprs[exprs$Entrez_Gene_Id %in% entrezid,2]
      return(mean(res, na.rm=T))
   })
   save(cancerExpr, file = paste0("data/rdata/", tissue, "/cancerExpr.RData"))
}

#####

