library(GenomicRanges)
library(rtracklayer)
dir.create("data/predictors/replication")

meta = read.table("data/rawdata/replication/metadata.tsv",
                  sep = "\t", header = T)
dumpVar = apply(meta, 1, function(x){
   if(x["view"] %in% c("Peaks", "Valleys")){
      temp = read.table(paste0("data/rawdata/replication/", 
                               x["fileName"]))
      gr = GRanges(seqnames=temp$V1, 
                   ranges=IRanges(start=temp$V2, end=temp$V3))
      save(gr, file = paste0("data/predictors/replication/", 
                             x["cell"], "_",x["view"], ".Rdata"))
   }else{
     gr = import(paste0("data/rawdata/replication/", 
                        x["fileName"]))
     vals = gr$score
     gr$score = (vals-min(vals))/(max(vals)-min(vals))
     save(gr, file = paste0("data/predictors/replication/", 
                            x["cell"], "_",x["view"], ".Rdata"))
   }
})

# BG02ES:	ESC	
# BJ:     	immortalized foreskin fibroblast -->	skin
# GM06990:	B-Lymphocyte	
# GM12801:	B-Lymphocyte	
# GM12812:	B-Lymphocyte
# GM12813:	B-Lymphocyte	
# GM12878:	B-Lymphocyte	
# HeLa-S3:	Cervix epithelial adenocarcinoma 	
# HepG2:    hepatocellular carcinoma -->	liver
# HUVEC:    Human Umbilical Vein Endothelial Cells	
# IMR90: 	fibroblasts isolated from the normal lung tissue -->	lung
# K562:  	myelogenous leukemia cell line	
# MCF-7: 	breast cancer cell line	 --> breast
# NHEK:   	Primary Normal Human Epidermal Keratinocytes -->	skin
# SK-N-SH: 	neuroblastoma -->	brain
