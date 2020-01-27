# make sure to use files that combined replicates 1 and 2 over each individual file.
# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)



# ChiP-seq: CTCF, EP300, POLR2A, POLR2AphosphoS5, and histone modifications #####
meta = read.table(paste0("data/rawdata/DNAbinding/", tissue, "/metadata.tsv"),
                  header = T, sep = "\t", as.is = T)
meta = meta[meta$Output.type %in% c("fold change over control", "signal") & 
             meta$Assembly == "hg19" & 
             meta$File.Status == "released" , ]
meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "insufficient",invert = T),]
rownames(meta) = meta$File.accession

files = list.files(paste0("data/procData/", tissue, "/DNAbinding/"))
files = substring(files, 1,11)
filesByTarget = split(files, f = meta[files, "Experiment.target"])

library(corrplot)
dir.create(paste0("fig/", tissue, "/DNAbinding/"), showWarnings=F)
DNAbinding = sapply(names(filesByTarget), function(target){
   ids=filesByTarget[[target]]
   res = sapply(ids, function(id){
      t = read.table(paste0("data/procData/", tissue,
                            "/DNAbinding/", id,".out.bed"))
      return(t[,5])
   })
   
   png(paste0("fig/", tissue, "/DNAbinding/",target,"_corplot.png"),
       width = 500, height = 500)
   corrplot(cor(cbind(res,"mean" = rowMeans(res))), title = paste0("\n\n",target))
   dev.off()
   return(rowMeans(res))
})

colnames(DNAbinding) = do.call(rbind,strsplit(colnames(DNAbinding), 
                                              split = "-"))[,1]
save(DNAbinding,file = paste0("data/rdata/", tissue, "/DNAbinding.RData"))
rm(DNAbinding, files,meta)

# histone marks
marks = c("H3k27ac", "H3k4me1", "H3k4me3")
lines = c("Gm12878", "H1hesc", "Hsmm", "Huvec", "K562", "Nhek", "Nhlf")
histones = sapply(marks, function(m){
   res = sapply(lines, function(l){
      t = read.table(paste0("data/procData/", 
                            tissue, "/UCSC_tracks/wgEncodeBroadHistone",
                            l,m,
                            "StdSig.out.bed"))
      return(t[,5])
   })
   res = rowMeans(res)
   return(res)
})
colnames(histones) = paste("UCSC", colnames(histones), sep = "_")
save(histones, 
     file = paste0("data/rdata/", tissue, "/UCSC_histones.RData"))
#####


# TF binding sites (ETS family and other) #####
tfbs = read.table(paste0("data/procData/", tissue,
                         "/UCSC_tracks/wgEncodeRegTfbsClusteredV3.out"),
                  as.is = T, na.strings = ".")
ETS = c("ELF", "ELF1", "ELF2", "NERF", "ELF4", "MEF", 
        "ELG", "GABPA", "GABP", 
        "ERG", "FLI1", "FEV",
        "ERF", "PE-2", "ETV3", "PE1", 
        "ESE", "ELF3", "ESE1", "ESX", "ELF5", "ESE2", "ESE3", "EHF",
        "ETS", "ETS1", "ETS2", 
        "PDEF", "SPDEF", "PSE", 
        "PEA3", "ETV4", "E1AF", "ETV5", "ERM", "ETV1", "ER81",
        "ER71", "ETV2",
        "SPI", "SPI1", "PU.1", "SPIB", "SPIC",
        "TCF", "ELK1", "ELK4", "SAP1", "SAP-1", "ELK3", "NET", "SAP2",
        "TEL", "ETV6", "ETV7", "TEL2") 
ETS_BS = as.integer(sapply(strsplit(tfbs[,4],split = ",", fixed = T), 
                           function(x){any(x %in% ETS)}))
TF_BS = as.integer(sapply(strsplit(tfbs[,4],split = ",", fixed = T), 
                          function(x){any(!is.na(x))}))
save(ETS_BS,TF_BS, file = paste0("data/rdata/", tissue, "/Tfbs.RData"))
# to do: tissue-specific TF binding
#####


