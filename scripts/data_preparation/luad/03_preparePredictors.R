setwd("/cellnet/MutationModel/")
load("data/rdata/luad/Muts.RData")

# 5mers #####
pentamer = read.table("data/procData/5mers/5mers.txt", as.is = T)
pentamer = factor(pentamer[,2])
trimer = factor(substr(pentamer,start = 2, stop = 4))
septamer = read.table("data/procData/5mers/7mers.txt", as.is = T)
septamer = factor(septamer[,2])
precedingBase = factor(substr(trimer,start = 1,stop = 1))
followingBase = factor(substr(trimer,start = 3,stop = 3))
precedingDimer = factor(substr(pentamer,start = 1,stop = 2))
followingDimer = factor(substr(pentamer,start = 4,stop = 5))
save(pentamer,trimer,septamer,
     precedingBase, followingBase,
     precedingDimer, followingDimer,
     file = "data/rdata/luad/mers.RData")
rm(pentamer,trimer,septamer,
   precedingBase, followingBase,
   precedingDimer, followingDimer)
#####


# GC content #####
# gc content 5 bases
gc5base = read.table("data/procData/luad/UCSC_tracks/hg19.gc5Base.UCSC.out",
                     as.is = T, na.strings = ".")[,4]
GCcontent = sapply(c("10bp", "100bp", "1Mb"), function(f){
   print(f)
   t = read.table(paste0("data/procData/luad/GCcontent/MutsBed",f,".bed.out"), as.is = T)
   missing = which(!1:400966 %in% t[,2])
   if (length(missing) > 0) {
      t = rbind(t,cbind(0,which(!1:400966 %in% t[,2])))
   }
   t = t[order(t[,2]),]
   ref = read.table(paste0("data/procData/luad/GCcontent/MutsBed",f,".bed"), as.is = T)
   t[,1]/(ref[,3] - ref[,2])
})
GCcontent = cbind(gc5base/100, GCcontent)
colnames(GCcontent) = paste0("GCcontent_",c("5bp", "10bp", "100bp", "1kb"))
save(GCcontent, file = "data/rdata/luad/GCcontent.RData")
rm(GCcontent)
#####


# methylation #####
# see 03e2_DNAmethylation.R
#####


# conservation #####
conservation = read.table("data/procData/luad/UCSC_tracks/hg19.100way.phyloP100way.bigWig.out.bed", as.is = T)
conservation = conservation[,5]
save(conservation, file = "data/rdata/luad/conservation.RData")
#####


# GTEx eqtl positions ######
# whether variant lies within 10 bp of an eQTL in this tissue. column 4 is minimum p-value of all hits.
GTEx_eqtl = read.table("data/procData/luad/GTEx_eqtls.out", sep = "\t", na.strings = ".", as.is = T)
GTEx_eqtl = GTEx_eqtl[,5]
GTEx_eqtl = as.integer(!is.na(GTEx_eqtl))
save(GTEx_eqtl, file = "data/rdata/luad/GTEx_eqtl.RData")
rm(GTEx_eqtl)
#####

# mutational Signatures #####
# this doesn't make sense.
# s = read.table("data/rawdata/Alexandrov2013/signatures.txt", as.is = T, sep = "\t", header = T)
# load("data/rdata/luad/mers.RData")
# Muts$trimer = trimer
# s$Substitution.Type = substr(s$Substitution.Type,start = 3,stop = 3)
# 
# compl = c("T" = "A", "A" = "T", "C" = "G", "G" = "C")
# MutSigs = t(apply(Muts,1,function(x) {
#    tri = x["trimer"]
#    alt = x["alt"]
#    if (is.na(alt))
#       alt = c("A", "T", "G", "C")
#    matching = s[s$Trinucleotide == tri & s$Substitution.Type %in% alt,]
#    if (nrow(matching) < 1) {
#       tri = strsplit(tri,split = "")[[1]]
#       tri = paste(compl[tri], collapse = "")
#       alt = compl[alt]
#       matching = s[s$Trinucleotide == tri & s$Substitution.Type %in% alt,]
#    }
#    return(colMeans(matching[,-(1:3)]))
# }))
# save(MutSigs, file = "data/rdata/luad/MutSigs.RData")
# rm(compl, s,MutSigs)
#####


# ChiP-seq: CTCF, EP300, POLR2A, POLR2AphosphoS5, and histone modifications #####
meta = read.table("data/rawdata/DNAbinding/lung/metadata.tsv",
                  header = T, sep = "\t", as.is = T)
meta = meta[meta$Assembly == "hg19" & 
               meta$Output.type == "fold change over control" &
               meta$File.Status == "released" , ]
meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "insufficient",invert = T),]
files = split(x = meta$File.accession, f = meta$Experiment.target)


library(corrplot)
dir.create("fig/luad/DNAbinding/")
DNAbinding = sapply(names(files), function(n){
   IDs = files[[n]]
   res = sapply(IDs,function(f){
      t = read.table(paste0("data/procData/luad/DNAbinding/",f,".out.bed"))
      return(t[,5])
   })
   png(paste0("fig/luad/DNAbinding/",n,"_corplot.png"),
       width = 500, height = 500)
   par(mar = c(5,4,10,2) + 0.1)
   corrplot(cor(cbind(res,"mean" = rowMeans(res))), title = paste0("\n\n",n))
   dev.off()
   return(rowMeans(res))
})

colnames(DNAbinding) = do.call(rbind,strsplit(colnames(DNAbinding), 
                                              split = "-"))[,1]
save(DNAbinding,file = "data/rdata/luad/DNAbinding.RData")
rm(DNAbinding, files,meta)

# histone marks
marks = c("H3k27ac", "H3k4me1", "H3k4me3")
lines = c("Gm12878", "H1hesc", "Hsmm", "Huvec", "K562", "Nhek", "Nhlf")
for (l in lines) {
   print(l)
   tab = sapply(marks, function(m){
      t = read.table(paste0("data/procData/luad/UCSC_tracks/wgEncodeBroadHistone",
                            l,m,
                            "StdSig.bigWig.out.bed"))
      return(t[,5])
   })
   name = paste0(l,"_histoneMarks")
   assign(name, tab)
}
save(list = paste0(lines,"_histoneMarks"), file = "data/rdata/luad/UCSC_histones.RData")
#####



# nucleosome binding #####
Nsome = read.table("data/procData/luad/UCSC_tracks/wgEncodeSydhNsomeGm12878Sig.bigWig.out.bed")
Nucleosome = Nsome[,5]
save(Nucleosome, file = "data/rdata/luad/Nucleosome.RData")
rm(Nsome, Nucleosome)
#####


# TF binding sites (ETS family and other) #####
tfbs = read.table("data/procData/luad/UCSC_tracks/wgEncodeRegTfbsClusteredV3.out",
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
save(ETS_BS,TF_BS, file = "data/rdata/luad/Tfbs.RData")
# to do: tissue-specific TF binding
#####


# abnormal sequence structures  #####
# first, execute script 03b_sequenceStructures.sh
aPhasedRepeats = read.table("data/procData/luad/structures/aPhasedRepeats.tsv", 
                            sep = "\t", stringsAsFactors = F)
directRepeats = read.table("data/procData/luad/structures/directRepeats.tsv", 
                           sep = "\t", stringsAsFactors = F)
gQuadruplex = read.table("data/procData/luad/structures/gQuadruplex.tsv", 
                         sep = "\t", stringsAsFactors = F)
mirrorRepeats = read.table("data/procData/luad/structures/mirrorRepeats.tsv", 
                           sep = "\t", stringsAsFactors = F)
shortTandemRepeats = read.table("data/procData/luad/structures/shortTandemRepeats.tsv", 
                                sep = "\t", stringsAsFactors = F)
zDNAmotifs = read.table("data/procData/luad/structures/zDNAmotifs.tsv", 
                        sep = "\t", stringsAsFactors = F)
structures = data.frame("aPhasedRepeats" = as.integer(aPhasedRepeats[,4] != 0),
                        "directRepeats" = as.integer(directRepeats[,4] != 0),
                        "gQuadruplex" = as.integer(gQuadruplex[,4] != 0),
                        "mirrorRepeats" = as.integer(mirrorRepeats[,4] != 0),
                        "shortTandemRepeats" = as.integer(shortTandemRepeats[,4] != 0),
                        "zDNAmotifs" = as.integer(zDNAmotifs[,4] != 0))
save(structures, file = "data/rdata/luad/structures.RData")
rm(structures, aPhasedRepeats, directRepeats, gQuadruplex, mirrorRepeats, shortTandemRepeats, zDNAmotifs)
#####


# mappability and repeats #####
mappability = sapply(c("24mer", "40mer", "100mer"),function(x){
   t = read.table(paste0("data/procData/luad/UCSC_tracks/wgEncodeCrgMapabilityAlign",x,".bigWig.out.bed"))
   return(t[,5])
})
colnames(mappability) = paste0("mappability_", c("24mer", "40mer", "100mer"))
# TandemRepeatsFinder
Trf = read.table("data/procData/luad/UCSC_tracks/hg19.TandemRepeatsFinder.UCSC.out",
                 as.is = T, na.strings = ".")
Trf = as.integer(!is.na(Trf[,4]))
# repeatMasker
repeatMasker = read.table("data/procData/luad/UCSC_tracks/hg19.repeatMasker.UCSC.out",
                          as.is = T, na.strings = ".")
repeatMasker = as.integer(is.na(repeatMasker[,4]))

save(mappability, Trf,repeatMasker, file = "data/rdata/luad/mappability_repeats.RData")
#####



# load gtex data and extract lung expression  #####
gtex = read.table(file = "data/rawdata/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.gz",
                  header = T, skip = 2, sep = "\t", as.is = T)
row.names(gtex) = do.call(rbind,strsplit(x = gtex$gene_id,split = ".", fixed = T))[,1]
healthyExpr = gtex[Muts$baseGeneID,"Lung"]
healthyExpr = log2(healthyExpr + 1)
save(healthyExpr,file = "data/rdata/luad/healthyExpr.RData")
rm(gtex, healthyExpr)
#####


# load and process cancer expression from cBioportal #####
luad = read.table("data/rawdata/cancer_expr/luad_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt",
                  header = T, sep = "\t", as.is = T)
luad$Entrez_Gene_Id = as.character(luad$Entrez_Gene_Id)
exprs = luad[,3:ncol(luad)]
exprs = apply(exprs,1,median)
exprs = log2(exprs + 1)
exprs = data.frame("Entrez_Gene_Id" = luad$Entrez_Gene_Id, "exprs" = exprs, stringsAsFactors = F)
library(biomaRt)
ensembl = useMart("ensembl",dataset = "hsapiens_gene_ensembl")
IDmap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id"),
               mart = ensembl)
IDmap$entrezgene_id = as.character(IDmap$entrezgene_id)
cancerExpr = sapply(Muts$baseGeneID,function(x){
   entrezid = IDmap[IDmap$ensembl_gene_id == x,2]
   if (length(entrezid) != 1) {
      return(NA)
   } else if (is.na(entrezid)) {
      return(NA)
   } else if (!any(exprs[,1] == entrezid)) {
      return(NA)
   } else{
      return(exprs[exprs[,1] == entrezid,2])
   }
})
save(cancerExpr, file = "data/rdata/luad/cancerExpr.RData")
rm(luad,exprs, cancerExpr, IDmap, ensembl )
#####


# HiC data  #####
# PC1 used for compartment A/B prediction
HiC_compPCA_table = read.table("data/rawdata/HiC_Schmitt2016/HiC_ABcompartments_hg19_Schmitt2016_PCA.csv",
                               header = T,sep = ",",as.is = T)
HiC_compPCA_table$chr = paste0("chr", HiC_compPCA_table$chr)
# A/B labels of compartment prediction
HiC_compLabels_table = read.table("data/rawdata/HiC_Schmitt2016/HiC_ABcompartments_hg19_Schmitt2016_ABlabels.csv",
                                  header = T,sep = ",",as.is = T)
HiC_compLabels_table$chr = paste0("chr", HiC_compLabels_table$chr)
# FIRE scores (frequently interacting regions)
HiC_FIRE_table = read.table("data/rawdata/HiC_Schmitt2016/HiC_FIREscores_hg19_Schmitt2016_primaryCohort.csv",
                            header = T,sep = ",",as.is = T)
HiC_FIRE_table$chr = paste0("chr", HiC_FIRE_table$chr)
HiC_TADbound_table = read.table("data/rawdata/HiC_Schmitt2016/HiC_TADboundaries_hg19_Schmitt2016_Lung.csv",
                                header = F,sep = ",",as.is = T)
colnames(HiC_TADbound_table) = c("chr", "start", "end")
HiC = t(apply(Muts[,2:3],1, function(x){
   c = x[1]
   p = as.numeric(x[2])
   HiC_compPCA = HiC_compPCA_table[HiC_compPCA_table$chr == c & 
                                      HiC_compPCA_table$start < p & 
                                      HiC_compPCA_table$end > p, "LG"]
   if (length(HiC_compPCA) == 0) {HiC_compPCA = NA}
   HiC_compLabels = HiC_compLabels_table[HiC_compLabels_table$chr == c & 
                                            HiC_compLabels_table$start < p & 
                                            HiC_compLabels_table$end > p, "LG"]
   if (length(HiC_compLabels) == 0) {HiC_compLabels = NA}
   HiC_FIRE =  HiC_FIRE_table[HiC_FIRE_table$chr == c & 
                                 HiC_FIRE_table$start < p & 
                                 HiC_FIRE_table$end > p, "LG"]
   if (length(HiC_FIRE) == 0) {HiC_FIRE = NA}
   HiC_TADbound = nrow(HiC_TADbound_table[HiC_TADbound_table$chr == c & 
                                             HiC_TADbound_table$start < p & 
                                             HiC_TADbound_table$end > p, ])
   if (length(c(HiC_compPCA, HiC_compLabels, HiC_FIRE, HiC_TADbound)) > 4)
      print(x)
   return(c(HiC_compPCA = HiC_compPCA, 
            HiC_compLabels = HiC_compLabels,
            HiC_FIRE = HiC_FIRE, 
            HiC_TADbound = HiC_TADbound))
}))
HiC = data.frame(HiC_compPCA = as.integer(HiC[,"HiC_compPCA"]),
                 HiC_compLabels = as.integer(HiC[,"HiC_compLabels"] == "A"), 
                 HiC_FIRE = as.integer(HiC[,"HiC_FIRE"]),
                 HiC_TADbound = as.integer(HiC[,"HiC_TADbound"]))
save(HiC, file = "data/rdata/luad/HiC.RData")
rm(HiC_compPCA_table, HiC_compLabels_table, HiC_FIRE_table, HiC_TADbound_table, HiC)
#####


# replication timing and direction #####
replDir = read.table("data/procData/luad/replication/Asymtools.out", 
                     as.is = T, na.strings = ".")
colnames(replDir) = c("chr", "start", "end", "replTiming", "meanExpr",
                      "is_left", "is_right", "tx_plus", "tx_minus")
replDir$replTiming = as.numeric(replDir$replTiming)
replDir$is_left = as.integer(replDir$is_left)
replDir$is_right = as.integer(replDir$is_right)
replDirection = replDir[,c("replTiming", "is_left", "is_right")]
save(replDirection, file = "data/rdata/luad/replDirection.RData")


replDir_Koren = read.table("data/procData/luad/replication/Koren.out", as.is = T)
replDir_Koren$V8[replDir_Koren$V8 > 10] = 10
replDir_Koren$V8[replDir_Koren$V8 < (-10)] = -10
replDirection_Koren = replDir_Koren[,7:9]
colnames(replDirection_Koren) = c("Koren_replTiming", "Koren_replSlope", "Koren_replDirection")
save(replDirection_Koren, file = "data/rdata/luad/replDirection_Koren.RData")

# after applying 03c1_replicationTiming.sh and 03c2_replicationTiming.sh
files = list.files("data/procData/luad/replicationdomain/", pattern = "*out")
files = sapply(files,function(x){strsplit(x,split = "_")[[1]][1]})
replicationTiming = sapply(files,function(name){
   t = read.table(paste0("data/procData/luad/replicationdomain/",name,"_scaled.out"), 
                  stringsAsFactors = F, na.strings = ".")
   return(as.numeric(t[,4]))
})
replTiming = rowMeans(replicationTiming, na.rm = T)
save(replTiming, file = "data/rdata/luad/replTiming.RData")
rm(files, replTiming, replDir, replDir_Koren, replDirection, replDirection_Koren, replicationTiming)
#####



# DNA accessibility #####
# tissue-specific:
meta = read.table("data/rawdata/DNAaccessibility/lung/metadata.tsv", header = T, sep = "\t", stringsAsFactors = F)
meta = meta[meta$Assembly == "hg19" & meta$File.Status == "released",]
DNAaccessibility_lung = sapply(meta$File.accession, function(f){
   t = read.table(file = paste0("data/procData/luad/DNAaccessibility/",f,".bigWig.out.bed"),
                  as.is = T,sep = "\t")
   return(t[,5])
})
DNAaccessibility_lung = rowMeans(log(DNAaccessibility_lung + 1))

# DNAse from UCSC, from multiple cell lines combined
DNAaccessibility_UCSC = read.table("data/procData/luad/UCSC_tracks/wgEncodeRegDnaseClusteredV3.out",
                   as.is = T, na.strings = ".")[,4]
DNAaccessibility_UCSC[is.na(DNAaccessibility_UCSC)] = 0
save(DNAaccessibility_lung, DNAaccessibility_UCSC, file = "data/rdata/luad/DNAaccessibility.RData")
rm(DNAaccessibility_lung, meta, DNAaccessibility_UCSC)
#####

