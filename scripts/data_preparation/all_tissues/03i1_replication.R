koren_hg18 = read.table("data/rawdata/Koren_et_al_Table_S2.txt", 
                   sep = "\t", as.is = T)
colnames(koren_hg18) = c("chr", "pos", "repl_timing")

koren_bed = data.frame(chr = koren_hg18$chr, start = koren_hg18$pos-1, end = koren_hg18$pos,
                       repl_timing = koren_hg18$repl_timing, stringsAsFactors = F)
koren_bed$start = as.integer(koren_bed$start)
koren_bed$end = as.integer(koren_bed$end)
koren_bed$chr = paste0("chr",koren_bed$chr)
write.table(koren_bed, file = "data/procData/Koren.bed",
            col.names = F, row.names = F, quote = F,sep = "\t")
