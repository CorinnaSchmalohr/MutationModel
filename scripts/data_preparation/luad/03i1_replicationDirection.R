dir.create("data/procData/AsymTools")
library(R.matlab)
t = readMat("data/rawdata/AsymTools-1.0.3/reference/per_base_territories_20kb.mat")
t = t$W[,,1]
t = data.frame(t)
colnames(t) = c("chr", "start", "end", "A1", "A2", "A3", "A4",
                "repl_timing", "mean_expr", "is_left", "is_right", 
                "tx_plus", "tx_minus", "tx_A1", "txA2", "txA3", "txA4")

t$chr = paste0("chr", t$chr)
t$start = as.integer(t$start)
t$end = as.integer(t$end)
write.table(t, file = "data/procData/AsymTools/per_base_territories_20kb.bed", 
            col.names = F,row.names = F,quote = F, sep = "\t")

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

# lifted over to hg19 using UCSC tool: https://genome.ucsc.edu/cgi-bin/hgLiftOver
koren = read.table("data/procData/Koren_hg19lifted.bed", sep = "\t")
koren$V2 = NULL
colnames(koren) = c("chr", "pos", "repl_timing")

type2col = c("right" = "blue", "left" = "red", "unknown" = "grey")
dir.create("fig/replicationDirection/")
replDir = sapply(unique(koren$chr), function(cr){
   print(cr)
   dat = koren[koren$chr == cr,]
   s = splinefun(dat$pos, dat$repl_timing)
   d = s(dat$pos, deriv = 1)*1000000
   col = sapply(1:length(d), function(i){
      p = dat[i,"pos"]
      w = which(dat$pos > p+-200000 & dat$pos < p+200000 )
      x = mean(d[w])
      if(x>0.9){
         return("right")
      } else if(x<(-0.9)){
         return("left")
      }else{
         return("unknown")
      }
   })
   png(paste0("fig/replicationDirection/", cr, ".png"), 
       width = 1200, height = 500)
   plot(x = dat$pos, dat$repl_timing,
        type = "l", lty = 2)
   segments(x0 = dat$pos[1:(nrow(dat)-1)], x1 = dat$pos[2:nrow(dat)],
            y0 = dat$repl_timing[1:(nrow(dat)-1)], y1 = dat$repl_timing[2:nrow(dat)],
            col =type2col[col] )
   dev.off()
   return(list(slopes = d, dir = col))
})
koren$slope = unlist(replDir[seq(1,46,2)])
koren$dir = unlist(replDir[seq(2,46,2)])
koren$end = koren$pos
koren = koren[,c(1,2,6,3:5)]
koren$pos = koren$pos-1
write.table(koren,file = "data/procData/Koren_hg19lifted_withSlopes.bed",
            col.names = F,row.names = F,quote = F, sep = "\t")
