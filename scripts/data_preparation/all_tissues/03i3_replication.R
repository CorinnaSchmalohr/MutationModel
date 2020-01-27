koren = read.table("data/procData/Koren_hg19lifted.bed", sep = "\t", as.is = T)
koren$V2 = NULL
colnames(koren) = c("chr", "pos", "repl_timing")
chroms = paste0("chr", 1:22)
koren = koren[koren$chr %in% chroms,]
koren = koren[order(koren$chr, koren$pos),]
type2col = c("right" = "blue", "left" = "red", "unknown" = "grey")
dir.create("fig/replicationDirection/")
replDir = lapply(unique(koren$chr), function(cr){
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
slopes = unlist(sapply(replDir, function(x){x$slopes}))
dir = unlist(sapply(replDir, function(x){x$dir}))
koren$end = koren$pos
koren$slope = slopes
koren$dir = dir

koren = koren[,c("chr", "pos", "end", "repl_timing", "slope", "dir")]
koren$pos = koren$pos-1
write.table(koren,file = "data/procData/Koren_hg19lifted_withSlopes.bed",
            col.names = F,row.names = F,quote = F, sep = "\t")
