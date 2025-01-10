
load(file = "data/procData/chrLengths.RData")
# create Bins #####
sizes = c("bins1Mb" = 1000000, "bins100kb" = 100000, 
          "bins10kb" = 10000, "bins1kb" = 1000, "bins100bp" = 100)
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

