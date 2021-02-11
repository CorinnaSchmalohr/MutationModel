tissues = c("luad", "breast", "skin", "colon",
            "ovary", "kidney", "prostate")
sizes = c("bins1Mb" = 1000000, "bins100kb" = 100000,
          "bins10kb" = 10000) #, "bins1kb" = 1000,"bins100bp" = 100

for(size in names(sizes)){
   print(size)
   load(paste0("data/rdata/", size, "_allChrs.RData"))
   GCcontent = read.table(paste0("data/procData/",
                                 size,"_GCcontent.out"))
   missing = which(!1:nrow(bins) %in% GCcontent[,2])
   if (length(missing) > 0) {
      missing = cbind(V1 = NA, V2 = missing)
      GCcontent = rbind(GCcontent, missing)
      GCcontent = GCcontent[order(GCcontent[,2]),]
   }
   GCcontent = GCcontent[,1]/(bins[,"end"] - bins[,"start"])
   save(GCcontent, file = paste0("data/procData/",
                                 size,"_GCcontent.RData"))
}
