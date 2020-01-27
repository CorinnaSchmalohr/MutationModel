# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)


# replication timing and direction #####
repl = read.table(paste0("data/procData/", tissue, 
                         "/replication/Koren.out"), as.is = T)
repl$V8[repl$V8 > 10] = 10
repl$V8[repl$V8 < (-10)] = -10
replication = repl[,7:9]
colnames(replication) = c("replTiming", "replSlope", "replDirection")
save(replication, file = paste0("data/rdata/", tissue,
                                "/replication.RData"))
#####