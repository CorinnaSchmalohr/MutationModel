# give 1 to 7 to get the different tissues
args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("luad", "breast", "skin", "colon", "ovary", "kidney", "prostate")
tissue = tissues[args]
print(tissue)




septamer = read.table(paste0("data/procData/",tissue,"/context/context.txt"), as.is = T)
septamer = factor(septamer[,2])
pentamer = factor(substr(septamer,start = 2, stop = 6))
trimer = factor(substr(septamer,start = 3, stop = 5))
precedingBase = factor(substr(trimer,start = 1,stop = 1))
followingBase = factor(substr(trimer,start = 3,stop = 3))
save(pentamer,trimer,septamer,
     precedingBase, followingBase,
     file = paste0("data/rdata/",tissue,"/mers.RData"))
#####