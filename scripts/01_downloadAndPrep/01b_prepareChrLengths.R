# prepare chromosome lengths #####
chrs = sort(paste0("chr", c(1:22)))
chrLengths = read.table("data/rawdata/genome/GRCh37.primary_assembly.genome.fa.fai")
chrLengths = setNames(chrLengths[,2],chrLengths[,1])
chrLengths = chrLengths[chrs]
save(chrLengths, file = "data/processedData/chrLengths.RData")
#####
