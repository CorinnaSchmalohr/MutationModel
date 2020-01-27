dir.create("data/procData/lung/DNAbinding/")
meta = read.table("data/rawdata/DNAbinding/lung/metadata.tsv", header = T, sep = "\t", as.is = T)
fc = meta[meta$Output.type == "fold change over control" & 
             meta$Assembly == "hg19" & 
             meta$File.Status == "released" , ]
fc = fc[grep(x = fc$Audit.ERROR,pattern = "extreme",invert = T),]
fc = fc[grep(x = fc$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
fc = fc[grep(x = fc$Audit.NOT_COMPLIANT,pattern = "insufficient",invert = T),]
write.table(fc$File.accession,file = "data/procData/lung/DNAbinding/fc_files.txt", 
            quote = F,row.names = F, col.names = F)
