library(Biostrings)

# all possible pentamers
bases = c("A", "C", "G", "T")
fivemers = sort(apply(expand.grid(bases, bases, bases, bases, bases),
                   1,paste0, collapse = ""))
fivemers = cbind(fivemers,
              as.character(reverseComplement(DNAStringSet(fivemers))))
fivemers = t(apply(fivemers,1,function(x){
  temp = which(substr(x,3,3) %in% c("C", "T"))
  return(c(x[temp], x[-temp]))
}))
fivemers = unique(fivemers)
rownames(fivemers) = paste(fivemers[,1], fivemers[,2], sep = "/")
fivemers = fivemers[order(substr(fivemers[,1], 3,3), 
                    substr(fivemers[,1], 2,2), 
                    substr(fivemers[,1], 4,4), 
                    substr(fivemers[,1], 1,1), 
                    substr(fivemers[,1], 5,5)),]
pent2context = setNames(nm=c(fivemers[,1], fivemers[,2]),
                        object=c(fivemers[,1], fivemers[,1]))
save(fivemers, pent2context, file = "data/processedData/pentamers.RData")
#####