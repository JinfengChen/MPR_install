snpSet <- sort(rownames(na.omit(allele.MPR)))
ids <- match(snpSet,rownames(myBaseData))
table(is.na(geno.data <- base2Geno(myBaseData[ids,],allele.MPR[ids,])))

geno.data.cr <- geno.data
SNPbyChr <- split(1:nrow(geno.data),substr(rownames(geno.data),1,2))
for(chr in names(SNPbyChr)){
  cat("\r",chr)
  ids <- SNPbyChr[[chr]]
  geno.data.chr <- correctGeno(geno.data[ids,],correct.FUN=correctFUNHMM,hmmFUN=hmm.vitFUN.rils,win=10,geno.probability=c(0.4975, 0.4975,0.005),seqERR=0.0106,minInterval=1)
  geno.data.cr[ids,] <- geno.data.chr
}

## lump genotypes to bin
str(geno.data.subset <- geno.data.cr) ##[grep('^05',rownames(geno.data.cr)),])
geno.data.bin <- genoToBin(geno.data.subset,base.position=as.numeric(substr(rownames(geno.data.subset),1,10)),corrected=TRUE,correct.FUN=correctFUNHMM,size=250e3,num=5,win=10,seqERR=0.01,fillSmallNA=TRUE,border.method=c("none", "midpoint", "left", "right")[1],border.validate=F,minBinsize=5e3,heterozygote=FALSE) ## correctFUNHMM


str(geno.bin <- geno.data.bin[[2]])
