
library(Biobase)
library(MPR)


## load SNP alleles at putative SNP sites
data(snpData)
## Or: 
## snpData <- as.matrix(read.table("MPR_data/snpData_chr05.txt"))

## load overlapping SNP alleles from low-depth sequences of one parent
## you can also use genotype results from traditional genotyping studies,
## e.g. results of SSR/RFLP markers
data(markerData)
## Or: 
## markerData <- as.matrix(read.table("MPR_data/markerData_chr05.txt",row.names=1,header=T))[,1]

## load negative SNP dataset
data(fSNP) 
## Or: 
## fSNP <- scan("MPR_data/fSNP_chr05.txt",'')

## load postitive SNP dataset
data(checkData)
## Or:
## checkData <- as.matrix(read.table("MPR_data/checkData_chr05.txt"))


#################################################################################

## select 50 SNP sites to test inference of parental genotypes
set.seed(123);myBaseData <- snpData[sample(200,50),]

## random assignments of parental genotypes to alleles will result in
## a big number of recombinations 
allele.random <- base2Allele(myBaseData)

## a big number of recombinations
NumRecomEvents(myBaseData,allele.random)

## MPR inference with maximum step size of 5
allele.MPR <- localMPR(baseData=myBaseData,maxIterate=50,maxNStep=5,showDetail=TRUE)

## should be a small number compared with random assignments above 
NumRecomEvents(myBaseData,allele.MPR)

set.seed(123)

## select 30 markers randomly
markers <- sample(names(markerData)[10:50],20)

## select SNP sites which contain the 30 markers
ids <- match(markers,rownames(snpData))
str(myBaseData <- snpData[min(ids):max(ids),])

## construct positive and negative dataset
fsnp <- fSNP
tsnp <- rownames(checkData)
str(fsnp <- intersect(fsnp,rownames(myBaseData)))
str(tsnp <- intersect(tsnp[!tsnp%in%fsnp],rownames(myBaseData)))



## global MPR aiding with marker data 
allele.MPR <- globalMPRByMarkers(myBaseData,markers=markerData,numTry=3,numBaseStep=50,
             numBaseCandidateStep=100,numMarkerStep=10,useMedianToFindKnown=TRUE,
             maxIterate=150,maxNStep=3,scoreMin=0.8,verbose=TRUE)

## genotypes at some SNP site may be missing due to concordence less than 80% (scoreMin)
table(is.na(alleleA <- allele.MPR[,1]))

## check with positive dataset. looks good but contains false SNPs and may be with some errors
ids <- match(rownames(allele.MPR),rownames(checkData))
table(checkData[ids,1]==alleleA)

## positive data
table(tmp <- tsnp%in%names(na.omit(alleleA)));mean(tmp)
## negative data
table(tmp <- fsnp%in%names(na.omit(alleleA)));mean(tmp)

## then you need to refine the MPR results
set.seed(123);system.time(all.res <- globalMPRRefine(myBaseData,alleleA=na.omit(
             allele.MPR[,1]),numGroup=238,groupSort=TRUE,numPerm=20,numTry=3,
             numBaseStep=50,numBaseCandidateStep=100,numKnownStep=30,
             numKnownCandidateStep=50,useMedianToFindKnown=TRUE,maxIterate=150,
             maxNStep=3,scoreMin=0.8,saveMidData=TRUE,verbose=TRUE))

## summarize results using Bayesian inference
perm <- 10
res <- all.res$midData[[perm]]
table(geno.res <- genotypeCallsBayes(res$call,errorRate=1e-11,eps=1e-10,maxIterate=100,
             verbose=FALSE)$type)

## reconstruct parental genotypes
allele.MPR <- res$allele
allele.MPR[geno.res==2|rowMin(res$call)>=perm,] <- NA
allele.MPR[geno.res==3,] <- allele.MPR[geno.res==3,c(2,1)]

## predicted alleles of parent 1
table(is.na(alleleA <- allele.MPR[,1]))

## check again using positive and negative dataset
ids <- match(names(alleleA),rownames(checkData))
table(checkData[ids,1]==alleleA)

## the power of the MPR algorithm. lost 1 high quality SNP but removed >90% of inferior SNPs

## positive data 
table(tmp <- tsnp%in%names(na.omit(alleleA)));mean(tmp)
## FALSE  TRUE 
##     1  1227 
## [1] 0.9991857

## negative data
table(tmp <- fsnp%in%names(na.omit(alleleA)));mean(tmp)
## FALSE  TRUE 
##    36     3 
## [1] 0.07692308



################################################################################
## correct genotypes using HMM
## Example
O <- c(1,1,1,1,1,1,1,2,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,1,2,2,2,2,2,2,2,2,1,2,2,2,2,
     2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,1,1,2,2,1,2,1,2,1,2,1,2,1,2,1,1,1,1,1,1,1,2,1,1,1,1,
     1,1,1,1,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,1)
O.pos <- 30e3*1:length(O)

O.cr <- hmm.vitFUN.rils(geno=O,position=O.pos,geno.probability=c(0.4975, 0.4975,0.005),transitionFUN =
                phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01))

rbind(raw=paste(O,collapse=''),correct=paste(O.cr,collapse=''))
##         [,1]                                                                                                           
## raw     "1111111211111121111111222222221222222221222221222222222222221211221212121212111111121111111111121111111211111"
## correct "1111111111111111111111222222222222222222222222222222222222223333333333333333111111111111111111111111111111111"

## real data
snpSet <- sort(rownames(na.omit(allele.MPR)))
ids <- match(snpSet,rownames(myBaseData))
table(is.na(geno.data <- base2Geno(myBaseData[ids,],allele.MPR[ids,])))

geno.data.cr <- geno.data
SNPbyChr <- split(1:nrow(geno.data),substr(rownames(geno.data),1,2))
for(chr in names(SNPbyChr)){
  cat("\r",chr)
  ids <- SNPbyChr[[chr]]
  geno.data.chr <- correctGeno(geno.data[ids,],correct.FUN=correctFUNHMM,hmmFUN=hmm.vitFUN.rils,
                               geno.probability=c(0.4975, 0.4975,0.005),transitionFUN=phy2get.haldane.rils,
                               emissionFUN=makeEmissionFUN(errorRate=0.0106))
  geno.data.cr[ids,] <- geno.data.chr
}


## lump genotypes to bin
str(geno.data.subset <- geno.data.cr[rowSums(!is.na(geno.data.cr))>0,])
snpSite <- as.numeric(substr(rownames(geno.data.subset),1,10))
geno.data.bin <- genoToBin(geno.data.subset,base.position=snpSite,corrected=TRUE,size=250e3,num=5,fillSmallNA=TRUE,minBinsize=5e3,heterozygote=FALSE)
## genptype bins with missing data
str(geno.bin <- geno.data.bin[[2]])

#################################################################################
library(qtl)
phenoData<-as.matrix(read.table("MPR_data/phenoData.txt",header=T,row.name=1))

## Attention!!!!! the name of SNP MUST BE started with ten numeric characters, the first two characters for chromosome and the other eight charactes for physical position of the SNP, otherwise you have to modify the code of geno2Cross (the row assigning value to myGeno.site)
summary(myCrossData <- geno2Cross(geno.bin,phenoData))
cross.data <- myCrossData

myMap <- est.map(cross.data, error.prob=0.00, map.function=c("haldane","kosambi","c-f","morgan")[1], m=0, p=0, maxit=1000, tol=1e-6, sex.sp=FALSE, verbose=FALSE, omit.noninformative=TRUE)
summary(myMap)
plot.map(myMap)

cross.data <- replace.map(cross.data, myMap)
cross.data <- calc.genoprob(cross.data,error.prob=0.00)
cross.fill <- fill.geno(cross.data,method=c("imp","argmax")[2], error.prob=0.00,map.function="haldane");
cross.fill <- calc.genoprob(cross.fill,error.prob=0.00)

id.pheno <- c("Grain_width") 
## res <- cim(cross.data,pheno.col=id.pheno,n.marcovar=3,window=10,method="em",imp.method=c("imp", "argmax")[2], error.prob=0.01)
res <- scanone(cross.fill,pheno.col=id.pheno,method="em")
summary(res,threshold=4, lodcolumn=1)
i <- which.max(res$lod)
res[(i-3):(i+3),]
plot.scanone(res,chr=5) ## vtype="S",chr=5,xlim=c(25,30)

str(geno.fill <- t(pull.geno(cross.fill))-1)
colnames(geno.fill) <- colnames(geno.bin)

## Final Bin Map
str(geno.uniq <- mergeBinMap(geno.fill))

#################################################################################
X11(width=10,height=8)
par(mfrow=c(2,1),mar=c(4.1,4.1,2.1,1.1))
sel.cols <- 1:10

## comparisons
## raw
geno.colors <- geno.data[,sel.cols]+1;geno.colors[is.na(geno.colors)] <- 0
plot(matrix(1:nrow(geno.data),nrow=nrow(geno.data),ncol=length(sel.cols)),
     matrix(rep(sel.cols,each=nrow(geno.data)),nrow=nrow(geno.data),ncol=length(sel.cols)),col=geno.colors,pch="|",
     xlab='',ylab="RIL index",mar=c(1.1,2.1,2.1,1.1))
## corrected
geno.colors <- geno.data.cr[,sel.cols]+1;geno.colors[is.na(geno.colors)] <- 0
plot(matrix(1:nrow(geno.data),nrow=nrow(geno.data),ncol=length(sel.cols)),
     matrix(rep(sel.cols,each=nrow(geno.data)),nrow=nrow(geno.data),ncol=length(sel.cols)),col=geno.colors,pch="|",
     xlab="SNP index",ylab="RIL index",mar=c(4.1,2.1,2.1,1.1))

dev.off()

sessionInfo()
## R version 2.5.0 (2007-04-23) 
## i686-pc-linux-gnu 

## locale:
## LC_CTYPE=en_US.UTF-8;LC_NUMERIC=C;LC_TIME=en_US.UTF-8;LC_COLLATE=en_US.UTF-8;LC_MONETARY=en_US.UTF-8;LC_MESSAGES
## =en_US.UTF-8;LC_PAPER=en_US.UTF-8;LC_NAME=C;LC_ADDRESS=C;LC_TELEPHONE=C;LC_MEASUREMENT=en_US.UTF-8;LC_IDENTIFICATION=C

## attached base packages:
## [1] "splines"   "tools"     "stats"     "graphics"  "grDevices" "utils"    
## [7] "datasets"  "methods"   "base"     

## other attached packages:
##        MPR genefilter   survival    Biobase 
##      "0.1"   "1.14.1"     "2.34"   "1.14.1"
