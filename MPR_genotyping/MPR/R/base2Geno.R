`base2Geno` <-
function (baseData = NULL, allele.matrix = NULL) 
{
    if (nrow(baseData) == ncol(allele.matrix)) 
        allele.matrix <- t(allele.matrix)
    if (nrow(baseData) != nrow(allele.matrix)) 
        stop("nrow(baseData)!=nrow(allele.matrix), allele.matrix error!!!")
    genoData <- baseData
    genoData[baseData == allele.matrix[, 1]] <- 0
    genoData[baseData == allele.matrix[, 2]] <- 1
    genoData[genoData != 1 & genoData != 0] <- NA
    genoData <- matrix(as.numeric(genoData), ncol = ncol(genoData))
    dimnames(genoData) <- dimnames(baseData)
    genoData
}
