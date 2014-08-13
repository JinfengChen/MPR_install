`localMPR` <-
function (baseData, allele.matrix = NULL, maxIterate = 50, returnNumIterate = FALSE, 
    maxNStep = 5, verbose = FALSE, strEND = "\n", ...) 
{
    if (is.null(allele.matrix)) 
        allele.matrix <- base2Allele(baseData)
    if (nrow(baseData) == ncol(allele.matrix)) 
        allele.matrix <- t(allele.matrix)
    if (nrow(baseData) != nrow(allele.matrix)) 
        stop("nrow(baseData)!=nrow(allele.matrix), allele.matrix error!!!")
    newALLELE <- ALLELE <- allele.matrix
    genoData <- base2Geno(baseData, allele.matrix)
    newRSum <- oriRSum <- length(baseData)
    numStepSize <- 1
    numIterate <- 0
    while (numStepSize <= maxNStep) {
        numIterate <- numIterate + 1
        if (verbose) 
            cat("\r", numIterate, "\t", newRSum, "\t", numStepSize, 
                "\n")
        oriRSum <- newRSum
        loopRes <- .loopMPR(genoData, allele.matrix = ALLELE, 
            numStep = numStepSize)
        newRSum <- loopRes[[1]]
        if (oriRSum > newRSum) {
            ALLELE <- loopRes[[2]]
            genoData <- loopRes[[3]]
            numStepSize <- 1
        }
        else {
            numStepSize <- numStepSize + 1
        }
    }
    if (verbose) 
        cat("\tDone.", strEND)
    rownames(ALLELE) <- rownames(baseData)
    if (returnNumIterate) 
        list(allele = ALLELE, num_iterate = numIterate, numR = newRSum)
    else ALLELE
}
