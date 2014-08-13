`globalMPRByMarkers` <-
function (baseData, markers = NULL, alleleA = NULL, numTry = 3, 
    numBaseStep = 50, numBaseCandidateStep = numBaseStep * 2, 
    numKnownStep = pmax(numBaseStep/5, 10), numKnownCandidateStep = numKnownStep * 
        1.5, useMedianToFindKnown = TRUE, maxIterate = 150, maxNStep = 3, 
    scoreMin = 0.8, verbose = FALSE, strSTART = "\r", strEND = "", 
    ...) 
{
    if (is.null(alleleA)) 
        alleleA <- markers
    ALLELE.mat <- matrix(NA, nrow = nrow(baseData), ncol = 2)
    rownames(ALLELE.mat) <- rownames(baseData)
    alleleA.base <- alleleA[match(rownames(ALLELE.mat), names(alleleA))]
    alleleA.ids <- which(!is.na(alleleA.base))
    if (numKnownCandidateStep > length(alleleA.ids)) 
        numKnownCandidateStep <- length(alleleA.ids)
    if (numKnownStep > numKnownCandidateStep) 
        numKnownStep <- numKnownCandidateStep
    j <- 0
    ids.RILrows <- which(rowSums(!is.na(cbind(baseData))) > 0)
    rowN <- length(ids.RILrows)
    ids.times <- rep(0, rowN)
    ids.ok <- rep(0, rowN)
    ids.candidate <- na.omit(which(ids.times < numTry & ids.ok == 
        0)[1:numBaseCandidateStep])
    n <- length(ids.candidate)
    while (n > 1) {
        if (length(ids.candidate) > numBaseStep) {
            filter.dis <- ids.candidate < (median(ids.candidate) - 
                numBaseCandidateStep)
            ids.dis <- ids.candidate[filter.dis]
            if (length(ids.dis) < numBaseStep) 
                ids.candidate <- c(ids.dis, sample(ids.candidate[!filter.dis], 
                  numBaseStep - length(ids.dis)))
            else ids.candidate <- ids.dis
        }
        ids.times[ids.candidate] <- ids.times[ids.candidate] + 
            1
        ids <- ids.RILrows[ids.candidate]
        ids.point <- ifelse(useMedianToFindKnown == TRUE, median(ids), 
            ids[1])
        is.known <- !is.na(alleleA.base[ids])
        if (sum(is.known) < numKnownStep) {
            ids.known <- na.omit(alleleA.ids[order(abs(alleleA.ids - 
                ids.point))[sample(numKnownCandidateStep, numKnownStep)]])
            if (length(ids.known) > (numKnownStep - sum(is.known))) 
                ids.know <- ids.known[1:(numKnownStep - sum(is.known))]
            ids <- unique(c(ids, ids.known))
        }
        ids <- sort(ids)
        is.known <- !is.na(alleleA.base[ids])
        iResult <- localMPR(baseData[ids, ], maxIterate = maxIterate, 
            maxNStep = maxNStep, returnNumIterate = TRUE, verbose = 0)
        allele.matrix <- iResult[["allele"]]
        a <- allele.matrix[is.known, ]
        b <- alleleA.base[ids[is.known]]
        a1 <- colSums(a == b, na.rm = T)/length(b)
        if (sum(a1 >= scoreMin) > 0) {
            if (a1[1] < a1[2]) 
                allele.matrix <- allele.matrix[, c(2, 1)]
            j <- j + 1
            a <- ALLELE.mat[ids, ]
            ids.na <- rowSums(is.na(a)) > 0
            if (sum(ids.na) > 0) 
                ALLELE.mat[ids[ids.na], ] <- allele.matrix[ids.na, 
                  ]
            ids.ok[ids.candidate] <- 1
        }
        ids.all <- which(ids.times < numTry & ids.ok == 0)
        ids.candidate <- na.omit(ids.all[1:numBaseCandidateStep])
        n <- length(ids.candidate)
        if (verbose) 
            cat(strSTART, length(ids.all), j, strEND, sep = "\t")
    }
    invisible(ALLELE.mat)
}
