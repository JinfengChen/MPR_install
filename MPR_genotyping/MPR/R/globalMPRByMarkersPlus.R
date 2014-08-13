`globalMPRByMarkersPlus` <-
function (baseData, markers = NULL, alleleA = NULL, numGroup = 1, 
    groupSort = FALSE, numPerm = 1, numTry = 3, numBaseStep = 50, 
    numBaseCandidateStep = numBaseStep * 2, numKnownStep = pmax(numBaseStep/5, 
        10), numKnownCandidateStep = numKnownStep * 2, useMedianToFindKnown = TRUE, 
    maxIterate = 150, maxNStep = 3, scoreMin = 0.8, saveMidData = FALSE, 
    verbose = FALSE, strSTART = "\r", strEND = "", ...) 
{
    if (is.null(alleleA)) 
        alleleA <- markers
    ALLELE.mat <- base2Allele(baseData)
    ALLELE.num <- matrix(0, nrow = nrow(ALLELE.mat), ncol = ncol(ALLELE.mat))
    dimnames(ALLELE.num) <- dimnames(ALLELE.mat)
    alleleA.base <- alleleA[match(rownames(ALLELE.mat), names(alleleA))]
    alleleA.ids <- which(!is.na(alleleA.base))
    if (numKnownCandidateStep > length(alleleA.ids)) 
        numKnownCandidateStep <- length(alleleA.ids)
    if (numKnownStep > numKnownCandidateStep) 
        numKnownStep <- numKnownCandidateStep
    j <- 0
    midData <- NULL
    for (perm in 1:numPerm) {
        if (numGroup >= 2) {
            if (groupSort) 
                groupRIL <- split(order(colSums(!is.na(baseData)), 
                  decreasing = TRUE), cut(1:ncol(baseData), numGroup))
            else groupRIL <- split(sample(ncol(baseData)), cut(1:ncol(baseData), 
                numGroup))
            names(groupRIL) <- 1:numGroup
        }
        else {
            groupRIL <- list("1" = 1:ncol(baseData))
        }
        groupRIL <- groupRIL[sapply(groupRIL, length) > 0]
        for (ids.cols in groupRIL) {
            ids.RILrows <- which(rowSums(!is.na(cbind(baseData[, 
                ids.cols]))) > 0)
            rowN <- length(ids.RILrows)
            ids.times <- rep(0, rowN)
            ids.ok <- rep(0, rowN)
            ids.candidate <- na.omit(which(ids.times < numTry & 
                ids.ok == 0)[1:numBaseCandidateStep])
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
                ids.point <- ifelse(useMedianToFindKnown == TRUE, 
                  median(ids), ids[1])
                is.known <- !is.na(alleleA.base[ids])
                if (sum(is.known) < numKnownStep) {
                  ids.known <- na.omit(alleleA.ids[order(abs(alleleA.ids - 
                    ids.point))[sample(numKnownCandidateStep, 
                    numKnownStep)]])
                  if (length(ids.known) > (numKnownStep - sum(is.known))) 
                    ids.know <- ids.known[1:(numKnownStep - sum(is.known))]
                  ids <- unique(c(ids, ids.known))
                }
                ids <- sort(ids)
                is.known <- !is.na(alleleA.base[ids])
                iResult <- localMPR(baseData[ids, ], allele.matrix = ALLELE.mat[ids, 
                  ], maxIterate = maxIterate, maxNStep = maxNStep, 
                  returnNumIterate = TRUE, verbose = 0)
                allele.matrix <- iResult[["allele"]]
                a <- allele.matrix[is.known, ]
                b <- alleleA.base[ids[is.known]]
                a1 <- colSums(a == b, na.rm = T)/length(b)
                if (sum(a1 >= scoreMin) > 0) {
                  if (a1[1] < a1[2]) 
                    allele.matrix <- allele.matrix[, c(2, 1)]
                  j <- j + 1
                  a <- ALLELE.mat[ids, ]
                  a1 <- rowSums(a == allele.matrix) == 2
                  ALLELE.num[ids, ] <- ALLELE.num[ids, ] + cbind(a1, 
                    !a1)
                  ids.ok[ids.candidate] <- 1
                }
                ids.all <- which(ids.times < numTry & ids.ok == 
                  0)
                ids.candidate <- na.omit(ids.all[1:numBaseCandidateStep])
                n <- length(ids.candidate)
                if (verbose) 
                  cat(strSTART, "perm", perm, "unprocessed", 
                    length(ids.all), j, strEND, sep = "\t")
            }
        }
        filter.exchange <- ALLELE.num[, 1] < ALLELE.num[, 2]
        tmp <- ALLELE.mat
        tmp[filter.exchange, ] <- tmp[filter.exchange, c(2, 1)]
        ALLELE.mat <- tmp
        tmp <- ALLELE.num
        tmp[filter.exchange, ] <- tmp[filter.exchange, c(2, 1)]
        ALLELE.num <- tmp
        if (saveMidData) 
            midData[[perm]] <- list(allele = ALLELE.mat, call = ALLELE.num)
    }
    invisible(list(allele = ALLELE.mat, call = ALLELE.num, midData = midData))
}
