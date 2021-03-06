\name{globalMPRByMarkersPlus}
\alias{globalMPRByMarkersPlus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
globalMPRByMarkersPlus(baseData, markers = NULL, alleleA = NULL, numGroup = 1, groupSort = FALSE, numPerm = 1, numTry = 3, numBaseStep = 50, numBaseCandidateStep = numBaseStep * 2, numKnownStep = pmax(numBaseStep/5, 10), numKnownCandidateStep = numKnownStep * 2, useMedianToFindKnown = TRUE, maxIterate = 150, maxNStep = 3, scoreMin = 0.8, saveMidData = FALSE, verbose = FALSE, strSTART = "\r", strEND = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseData}{ ~~Describe \code{baseData} here~~ }
  \item{markers}{ ~~Describe \code{markers} here~~ }
  \item{alleleA}{ ~~Describe \code{alleleA} here~~ }
  \item{numGroup}{ ~~Describe \code{numGroup} here~~ }
  \item{groupSort}{ ~~Describe \code{groupSort} here~~ }
  \item{numPerm}{ ~~Describe \code{numPerm} here~~ }
  \item{numTry}{ ~~Describe \code{numTry} here~~ }
  \item{numBaseStep}{ ~~Describe \code{numBaseStep} here~~ }
  \item{numBaseCandidateStep}{ ~~Describe \code{numBaseCandidateStep} here~~ }
  \item{numKnownStep}{ ~~Describe \code{numKnownStep} here~~ }
  \item{numKnownCandidateStep}{ ~~Describe \code{numKnownCandidateStep} here~~ }
  \item{useMedianToFindKnown}{ ~~Describe \code{useMedianToFindKnown} here~~ }
  \item{maxIterate}{ ~~Describe \code{maxIterate} here~~ }
  \item{maxNStep}{ ~~Describe \code{maxNStep} here~~ }
  \item{scoreMin}{ ~~Describe \code{scoreMin} here~~ }
  \item{saveMidData}{ ~~Describe \code{saveMidData} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
  \item{strSTART}{ ~~Describe \code{strSTART} here~~ }
  \item{strEND}{ ~~Describe \code{strEND} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
