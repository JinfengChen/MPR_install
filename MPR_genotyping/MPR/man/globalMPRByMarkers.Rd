\name{globalMPRByMarkers}
\alias{globalMPRByMarkers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
globalMPRByMarkers(baseData, markers = NULL, alleleA = NULL, numTry = 3, numBaseStep = 50, numBaseCandidateStep = numBaseStep * 2, numKnownStep = pmax(numBaseStep/5, 10), numKnownCandidateStep = numKnownStep * 1.5, useMedianToFindKnown = TRUE, maxIterate = 150, maxNStep = 3, scoreMin = 0.8, verbose = FALSE, strSTART = "\r", strEND = "", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseData}{ ~~Describe \code{baseData} here~~ }
  \item{markers}{ ~~Describe \code{markers} here~~ }
  \item{alleleA}{ ~~Describe \code{alleleA} here~~ }
  \item{numTry}{ ~~Describe \code{numTry} here~~ }
  \item{numBaseStep}{ ~~Describe \code{numBaseStep} here~~ }
  \item{numBaseCandidateStep}{ ~~Describe \code{numBaseCandidateStep} here~~ }
  \item{numKnownStep}{ ~~Describe \code{numKnownStep} here~~ }
  \item{numKnownCandidateStep}{ ~~Describe \code{numKnownCandidateStep} here~~ }
  \item{useMedianToFindKnown}{ ~~Describe \code{useMedianToFindKnown} here~~ }
  \item{maxIterate}{ ~~Describe \code{maxIterate} here~~ }
  \item{maxNStep}{ ~~Describe \code{maxNStep} here~~ }
  \item{scoreMin}{ ~~Describe \code{scoreMin} here~~ }
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
