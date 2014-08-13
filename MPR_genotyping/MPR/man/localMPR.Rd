\name{localMPR}
\alias{localMPR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
localMPR(baseData, allele.matrix = NULL, maxIterate = 50, returnNumIterate = FALSE, maxNStep = 5, verbose = FALSE, strEND = "\n", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseData}{ ~~Describe \code{baseData} here~~ }
  \item{allele.matrix}{ ~~Describe \code{allele.matrix} here~~ }
  \item{maxIterate}{ ~~Describe \code{maxIterate} here~~ }
  \item{returnNumIterate}{ ~~Describe \code{returnNumIterate} here~~ }
  \item{maxNStep}{ ~~Describe \code{maxNStep} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
