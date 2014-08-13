\name{base2Geno}
\alias{base2Geno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
base2Geno(baseData = NULL, allele.matrix = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseData}{ ~~Describe \code{baseData} here~~ }
  \item{allele.matrix}{ ~~Describe \code{allele.matrix} here~~ }
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
