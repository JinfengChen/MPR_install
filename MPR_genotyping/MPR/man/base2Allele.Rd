\name{base2Allele}
\alias{base2Allele}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
base2Allele(baseData = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{baseData}{ ~~Describe \code{baseData} here~~ }
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
function (baseData = NULL) 
{
    allele.matrix <- t(apply(baseData, 1, function(x) {
        x <- unique(x[!is.na(x)])
        if (length(x) == 2) 
            return(x)
        else warning("SNP sit is not biallelic!")
        c(NA, NA)
    }))
    colnames(allele.matrix) <- c("P1", "P2")
    allele.matrix
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
