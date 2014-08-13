\name{correctFUNHMM}
\alias{correctFUNHMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correctFUNHMM(x, base.position = as.numeric(sub("[^0-9]*([0-9]*)[^0-9]*", "\\1", names(x))), hmmFUN = hmm.vitFUN.rils, geno.probability = c(0.495, 0.495, 0.01), transitionFUN = phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{base.position}{ ~~Describe \code{base.position} here~~ }
  \item{hmmFUN}{ ~~Describe \code{hmmFUN} here~~ }
  \item{geno.probability}{ ~~Describe \code{geno.probability} here~~ }
  \item{transitionFUN}{ ~~Describe \code{transitionFUN} here~~ }
  \item{emissionFUN}{ ~~Describe \code{emissionFUN} here~~ }
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
function (x, base.position = as.numeric(sub("[^0-9]*([0-9]*)[^0-9]*", 
    "\\1", names(x))), hmmFUN = hmm.vitFUN.rils, geno.probability = c(0.495, 
    0.495, 0.01), transitionFUN = phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01), 
    ...) 
{
    x.cr <- hmmFUN(geno = x + 1, position = base.position, geno.probability = geno.probability, 
        transitionFUN = transitionFUN, emissionFUN = emissionFUN, 
        ...) - 1
    x.cr[x.cr == 2] <- .MPR_hetero_
    x.cr
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
