\name{correctGeno}
\alias{correctGeno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
correctGeno(geno.data, base.position = as.numeric(sub("[^0-9]*([0-9]*)[^0-9]*", "\\1", rownames(geno.data))), correct.FUN = correctFUNHMM, minInterval = 1, verbose = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno.data}{ ~~Describe \code{geno.data} here~~ }
  \item{base.position}{ ~~Describe \code{base.position} here~~ }
  \item{correct.FUN}{ ~~Describe \code{correct.FUN} here~~ }
  \item{minInterval}{ ~~Describe \code{minInterval} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
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
function (geno.data, base.position = as.numeric(sub("[^0-9]*([0-9]*)[^0-9]*", 
    "\\1", rownames(geno.data))), correct.FUN = correctFUNHMM, 
    minInterval = 1, verbose = TRUE, ...) 
{
    i <- 0
    geno.data.cr <- apply(geno.data, 2, function(x, ...) {
        if (verbose) 
            cat("\r", i <<- i + 1)
        x.nna <- which(!is.na(x))
        ids <- sort(unique(x.nna[c(1, which(diff(x[x.nna]) != 
            0 | diff(base.position[x.nna]) >= minInterval) + 
            1)]))
        x.cr <- correct.FUN(x[ids], base.position = base.position[ids], 
            ...)
        x[x.nna] <- NA
        x[ids] <- x.cr
        x
    }, ...)
    if (verbose) 
        cat("\tDone.\n")
    geno.data.cr
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
