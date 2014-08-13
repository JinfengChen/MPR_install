\name{makeEmissionFUN}
\alias{makeEmissionFUN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeEmissionFUN(errorRate = 0.01)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{errorRate}{ ~~Describe \code{errorRate} here~~ }
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
function (errorRate = 0.01) 
{
    E <- log(errorRate)
    E2 <- log(1 - errorRate)
    E3 <- log(0.5)
    function(h, x, n) {
        if (h != 3) 
            return(ifelse(h == x, E2, E))
        else return(n * E3)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
