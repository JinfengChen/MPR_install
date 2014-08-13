\name{mergeBinMap}
\alias{mergeBinMap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mergeBinMap(geno.fill)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno.fill}{ ~~Describe \code{geno.fill} here~~ }
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
function (geno.fill) 
{
    pre <- NULL
    pre.nna <- NULL
    IsUniq <- rep(1, nrow(geno.fill))
    for (j in 1:nrow(geno.fill)) {
        cur.nna <- which(!is.na(geno.fill[j, ]))
        cur <- geno.fill[j, cur.nna]
        if (identical(pre.nna, cur.nna) && identical(pre, cur)) 
            IsUniq[j] <- 0
        else {
            pre <- cur
            pre.nna <- cur.nna
        }
    }
    geno.fill[unique(c(which(IsUniq == 1)[-1] - 1, nrow(geno.fill))), 
        ]
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
