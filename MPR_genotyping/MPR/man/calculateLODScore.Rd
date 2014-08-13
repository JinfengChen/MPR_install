\name{calculateLODScore}
\alias{calculateLODScore}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
calculateLODScore(geno.data, pheno.data, mapInfo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno.data}{ ~~Describe \code{geno.data} here~~ }
  \item{pheno.data}{ ~~Describe \code{pheno.data} here~~ }
  \item{mapInfo}{ ~~Describe \code{mapInfo} here~~ }
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
function (geno.data, pheno.data, mapInfo) 
{
    geno.data.tmp <- geno.data > 0.5
    geno.data.tmp[geno.data.tmp == FALSE] <- NA
    ph1 <- matrix(t(matrix(pheno.data, nrow = length(pheno.data), 
        ncol = nrow(geno.data)))[geno.data.tmp], ncol = ncol(geno.data))
    geno.data.tmp <- geno.data < 0.5
    geno.data.tmp[geno.data.tmp == FALSE] <- NA
    ph0 <- matrix(t(matrix(pheno.data, nrow = length(pheno.data), 
        ncol = nrow(geno.data)))[geno.data.tmp], ncol = ncol(geno.data))
    g1 <- rowMeans(ph1, na.rm = T)
    n1 <- rowSums(!is.na(ph1))
    s1 <- rowVars(ph1, na.rm = T)
    g0 <- rowMeans(ph0, na.rm = T)
    n0 <- rowSums(!is.na(ph0))
    s0 <- rowVars(ph0, na.rm = T)
    n <- rowSums(!is.na(geno.data))
    ss <- (s1 * (n1 - 1) + s0 * (n0 - 1))/(n1 + n0 - 2)
    t <- (g1 - g0)/sqrt(ss/n1 + ss/n0)
    res <- data.frame(mapInfo, "-log10 P" = -log10(pt(abs(t), 
        df = pmax(1, n1 + n0 - 2), lower.tail = F, log.p = FALSE) * 
        2), t = t, d = g1 - g0, n0 = n0, n1 = n1, check.names = FALSE)
    class(res) <- c("scanone", "data.frame")
    res
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
