\name{genotypeCallsBayes}
\alias{genotypeCallsBayes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genotypeCallsBayes(ALLELE.num, errorRate = 5e-04, eps = 1e-10, maxIterate = 100, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ALLELE.num}{ ~~Describe \code{ALLELE.num} here~~ }
  \item{errorRate}{ ~~Describe \code{errorRate} here~~ }
  \item{eps}{ ~~Describe \code{eps} here~~ }
  \item{maxIterate}{ ~~Describe \code{maxIterate} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
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
function (ALLELE.num, errorRate = 5e-04, eps = 1e-10, maxIterate = 100, 
    verbose = FALSE) 
{
    P0.n <- 0.5
    P0.p1 <- P0.p2 <- (1 - P0.n)/2
    ALLELE.prob <- NULL
    P.n <- 1/3
    P.p1 <- P.p2 <- (1 - P.n)/2
    numIterate <- 1
    E <- errorRate
    while (abs(P0.n - P.n) > eps && numIterate <= maxIterate) {
        P0.p1 <- P.p1
        P0.p2 <- P.p2
        P0.n <- P.n
        n <- rowSums(ALLELE.num)
        k <- ALLELE.num[, 1]
        ALLELE.prob <- cbind((1 - E)^k * E^(n - k), 0.5^n, (1 - 
            E)^(n - k) * E^k) \%*\% diag(c(P0.p1, P0.n, P0.p2))
        ALLELE.type <- rep(NA, nrow(ALLELE.prob))
        ALLELE.maxprob <- rowMax(ALLELE.prob)
        ALLELE.type[ALLELE.prob[, 1] == ALLELE.maxprob] <- 1
        ALLELE.type[ALLELE.prob[, 3] == ALLELE.maxprob] <- 3
        ALLELE.type[ALLELE.prob[, 2] == ALLELE.maxprob | ALLELE.prob[, 
            1] == ALLELE.prob[, 3]] <- 2
        P.n <- mean(ALLELE.type == 2)
        P.p1 <- P.p2 <- (1 - P.n)/2
        if (verbose) 
            cat(P.p1, P.n, P.p2, table(filter.markers <- ALLELE.type != 
                2), "\n", sep = "\t")
        numIterate <- numIterate + 1
    }
    ALLELE.prob <- ALLELE.prob/rowSums(ALLELE.prob, na.rm = TRUE)
    list(prop = c(P.p1, P.n, P.p2), prob = ALLELE.prob, type = ALLELE.type)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
