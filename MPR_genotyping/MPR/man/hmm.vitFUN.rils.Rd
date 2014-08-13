\name{hmm.vitFUN.rils}
\alias{hmm.vitFUN.rils}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
hmm.vitFUN.rils(geno, position, geno.probability, transitionFUN = phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{ ~~Describe \code{geno} here~~ }
  \item{position}{ ~~Describe \code{position} here~~ }
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
function (geno, position, geno.probability, transitionFUN = phy2get.haldane.rils, 
    emissionFUN = makeEmissionFUN(errorRate = 0.01), ...) 
{
    n.obs <- length(geno)
    n.state <- length(geno.probability)
    psi <- delta <- matrix(0, nrow = n.state, ncol = n.obs)
    n.con <- geno.cr <- numeric(n.obs)
    geno.dis <- abs(diff(position))
    n.con[1] <- 1
    g <- geno[1]
    for (i in 2:n.obs) {
        n.con[i] <- ifelse(geno[i] == g, n.con[i - 1] + 1, 1)
        g <- geno[i]
    }
    for (i in 1:n.state) delta[i, 1] <- log(geno.probability[i]) + 
        emissionFUN(i, geno[1], n.con[1])
    preProb <- numeric(n.state)
    for (t in 2:n.obs) {
        for (j in 1:n.state) {
            for (i in 1:n.state) preProb[i] <- delta[i, t - 1] + 
                log(transitionFUN(i, j, geno.dis[t - 1]))
            psi[j, t] <- which.max(preProb)
            delta[j, t] <- max(preProb) + emissionFUN(j, geno[t], 
                n.con[t])
        }
    }
    geno.cr[n.obs] <- which.max(delta[, n.obs])
    for (t in seq(n.obs - 1, 1, by = -1)) geno.cr[t] <- psi[geno.cr[t + 
        1], t + 1]
    geno.cr
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
