\name{geno2Cross}
\alias{geno2Cross}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
geno2Cross(geno.data, pheno.data, cm_unit = 250000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno.data}{ ~~Describe \code{geno.data} here~~ }
  \item{pheno.data}{ ~~Describe \code{pheno.data} here~~ }
  \item{cm_unit}{ ~~Describe \code{cm_unit} here~~ }
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
function (geno.data, pheno.data, cm_unit = 250000) 
{
    myGeno.data <- geno.data
    geno.value <- sort(unique(na.omit(c(myGeno.data))))
    if (identical(as.numeric(geno.value), c(0, 1))) {
        myGeno.data <- myGeno.data + 1
        geno.value <- sort(unique(na.omit(c(myGeno.data))))
    }
    if (!identical(as.numeric(geno.value), c(1, 2))) 
        stop("the value of geno.data should be 1, 2, or NA.")
    ids.RILs <- match(colnames(myGeno.data), rownames(pheno.data))
    myGeno.site <- data.frame(Chr = as.numeric(substr(rownames(myGeno.data), 
        1, 2)), Position = as.numeric(substr(rownames(myGeno.data), 
        3, 10)))
    myCrossData <- list()
    myCrossData$geno <- lapply(split(1:nrow(myGeno.data), myGeno.site$Chr), 
        function(ids) {
            myMarkers <- myGeno.site$Position[ids]/cm_unit
            names(myMarkers) <- rownames(myGeno.data)[ids]
            myMap <- list(data = t(myGeno.data[ids, ]), map = myMarkers)
            class(myMap) <- "A"
            myMap
        })
    myCrossData$pheno <- as.data.frame(pheno.data[ids.RILs, ])
    class(myCrossData) <- c("riself", "cross")
    attr(myCrossData, "alleles") <- c("A", "B")
    myCrossData
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
