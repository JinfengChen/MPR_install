\name{mergeBlocks}
\alias{mergeBlocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mergeBlocks(blocks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blocks}{ ~~Describe \code{blocks} here~~ }
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
function (blocks) 
{
    if (nrow(blocks) < 2) 
        return(blocks)
    type.val <- blocks[, "type"]
    t <- length(type.val)
    type.val[is.na(blocks[, "type"])] <- as.numeric(max(blocks[, 
        "type"], na.rm = TRUE)) + 1
    block.start.ids <- sort(unique(c(1, which(type.val[-1] != 
        type.val[-t]) + 1)))
    block.start.ids <- block.start.ids[block.start.ids <= nrow(blocks)]
    block.block.num <- length(block.start.ids)
    blocks[-block.start.ids, "start"] <- NA
    blocks[block.start.ids, c("end", "num")] <- t(sapply(1:block.block.num, 
        function(i) {
            ids <- block.start.ids[i]:ifelse(i >= block.block.num, 
                nrow(blocks), block.start.ids[i + 1] - 1)
            c(blocks[ids[length(ids)], "end"], sum(blocks[ids, 
                "num"], na.rm = TRUE))
        }))
    blocks <- rbind(blocks[!is.na(blocks[, "start"]), ])
    blocks[, "size"] <- blocks[, "end"] - blocks[, "start"] + 
        1
    blocks
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
