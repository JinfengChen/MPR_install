\name{genoToBin}
\alias{genoToBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
genoToBin(genoData, base.position = 1:nrow(genoData), corrected = FALSE, correct.FUN = correctFUNHMM, size = 250000, num = 5, fillSmallNA = TRUE, minBinsize = 0, seqERR = 0.01, heterozygote = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genoData}{ ~~Describe \code{genoData} here~~ }
  \item{base.position}{ ~~Describe \code{base.position} here~~ }
  \item{corrected}{ ~~Describe \code{corrected} here~~ }
  \item{correct.FUN}{ ~~Describe \code{correct.FUN} here~~ }
  \item{size}{ ~~Describe \code{size} here~~ }
  \item{num}{ ~~Describe \code{num} here~~ }
  \item{fillSmallNA}{ ~~Describe \code{fillSmallNA} here~~ }
  \item{minBinsize}{ ~~Describe \code{minBinsize} here~~ }
  \item{seqERR}{ ~~Describe \code{seqERR} here~~ }
  \item{heterozygote}{ ~~Describe \code{heterozygote} here~~ }
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
function (genoData, base.position = 1:nrow(genoData), corrected = FALSE, 
    correct.FUN = correctFUNHMM, size = 250000, num = 5, fillSmallNA = TRUE, 
    minBinsize = 0, seqERR = 0.01, heterozygote = FALSE, ...) 
{
    SNPbyChr <- split(1:nrow(genoData), substr(rownames(genoData), 
        1, 2))
    i.chr <- 0
    res <- lapply(SNPbyChr, function(ids) {
        i.line <- 0
        i.chr <<- i.chr + 1
        blocks <- apply(genoData[ids, ], 2, function(x) {
            cat("chr: ", i.chr, "\tline: ", i.line <<- i.line + 
                1, "\t", colnames(genoData)[i.line], "\r")
            x.nna <- !is.na(x)
            if (corrected == FALSE) 
                x.correct <- correct.FUN(x[x.nna], base.position[ids][x.nna], 
                  ...)
            else x.correct <- x[x.nna]
            blocks.mat <- findBlockAndFilter(x.correct, base.position = base.position[ids][x.nna], 
                size = size, num = num, fillSmallNA = fillSmallNA)
            blocks.mat[nrow(blocks.mat), "end"] <- max(base.position[ids])
            blocks.mat[nrow(blocks.mat), "size"] <- blocks.mat[nrow(blocks.mat), 
                "end"] - blocks.mat[nrow(blocks.mat), "start"] + 
                1
            if (heterozygote == FALSE) 
                blocks.mat[blocks.mat[, "type"] == .MPR_hetero_, 
                  "type"] <- NA
            blocks.mat <- mergeBlocks(blocks.mat)
            t(blocks.mat)
        })
        blocks.mat <- matrix(unlist(blocks, recursive = TRUE, 
            use.names = FALSE), ncol = 5, byrow = TRUE)
        bin.border <- sort(unique(as.numeric(blocks.mat[, 2])))
        cat("\rchr: ", i.chr, "\tTotal", length(bin.border), 
            "borders. ")
        bin.border <- sort(unique(bin.border[c(1, which(diff(bin.border) >= 
            minBinsize) + 1)]))
        cat(length(bin.border), "borders after filtering out bins less than", 
            round(minBinsize/1000, 1), "kb.\n")
        geno.bin <- sapply(blocks, function(blocks.line) {
            blocks.end <- rbind(blocks.line)[2, ]
            ids <- match(bin.border, blocks.end)
            ids[is.na(ids)] <- findInterval(bin.border[is.na(ids)], 
                blocks.end, rightmost.closed = FALSE) + 1
            filter <- ids > ncol(blocks.line)
            ids[filter] <- ncol(blocks.line)
            rbind(blocks.line)[5, ids]
        })
        rownames(geno.bin) <- sprintf("\%010d", bin.border)
        list(block = blocks, bin = geno.bin, border = bin.border)
    })
    geno.bin <- NULL
    for (i in 1:length(res)) geno.bin <- rbind(geno.bin, res[[i]][[2]])
    cat("Done.\n")
    list(block = res, bin = geno.bin, border = as.numeric(rownames(geno.bin)))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
