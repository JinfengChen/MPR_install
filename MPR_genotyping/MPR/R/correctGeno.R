`correctGeno` <-
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
