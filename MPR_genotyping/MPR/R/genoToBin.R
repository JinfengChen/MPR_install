`genoToBin` <-
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
        rownames(geno.bin) <- sprintf("%010d", bin.border)
        list(block = blocks, bin = geno.bin, border = bin.border)
    })
    geno.bin <- NULL
    for (i in 1:length(res)) geno.bin <- rbind(geno.bin, res[[i]][[2]])
    cat("Done.\n")
    list(block = res, bin = geno.bin, border = as.numeric(rownames(geno.bin)))
}
