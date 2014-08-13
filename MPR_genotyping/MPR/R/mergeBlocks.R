`mergeBlocks` <-
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
