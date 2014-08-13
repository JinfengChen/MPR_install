`findBlockAndFilter` <-
function (x, base.position = 1:length(x), size = sum(blocks[, 
    "size"], na.rm = TRUE)/100, num = sum(blocks[, "num"], na.rm = TRUE)/100, 
    extendNA = TRUE, fillSmallNA = FALSE) 
{
    x <- as.numeric(x)
    t <- length(x)
    block.num <- 1
    block_start <- base.position[1]
    block_end <- base.position[t]
    blockNumItem <- t
    blockType <- x[1]
    x.val <- x
    x.val[is.na(x)] <- max(x, na.rm = TRUE) + 1
    x.diff <- sort(unique(c(1, which(x.val[-1] != x.val[-t]) + 
        1)))
    x.diff <- x.diff[x.diff <= t]
    block.num <- length(x.diff)
    if (block.num <= 1) 
        return(blocks <- cbind(start = block_start, end = block_end, 
            size = block_end - block_start + 1, num = blockNumItem, 
            type = blockType))
    block_start <- base.position[x.diff]
    block_end <- base.position[c(x.diff[-1] - 1, t)]
    blockNumItem <- c(diff(x.diff), t - x.diff[block.num] + 1)
    blockType <- x[x.diff]
    blocks <- cbind(start = block_start, end = block_end, size = block_end - 
        block_start + 1, num = blockNumItem, type = blockType)
    if (extendNA) {
        block.ids <- which(is.na(blocks[, "type"]) & blocks[, 
            "size"] > -1)
        block.ids.next <- block.ids[block.ids < nrow(blocks)]
        block.ids.pre <- block.ids[block.ids > 1]
        if (length(block.ids.next) > 0) 
            blocks[block.ids.next, "end"] <- blocks[block.ids.next + 
                1, "start"] - 1
        if (length(block.ids.pre) > 0) 
            blocks[block.ids.pre, "start"] <- blocks[block.ids.pre - 
                1, "end"] + 1
        blocks[, "size"] <- blocks[, "end"] - blocks[, "start"] + 
            1
    }
    blocks.margin <- block_start[-1] - block_end[-length(block_end)] - 
        1
    blocks.extentSize <- as.numeric(blocks[, "size"]) + c(blocks.margin, 
        0) + c(0, blocks.margin)
    filter.block <- (blocks.extentSize < size & as.numeric(blocks[, 
        "num"]) < num)
    if (sum(filter.block, na.rm = T) > 0) {
        blocks[filter.block, "type"] <- NA
    }
    blocks <- mergeBlocks(blocks)
    if (nrow(blocks) == 1) 
        return(blocks)
    if (fillSmallNA) {
        blocks.margin <- blocks[-1, "start"] - blocks[-nrow(blocks), 
            "end"] - 1
        blocks.extentSize <- as.numeric(blocks[, "size"]) + c(blocks.margin, 
            0) + c(0, blocks.margin)
        filter.block <- blocks.extentSize < size
        block.na.ids <- which(is.na(blocks[, "type"]) & filter.block)
        block.na.num <- length(block.na.ids)
        if (block.na.num > 0) {
            blocks[block.na.ids, "type"] <- sapply(block.na.ids, 
                function(i) {
                  if (i == 1) 
                    return(blocks[2, "type"])
                  if (i == nrow(blocks)) 
                    return(blocks[nrow(blocks) - 1, "type"])
                  if (blocks[i - 1, "type"] == blocks[i + 1, 
                    "type"]) 
                    return(blocks[i - 1, "type"])
                  blocks[i, "type"]
                })
        }
    }
    blocks <- mergeBlocks(blocks)
    if (nrow(blocks) == 1) 
        return(blocks)
    block.ids <- which(is.na(blocks[, "type"]))
    block.ids.next <- block.ids[block.ids < nrow(blocks)]
    block.ids.pre <- block.ids[block.ids > 1]
    if (length(block.ids.next) > 0) 
        blocks[block.ids.next, "end"] <- blocks[block.ids.next + 
            1, "start"] - 1
    if (length(block.ids.pre) > 0) 
        blocks[block.ids.pre, "start"] <- blocks[block.ids.pre - 
            1, "end"] + 1
    tmp.blocks <- blocks
    for (i in 2:nrow(blocks)) if (blocks[i, 1] - blocks[i - 1, 
        2] > 1) 
        tmp.blocks <- rbind(tmp.blocks, c(blocks[i - 1, "end"] + 
            1, blocks[i, "start"] - 1, blocks[i, "start"] - blocks[i - 
            1, "end"] - 1, 0, NA))
    blocks <- tmp.blocks[order(tmp.blocks[, "end"]), ]
    mergeBlocks(blocks)
}
