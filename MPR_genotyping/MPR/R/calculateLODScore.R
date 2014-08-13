`calculateLODScore` <-
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
