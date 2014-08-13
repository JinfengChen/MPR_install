`base2Allele` <-
function (baseData = NULL) 
{
    allele.matrix <- t(apply(baseData, 1, function(x) {
        x <- unique(x[!is.na(x)])
        if (length(x) == 2) 
            return(x)
        else warning("SNP sit is not biallelic!")
        c(NA, NA)
    }))
    colnames(allele.matrix) <- c("P1", "P2")
    allele.matrix
}
