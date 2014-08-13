`NumRecomEvents` <-
function (baseData, allele.matrix, genoData = NULL) 
{
    if (is.null(genoData)) {
        genoData <- base2Geno(baseData, allele.matrix)
    }
    y <- !is.na(genoData)
    idsBorder <- cumsum(colSums(y))
    idsBorder <- idsBorder[idsBorder < sum(y)]
    sum(diff(as.numeric(genoData[y]))[-idsBorder] != 0, na.rm = TRUE)
}
