`.loopMPR` <-
function (genoData, allele.matrix, numStep = 1, ...) 
{
    rawGenoSum <- NumRecomEvents(genoData = genoData)
    for (x in 1:(nrow(allele.matrix) - numStep + 1)) {
        ids <- x:(x + numStep - 1)
        new.genoData <- genoData
        new.genoData[ids, ] <- abs(new.genoData[ids, ] - 1)
        newGenoSum <- NumRecomEvents(genoData = new.genoData)
        if (newGenoSum < rawGenoSum) {
            allele.matrix[ids, ] <- allele.matrix[ids, c(2, 1)]
            genoData <- new.genoData
            rawGenoSum <- newGenoSum
        }
    }
    list(rawGenoSum, allele.matrix, genoData)
}
`.MPR_hetero_` <-
0.5
`.MPR_p0_` <-
0
`.MPR_p1_` <-
1
`._rice_phy2get_factor_` <-
244000
