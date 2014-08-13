`correctFUNHMM` <-
function (x, base.position = as.numeric(sub("[^0-9]*([0-9]*)[^0-9]*", 
    "\\1", names(x))), hmmFUN = hmm.vitFUN.rils, geno.probability = c(0.495, 
    0.495, 0.01), transitionFUN = phy2get.haldane.rils, emissionFUN = makeEmissionFUN(errorRate = 0.01), 
    ...) 
{
    x.cr <- hmmFUN(geno = x + 1, position = base.position, geno.probability = geno.probability, 
        transitionFUN = transitionFUN, emissionFUN = emissionFUN, 
        ...) - 1
    x.cr[x.cr == 2] <- .MPR_hetero_
    x.cr
}
