`genotypeCallsBayes` <-
function (ALLELE.num, errorRate = 5e-04, eps = 1e-10, maxIterate = 100, 
    verbose = FALSE) 
{
    P0.n <- 0.5
    P0.p1 <- P0.p2 <- (1 - P0.n)/2
    ALLELE.prob <- NULL
    P.n <- 1/3
    P.p1 <- P.p2 <- (1 - P.n)/2
    numIterate <- 1
    E <- errorRate
    while (abs(P0.n - P.n) > eps && numIterate <= maxIterate) {
        P0.p1 <- P.p1
        P0.p2 <- P.p2
        P0.n <- P.n
        n <- rowSums(ALLELE.num)
        k <- ALLELE.num[, 1]
        ALLELE.prob <- cbind((1 - E)^k * E^(n - k), 0.5^n, (1 - 
            E)^(n - k) * E^k) %*% diag(c(P0.p1, P0.n, P0.p2))
        ALLELE.type <- rep(NA, nrow(ALLELE.prob))
        ALLELE.maxprob <- rowMax(ALLELE.prob)
        ALLELE.type[ALLELE.prob[, 1] == ALLELE.maxprob] <- 1
        ALLELE.type[ALLELE.prob[, 3] == ALLELE.maxprob] <- 3
        ALLELE.type[ALLELE.prob[, 2] == ALLELE.maxprob | ALLELE.prob[, 
            1] == ALLELE.prob[, 3]] <- 2
        P.n <- mean(ALLELE.type == 2)
        P.p1 <- P.p2 <- (1 - P.n)/2
        if (verbose) 
            cat(P.p1, P.n, P.p2, table(filter.markers <- ALLELE.type != 
                2), "\n", sep = "\t")
        numIterate <- numIterate + 1
    }
    ALLELE.prob <- ALLELE.prob/rowSums(ALLELE.prob, na.rm = TRUE)
    list(prop = c(P.p1, P.n, P.p2), prob = ALLELE.prob, type = ALLELE.type)
}
