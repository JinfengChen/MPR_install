`geno2Cross` <-
function (geno.data, pheno.data, cm_unit = 250000) 
{
    myGeno.data <- geno.data
    geno.value <- sort(unique(na.omit(c(myGeno.data))))
    if (identical(as.numeric(geno.value), c(0, 1))) {
        myGeno.data <- myGeno.data + 1
        geno.value <- sort(unique(na.omit(c(myGeno.data))))
    }
    if (!identical(as.numeric(geno.value), c(1, 2))) 
        stop("the value of geno.data should be 1, 2, or NA.")
    ids.RILs <- match(colnames(myGeno.data), rownames(pheno.data))
    myGeno.site <- data.frame(Chr = as.numeric(substr(rownames(myGeno.data), 
        1, 2)), Position = as.numeric(substr(rownames(myGeno.data), 
        3, 10)))
    myCrossData <- list()
    myCrossData$geno <- lapply(split(1:nrow(myGeno.data), myGeno.site$Chr), 
        function(ids) {
            myMarkers <- myGeno.site$Position[ids]/cm_unit
            names(myMarkers) <- rownames(myGeno.data)[ids]
            myMap <- list(data = t(myGeno.data[ids, ]), map = myMarkers)
            class(myMap) <- "A"
            myMap
        })
    myCrossData$pheno <- as.data.frame(pheno.data[ids.RILs, ])
    class(myCrossData) <- c("riself", "cross")
    attr(myCrossData, "alleles") <- c("A", "B")
    myCrossData
}
