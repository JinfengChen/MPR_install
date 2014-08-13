`mergeBinMap` <-
function (geno.fill) 
{
    pre <- NULL
    pre.nna <- NULL
    IsUniq <- rep(1, nrow(geno.fill))
    for (j in 1:nrow(geno.fill)) {
        cur.nna <- which(!is.na(geno.fill[j, ]))
        cur <- geno.fill[j, cur.nna]
        if (identical(pre.nna, cur.nna) && identical(pre, cur)) 
            IsUniq[j] <- 0
        else {
            pre <- cur
            pre.nna <- cur.nna
        }
    }
    geno.fill[unique(c(which(IsUniq == 1)[-1] - 1, nrow(geno.fill))), 
        ]
}
