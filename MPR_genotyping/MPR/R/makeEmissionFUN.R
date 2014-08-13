`makeEmissionFUN` <-
function (errorRate = 0.01) 
{
    E <- log(errorRate)
    E2 <- log(1 - errorRate)
    E3 <- log(0.5)
    function(h, x, n) {
        if (h != 3) 
            return(ifelse(h == x, E2, E))
        else return(n * E3)
    }
}
