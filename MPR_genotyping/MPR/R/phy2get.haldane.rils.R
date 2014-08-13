`phy2get.haldane.rils` <-
function (a, b, l) 
{
    d <- l/(._rice_phy2get_factor_ * 100)
    p <- (1 - exp(-2 * d))
    p <- p/(1 + p)
    ifelse(a == b, 1 - p, p)
}
