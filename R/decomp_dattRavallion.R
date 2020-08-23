#' Makes a Lorenz Curve Function out of a income vector
#'
#' @param z A numeric income value representing the poverty line
#' @param mu0 A numeric value for the mean of the group/time 0
#' @param mu1 A numeric value for the mean of the group/time 1
#' @param lorenz0 A vector valued function which takes cumulative probabilities as input and returns points at the Lorenz Curve for group/time 0
#' @param lorenz1 A vector valued function which takes cumulative probabilities as input and returns points at the Lorenz Curve for group/time 1
#' 
#' @return Returns a 1x5 tibble with:
#'     - p0_t0: Poverty rate (P0) for group/time 0
#'     - p0_t1: Poverty rate (P0) for group/time 1
#'     - growth: The difference in poverty rates due to the Growth Effect
#'     - redistr: The difference in poverty rates due to the Redistribution Effect
#'     - povertyDiff: The total difference in poverty rates between groups/times 0 and 1 (Growth Effect + Redistribution Effect). 
#'     
#' @import tibble
#' 
#' @export
decomp_dattRavallion <- function(z, mu0, mu1, lorenz0, lorenz1){
        
        p0 = inequalityTools:::p0
        
        growth <- (1/2)*((p0(z, mu1, lorenz0) - p0(z, mu0, lorenz0)) +
                                 (p0(z, mu1, lorenz1) - p0(z, mu0, lorenz1)))
        
        redistr <- (1/2)*((p0(z, mu0, lorenz1) - p0(z, mu0, lorenz0)) +
                                  (p0(z, mu1, lorenz1) - p0(z, mu1, lorenz0)))  
        
        tibble(p0_t0 = p0(z, mu0, lorenz0),
               p0_t1 = p0(z, mu1, lorenz1),
               growth, 
               redistr, 
               povertyDiff = growth + redistr)
}