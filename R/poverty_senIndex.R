#' Calculates Amartya Sen's Index of poverty
#'
#' @param x A vector of a income source
#' @param z A numeric income value representing the poverty line
#' @param w (optional) A vector of sample weights
#' @param alpha An arbitrary inequality aversion parameter: alpha = 0 returns the poverty rate; alpha = 1 returns the average income gap. Default to 0.
#' 
#' @return Returns the numeric value of the FGT measure of poverty
#' 
#' @import tibble
#' @import dplyr
#' 
#' @references 
#' Sen, Amartya. (1976) "Poverty: An ordinal approach to measurement." 
#'     Econometrica: Journal of the Econometric Society: 219-231.
#' 
#' @export 
poverty_senIndex <- function(x, z, w = NULL){
        H <- fgt(x = x, w = w, z = z, alpha = 0)
        I <- fgt(x = x, w = w, z = z, alpha = 1)
        
        if(length(z) == 1){
                z = rep(z, length(x))
        }
        
        data <- tibble(x,w,z) %>% 
                filter(complete.cases(.)) %>%
                filter(x < z)
        
        rm(x, w, z)
        
        G <- calcSGini(x = data$x, w = data$w)$ineq$index[[1]]
        q <- sum(data$w)
        
        # Simplified formula (for large n)
        # H*(I + (1 - I)*G)
        
        # Complete formula
        H*(1 - (1 - I)*(1 - G*(q/(q+1))))
}