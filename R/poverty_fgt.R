#' Calculates the Foster-Greer-Thorbecke (FGT) family of poverty measures
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
#' @export
poverty_fgt <- function(x, z, w = NULL, alpha = 0){
        if(is.null(w)){
                w = rep(1, length(x))
        }
        
        if(length(z) == 1){
                z = rep(z, length(x))
        }
        
        data <- tibble(x,w,z) %>% filter(complete.cases(.))
        
        rm(x, w, z)
        
        data <- data %>%
                mutate(g   = ifelse(x < z, ((z - x)/z), 0),
                       fgt = ifelse(x < z, g^alpha, 0)) %>%
                summarise(fgt = sum(w*fgt),
                          n   = sum(w))
        
        data$fgt/data$n
}