#' Computes the mean income for a strata defined by a lower bound and an upper bound in the ordered cumulative income distribution
#'
#' @param x A vector of incomes 
#' @param lowerBound A value between 0 and 1, indicating the lower bound of the strata
#' @param upperBound A value between 0 and 1, indicating the upper bound of the strata
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns the mean (a numeric value) within the strata
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export
meanBystrata <- function(x, lowerBound = 0, upperBound = 1, w = NULL){
        
        if(is.null(w)){
                w = rep(1, length(x))
        }
        
        data <- tibble(obs_order, x, w) %>%
                filter(complete.cases(.)) %>%
                arrange(x) %>%
                mutate(p_cum = cumsum(w)/sum(w))
        
        if(lowerBound == 0){
                subData <- data %>%
                        filter(p_cum >= lowerBound, p_cum <= upperBound)
        }else{
                subData <- data %>%
                        filter(p_cum > lowerBound, p_cum <= upperBound)
        }
        
        sum(subData$x*subData$w)/sum(subData$w)
}