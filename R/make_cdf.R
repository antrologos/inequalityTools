#' Makes a Cumulative Distribution Function (CDF) out of a income vector
#'
#' @param x A vector of a income source
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns a function which takes a vector of incomes as inputs (p) and gives points at the CDF as outputs
#' 
#' @import dplyr
#' @import tibble
#' @import spatstat
#'   
#' @export
make_cdf <- function(x, w = NULL){
        
        if(is.null(w)){
                w = rep(1, length(x))
        }
        
        data = tibble(x, w) %>%
                filter(complete.cases(.)) %>%
                group_by(x) %>%
                summarise(w = sum(w)) %>%
                ungroup()
        
        spatstat::ewcdf(x = data$x, 
                        weights = data$w, 
                        normalise = T)
}


#' Makes a Cumulative Distribution Function (CDF) out of a quantile function
#'
#' @param qf A quantile function
#' 
#' @return Returns a function which takes a vector of incomes as inputs (p) and gives points at the CDF as outputs
#' 
#' @export
make_cdf_fromQuantile <- function(qf){
        
        cdf = inequalityTools:::inverse(f = qf, lower = 0, upper = 1, extendInt = "yes")
        Vectorize(cdf)
        
}
