#' Makes a Lorenz Curve Function out of a income vector
#'
#' @param x A vector of incomes 
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Lorenz Curve as outputs
#' 
#' @import tibble
#' @import dplyr
#' @import Hmisc
#' 
#' @export
make_quantileFunc <- function(x, w = NULL){

        if(is.null(w)){
                w = rep(1, length(x))
        }

        t = tibble(x, w) %>%
                filter(complete.cases(.)) %>%
                group_by(x) %>%
                summarise(w = sum(w))

        x = t$x
        w = t$w

        function(p){
                wtd.quantile(x = x, weights = w, probs = p)
        }
}
