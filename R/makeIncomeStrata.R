#' Cuts a continuous income variable into equal-sized bins
#'
#' @param x A vector of incomes 
#' @param n_strata A numeric value indicating how many bins will be created 
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns a numeric vector with indicating which income strata each observation belongs to.
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export
makeIncomeStrata <- function(x, n_strata = 10, w = NULL){
        
        obs_order = 1:length(x)
        
        if(is.null(w)){
                w = rep(1, length(x))
        }
        
        data <- tibble(obs_order, x, w) %>%
                filter(complete.cases(.)) %>%
                arrange(x) %>%
                mutate(p_cum  = cumsum(w)/sum(w),
                       strata = cut(p_cum,
                                    breaks = seq(0, 1, 1/n_strata),
                                    labels = 1:n_strata,
                                    include.lowest = T,
                                    ordered_result = T)) %>%
                arrange(obs_order)
        
        data <- left_join(tibble(obs_order),
                          data %>% select(obs_order, strata))
        
        as.numeric(data$strata)
}