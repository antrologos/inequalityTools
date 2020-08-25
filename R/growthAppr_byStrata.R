#' Computes the Growth Appropriation between two periods for equal-sized income strata/bins
#'
#' @param x0 A vector of incomes at time 0
#' @param x1 A vector of incomes at time 1
#' @param n_strata A numeric value indicating how many bins will be created 
#' @param w0 (optional) A vector of sample weights at time 0
#' @param w1 (optional) A vector of sample weights at time 1
#' @param correct_for_negativeGrowth A logical value. Correct for negative growth. Default to TRUE.
#' 
#' @return Returns a tibble with the following columns:
#'     - Strata: Number of the ordered income strata
#'     - Pop. Share: Population share at the strata
#'     - Mean 0: Mean income of the strata at time 0
#'     - Mean 1: Mean income of the strata at time 1
#'     - Growth Appropriation: Growth appropriation of the strata between times 0 and 1
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export
growthAppr_byStrata <- function(x0, x1, n_strata = 10, w0 = NULL, w1 = NULL){
        
        if(is.null(w0)){
                w0 = rep(1, length(x0))
        }  
        
        if(is.null(w1)){
                w1 = rep(1, length(x1))
        }  
        
        mu0 = wtd.mean(x0, w0)
        mu1 = wtd.mean(x1, w1)
        
        if(correct_for_negativeGrowth == T){
                diff_mu = abs(mu1 - mu0)
        }else{
                diff_mu = mu1 - mu0
        }
        
        d0 = tibble(x = x0,
                    w = w0) %>%
                filter(complete.cases(.)) %>%
                mutate(strata = makeIncomeStrata(x = x, w = w, n_strata = n_strata)) %>%
                group_by(strata) %>%
                summarise(mean_0 = wtd.mean(x, w))
        
        d1 = tibble(x = x1,
                    w = w1) %>%
                filter(complete.cases(.)) %>%
                mutate(strata = makeIncomeStrata(x = x, w = w, n_strata = n_strata)) %>%
                group_by(strata) %>%
                summarise(mean_1 = wtd.mean(x, w))
        
        d = left_join(d0, d1) %>%
                arrange(strata) %>%
                mutate(p   = 1/n_strata, 
                       gac = (p*(mean_1 - mean_0))/(diff_mu),
                       gac = gac/sum(gac)) %>%
                rename(Strata = strata,
                       `Pop. Share` = p,
                       `Mean 0` = mean_0,
                       `Mean 1` = mean_1,
                       `Growth Appropriation` = gac) %>%
                select(Strata, 
                       `Pop. Share`,
                       `Mean 0`,
                       `Mean 1`,
                       `Growth Appropriation`)
        
        d
}
