#' Computes the Concentration Coefficient for the income growth between two periods
#'
#' @param x0 A vector of incomes at time 0
#' @param x1 A vector of incomes at time 1
#' @param w0 (optional) A vector of sample weights at time 0
#' @param w1 (optional) A vector of sample weights at time 1
#' 
#' @return Returns the Growth Concentration Coenfficient (numeric value)
#' 
#' @export
calc_growthConcCoef = function(x0, x1, w0 = NULL, w1 = NULL){
        
        if(is.null(w0)){
                w0 = rep(1, length(x0))
        }  
        
        if(is.null(w1)){
                w1 = rep(1, length(x1))
        }  
        
        mu0 = wtd.mean(x0, w0)
        mu1 = wtd.mean(x1, w1)
        
        g0 = calc_giniCoef(x0, w0)
        g1 = calc_giniCoef(x1, w1)
        
        1 - (mu1*(1 - g1) - mu0*(1 - g0))/(mu1 - mu0)        
}

