#' Calculates the Gini Coefficient
#'
#' @param x A vector of a income source
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns the numeric value of the Gini Coefficient
#' 
#' @import pracma
#' 
#' @export
calc_giniCoef <- function(x, w = NULL){
        
        info <- inequalityTools:::concentration_info(x, x, w)
        
        area = pracma::trapz(info$pcum_w, info$pcum_x)
        
        G = 2*(0.5 - area)
        
        G  
}