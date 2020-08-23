#' Calculates the Concentration Coefficient
#'
#' @param x A vector of a income source
#' @param y The vector of total income (according to which the income source x will be sorted)
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns the numeric value of the Concentration Coefficient
#' 
#' @import pracma
#' 
#' @export
calc_concCoef <- function(x, y, w = NULL){
        
        info <- concentration_info(x, y, w)
        
        area = pracma::trapz(info$pcum_w, info$pcum_x)
        
        C = 2*(0.5 - area)
        
        C  
}

#' Makes a Concentration Curve Function
#'
#' @param x A vector of a income source
#' @param y The vector of total income (according to which the income source x will be sorted)
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Concentration Curve as outputs
#' 
#' @import pracma
#' 
#' @export
make_concentrationCurve <- function(x, y, w = NULL){
        
        info <- concentration_info(x, y, w)
        
        approxfun(x = info$pcum_w,
                  y = info$pcum_x,
                  method = "linear")
}