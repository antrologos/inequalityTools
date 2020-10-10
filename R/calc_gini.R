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
calc_gini <- function(x, w = NULL){
        
        info <- inequalityTools:::concentration_info(x, x, w)
        
        area = pracma::trapz(info$pcum_w, info$pcum_x)
        
        G = 2*(0.5 - area)
        
        G  
}


#' Calculates the Gini Coefficient from a quantile function
#'
#' @param quantile_func A quantile function
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns the numeric value of the Gini Coefficient
#' 
#' @import mvQuad
#' 
#' @export 
calc_gini_fromQuantile <- function(qf, gridIntegration = NULL){
        
        if(is.null(gridIntegration)){
                gridIntegration <- createNIGrid(dim=1, type="GLe", level=1000)
        }
        
        lorenz = make_lorenz_fromQuantile(gridIntegration = gridIntegration)
        
        mvQuad::rescale(gridIntegration, cbind(0, 1))
        areaUnder = mvQuad::quadrature(lorenz, gridIntegration)
        
        G = 1 - 2*areaUnder
        
        G  
}


#' Calculates the Gini Coefficient from a CDF function
#'
#' @param cdf A CDF function
#' @param boundedIntegral A logical value. Indicates if the numerical integral should be performed only over a given bounded/finite interval. In this case p_lowerBound_integration and p_upperBound_integration should be informed
#' @param p_lowerBound_integration Lower bound of the finite interval over which the integral will be calculated
#' @param p_upperBound_integration Upper bound of the finite interval over which the integral will be calculated
#' 
#' @return Returns the numeric value of the Gini Coefficient
#' 
#' @import pracma
#' 
#' @export 
calc_gini_fromCDF <- function(cdf, 
                              boundedIntegral = F,
                              p_lowerBound_integration = .004,
                              p_upperBound_integration = .997){
        
        if(boundedIntegral == T){
                qf = inequalityTools:::inverse(f = cdf, lower = 0, upper = 1e12, extendInt = "yes")
                
                mu <- pracma::integral(function(x) { 1 - cdf(x) }, 
                                       xmin = qf(p_lowerBound_integration), 
                                       xmax = qf(p_upperBound_integration))
                
                G  <- pracma::integral(function(x) { (1/mu)*(cdf(x)*(1 - cdf(x))) }, 
                                       xmin = qf(p_lowerBound_integration), 
                                       xmax = qf(p_upperBound_integration))
        }else{
                
                mu <- pracma::integral(function(x) { 1 - cdf(x) }, 
                                       xmin = 0, 
                                       xmax = Inf)
                
                G  <- pracma::integral(function(x) { (1/mu)*(cdf(x)*(1 - cdf(x))) }, 
                                       xmin = 0, 
                                       xmax = Inf)
        }
        
        G  
}



