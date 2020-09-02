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
calc_giniCoef_fromQuantile <- function(qf, gridIntegration = NULL){
        
        if(is.null(gridIntegration)){
                gridIntegration <- createNIGrid(dim=1, type="GLe", level=1000)
        }
        
        lorenz = make_lorenz_fromQuantile(gridIntegration = gridIntegration)
        
        
        mvQuad::rescale(gridIntegration, cbind(0, 1))
        areaUnder = mvQuad::quadrature(lorenz, gridIntegration)
        
        G = 1 - 2*areaUnder
        
        G  
}