#' Makes a Growth Appropriation Curve Function out of two income vectors
#'
#' @param x0 A vector of incomes for group/time 0
#' @param x1 A vector of incomes for group/time 1 
#' @param w0 (optional) A vector of sample weights for group/time 0
#' @param w1 (optional) A vector of sample weights for group/time 1
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Growth Appropriation Curve as outputs
#' 
#' @import mvQuad
#' 
#' @export
make_growthApprCurve = function(x0, x1, w0 = NULL, w1 = NULL, 
                                gridIntegration = NULL){
        
        if(is.null(w0)){
                w0 = rep(1, length(x0))
        }
        
        if(is.null(w1)){
                w1 = rep(1, length(x1))
        }
        
        if(is.null(gridIntegration)){
                gridIntegration = mvQuad::createNIGrid(dim = 1, type = "GLe", level = 2000)
        }
        
        qf0 = make_quantileFunc(x = x0, w = w0)
        qf1 = make_quantileFunc(x = x1, w = w1)
        
        mu0 = mvQuad::quadrature(qf0, gridIntegration)
        mu1 = mvQuad::quadrature(qf1, gridIntegration)
        
        function(p){
                (qf1(p) - qf0(p))/(mu1 - mu0)
        }
}


#' Makes a Growth Appropriation Curve Function out of two quantile functions
#'
#' @param qf0 A quantile function for group/time 0
#' @param qf1 A quantile function for group/time 1 
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Growth Appropriation Curve as outputs
#' 
#' @import mvQuad
#' 
#' @export
make_growthApprCurve_fromQuantile <- function(qf0, qf1, 
                                              gridIntegration = NULL){
        
        if(is.null(gridIntegration)){
                gridIntegration = mvQuad::createNIGrid(dim = 1, type = "GLe", level = 2000)
        }
        
        mu0 = mvQuad::quadrature(qf0, gridIntegration)
        mu1 = mvQuad::quadrature(qf1, gridIntegration)
        
        function(p){
                (qf1(p) - qf0(p))/(mu1 - mu0)
        }
}


