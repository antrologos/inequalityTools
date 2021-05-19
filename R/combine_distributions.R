#' Combines two different distributions that were fitted using inequalityTools 
#'
#' @description  Combines two distributions fitted with inequalityTools by weighted summing their PDFs. Then recalculates the other distribution functions (pdf, cdf, quantile function, and Lorenz Curve) and descriptives.
#'
#' @param fit1 A interpolated distribution fitted with tableFit_nonParametric or tableFit_gb2
#' @param fit2 Another interpolated distribution fitted with tableFit_nonParametric or tableFit_gb2
#' @param n1 The number of observations or proportion of the fit1.
#' @param n2 The number of observations or proportion of the fit2.
#' @param quadrature_precision Gauss-Legendre Quadrature level - argument to be passed to mvQuad::createNIGrid. Default to 5e3
#' @param points_to_adjust The number of points from which the cdf, quantile function and Lorenz Curve will be calculated. Default to 2000
#' @param upper_limit The maximum income value of the distribution support. Default to 1e5
#' @param upper_p_limit A cumulative probability value which will be taken as the maximum for the quantile function. After it, quantile function will be top coded. Default to .9995
#' 
#' @return Returns a list with the following components
#' \itemize{
#'   \descriptive A tibble with two columns: estMean (the distribution mean), gini (the Gini Coefficient, if \code{gini = TRUE}, NA otherwise)
#'   \pdf The PDF function. Takes values as input and returns probability densities.
#'   \cdf The CDF function. Takes values as input and returns cumulative probabilities.
#'   \quantile The Quantile function. Takes cumulative probabilities as input and returns quantile values as outputs.
#'   \lorenz The Lorenz Curve function. Takes cumulative probabilities as input and returns the proportion of income accumulated at the quantile.
#'}
#' 
#' @import tibble
#' @import binsmooth
#' @import mvQuad
#' @import pracma
#' 
#' @export
combine_distributions <- function(fit1, fit2, n1, n2, 
                                  quadrature_precision = 5e3,
                                  points_to_adjust = 2000,
                                  upper_limit   = 1e5,
                                  upper_p_limit = .9995){
        
        p1 = n1/(n1 + n2)
        p2 = 1 - p1
        
        pdf1 = fit1$pdf
        pdf2 = fit2$pdf
        
        pdf = function(x) p1*pdf1(x) + p2*pdf2(x)
        
        grid_tmp = mvQuad::createNIGrid(dim = 1, type = "GLe", level =  quadrature_precision)
        cdf_tmp   = Vectorize(function(x){
                mvQuad::rescale(grid_tmp, domain = matrix(c(0, x), ncol = 2))
                mvQuad::quadrature(pdf, grid_tmp)
        })
        quantile_tmp = Vectorize(inequalityTools:::inverse(f = cdf_tmp, lower = 0, upper = upper_limit, extendInt = "yes"))
        
        x = quantile_tmp(seq(0, upper_p_limit, length.out = points_to_adjust))
        cdf_points = sapply(x, cdf_tmp)
        
        selected   = c(TRUE, !(diff(cdf_points) < 0))
        cdf_points <- cdf_points[selected]
        x_points   <- x[selected]
        
        cdf      = approxfun(x = x_points,   y = cdf_points, yleft = NA, yright = 1,  method = "linear")
        quantile = approxfun(x = cdf_points, y = x_points,   yleft = NA, yright = max(x), method = "linear")
        lorenz   = inequalityTools::make_lorenz_fromQuantile(quantile)
        
        descriptive = tibble(estMean      = integrate(function(x) x*pdf(x), 0, Inf, subdivisions = 5000, rel.tol = .0005)[[1]],
                             gini         = 1 - 2*integrate(lorenz, 0, 1, subdivisions = 5000, rel.tol = .0005)[[1]])
        
        list(descriptive = descriptive,
             pdf      = pdf,
             cdf      = cdf,
             quantile = quantile,
             lorenz   = lorenz)
}


