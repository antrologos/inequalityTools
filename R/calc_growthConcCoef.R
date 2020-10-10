#' Computes the Concentration Coefficient for the income growth between two periods
#'
#' @param x0 A vector of incomes at time 0
#' @param x1 A vector of incomes at time 1
#' @param w0 (optional) A vector of sample weights at time 0
#' @param w1 (optional) A vector of sample weights at time 1
#' @param correct_for_negativeAppropriation A logical value. If \code{TRUE}, applies a alternative method that corrects for negative values in the Growth Concentration Curve. Default to \code{FALSE}.
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns the Growth Concentration Coefficient (numeric value)
#' 
#' @export
calc_growthConcCoef = function(x0, 
                               x1, 
                               w0 = NULL, 
                               w1 = NULL,
                               correct_for_negativeAppropriation = F,
                               gridIntegration = NULL
                               ){
        
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
        
        cc = 1 - (mu1*(1 - g1) - mu0*(1 - g0))/(mu1 - mu0)
        
        if(correct_for_negativeAppropriation == T){
                
                gcc = make_growthConcCurve(x0 = x0,
                                           x1 = x1,
                                           w0 = w0,
                                           w1 = w1,
                                           gridIntegration = gridIntegration)
                
                # 1
                f1_AreaBelowCurveAboveP_CurveAboveP = function(p) { 
                        ifelse((gcc(p) > p), -(gcc(p) - p), 0)
                }
                
                # 2
                f2_AreaAboveCurve_CurveAbove0BelowP = function(p) {
                        ifelse(gcc(p) <= p, ifelse(gcc(p) >= 0, p - gcc(p), p), 0)
                }
                
                # 3
                f3_AreaAboveCurve_CurveBelow0 = function(p) {
                        ifelse(gcc(p) < 0, -gcc(p), 0)
                }
                
                # 4
                f4_AreaBetween0andCurveOrP_nonNegativeCurve = function(p) {
                        ifelse(gcc(p) >= 0, ifelse(gcc(p) > p, p, gcc(p)), 0)
                }
                
                area_1 = pracma::integral(f1_AreaBelowCurveAboveP_CurveAboveP,         0, 1, no_intervals = 5000)
                area_2 = pracma::integral(f2_AreaAboveCurve_CurveAbove0BelowP,         0, 1, no_intervals = 5000)
                area_3 = pracma::integral(f3_AreaAboveCurve_CurveBelow0,               0, 1, no_intervals = 5000)
                area_4 = pracma::integral(f4_AreaBetween0andCurveOrP_nonNegativeCurve, 0, 1, no_intervals = 5000)
                
                corrected_cc = (area_1 + area_2 + area_3)/(area_1 + area_2 + area_3 + area_4)
                
                if(sign(cc) != sign(corrected_cc)){
                        corrected_cc = corrected_cc*sign(cc)
                }
                
                cc = corrected_cc
        }
        
        cc
}


#' Computes the Concentration Coefficient for the income growth between two periods using Quantile Functions as input
#'
#' @param qf0 A quantile function of the income distribution at time 0
#' @param qf1 A quantile function of the income distribution at time 1
#' @param correct_for_negativeAppropriation A logical value. If \code{TRUE}, applies a alternative method that corrects for negative values in the Growth Concentration Curve. Default to \code{FALSE}.
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns the Growth Concentration Coefficient (numeric value)
#' 
#' @export
calc_growthConcCoef_fromQuantile = function(qf0, qf1,
                                            correct_for_negativeAppropriation = F,
                                            gridIntegration = NULL){
        
        if(is.null(gridIntegration)){
                gridIntegration = mvQuad::createNIGrid(dim = 1, type = "GLe", level = 2000)
        }
        
        mu0 = pracma::integral(qf0, xmin = 0, xmax = 1, no_intervals = 5000)
        mu1 = pracma::integral(qf1, xmin = 0, xmax = 1, no_intervals = 5000)
        
        g0 = calc_gini_fromQuantile(qf0, gridIntegration = gridIntegration)
        g1 = calc_gini_fromQuantile(qf1, gridIntegration = gridIntegration)
        
        cc = 1 - (mu1*(1 - g1) - mu0*(1 - g0))/(mu1 - mu0)
        
        if(correct_for_negativeAppropriation == T){
                
                gcc = make_growthConcCurve_fromQuantile(qf0 = qf0,
                                                        qf1 = qf1,
                                                        gridIntegration = gridIntegration)
                
                # 1
                f1_AreaBelowCurveAboveP_CurveAboveP = function(p) { 
                        ifelse((gcc(p) > p), -(gcc(p) - p), 0)
                }
                
                # 2
                f2_AreaAboveCurve_CurveAbove0BelowP = function(p) {
                        ifelse(gcc(p) <= p, ifelse(gcc(p) >= 0, p - gcc(p), p), 0)
                }
                
                # 3
                f3_AreaAboveCurve_CurveBelow0 = function(p) {
                        ifelse(gcc(p) < 0, -gcc(p), 0)
                }
                
                # 4
                f4_AreaBetween0andCurveOrP_nonNegativeCurve = function(p) {
                        ifelse(gcc(p) >= 0, ifelse(gcc(p) > p, p, gcc(p)), 0)
                }
                
                area_1 = pracma::integral(f1_AreaBelowCurveAboveP_CurveAboveP,         0, 1, no_intervals = 5000)
                area_2 = pracma::integral(f2_AreaAboveCurve_CurveAbove0BelowP,         0, 1, no_intervals = 5000)
                area_3 = pracma::integral(f3_AreaAboveCurve_CurveBelow0,               0, 1, no_intervals = 5000)
                area_4 = pracma::integral(f4_AreaBetween0andCurveOrP_nonNegativeCurve, 0, 1, no_intervals = 5000)
                
                corrected_cc = (area_1 + area_2 + area_3)/(area_1 + area_2 + area_3 + area_4)
                
                if(sign(cc) != sign(corrected_cc)){
                        corrected_cc = corrected_cc*sign(cc)
                }
                
                cc = corrected_cc
        }
        
        cc
}