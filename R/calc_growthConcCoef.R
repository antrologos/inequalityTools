#' Computes the Concentration Coefficient for the income growth between two periods
#'
#' @param x0 A vector of incomes at time 0
#' @param x1 A vector of incomes at time 1
#' @param w0 (optional) A vector of sample weights at time 0
#' @param w1 (optional) A vector of sample weights at time 1
#' @param correct_for_negativeAppropriation A logical value. Correct for negative values in the Growth Appropriation Curve. Default to TRUE.
#' @param correct_for_negativeGrowth A logical value. Correct for negative growth. Default to TRUE.
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns the Growth Concentration Coenfficient (numeric value)
#' 
#' @export
calc_growthConcCoef = function(x0, 
                               x1, 
                               w0 = NULL, 
                               w1 = NULL,
                               
                               correct_for_negativeAppropriation = T,
                               correct_for_negativeGrowth = T,
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
                
                gcc_aboveCurve_above1 = function(p) {
                        ifelse(gcc(p) > 1, -(gcc(p) - 1), 0)
                }
                
                gcc_aboveCurve_below0 = function(p) {
                        ifelse(gcc(p) < 0, -gcc(p), 0)
                }
                
                gcc_aboveCurve_above0belowP = function(p) {
                        ifelse((gcc(p) <= p) & (gcc(p) >= 0), p - gcc(p), 0)
                }
                
                gcc_belowCurve_above0belowP = function(p) {
                        ifelse((gcc(p) <= p) & (gcc(p) >= 0), gcc(p), 0)
                }
                
                gcc_belowCurve_abovePbelow1 = function(p) {
                        ifelse((gcc(p) >= p) & (gcc(p) <= 1), -(gcc(p) - p), 0)
                }
                
                area_a1 = pracma::integral(gcc_aboveCurve_above1,       0, 1, no_intervals = 5000)
                area_a2 = pracma::integral(gcc_aboveCurve_below0,       0, 1, no_intervals = 5000)
                area_a3 = pracma::integral(gcc_aboveCurve_above0belowP, 0, 1, no_intervals = 5000)
                area_b1 = pracma::integral(gcc_belowCurve_above0belowP, 0, 1, no_intervals = 5000)
                area_b2 = pracma::integral(gcc_belowCurve_abovePbelow1, 0, 1, no_intervals = 5000)
                
                corrected_cc = (area_a1 + area_a2 + area_a3)/(area_a1 + area_a2 + area_a3 + area_b1 + area_b1 + area_b2)
                
                if(sign(cc) != sign(corrected_cc)){
                        corrected_cc = corrected_cc*sign(cc)
                }
                
                cc = corrected_cc
        }
        
        if(correct_for_negativeGrowth == T){
                return(sign(mu1 - mu0)*cc)
        }else{
                return(cc)
        }
}
