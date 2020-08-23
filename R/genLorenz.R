#' Makes a Generalized Lorenz Curve Function out of a income vector
#'
#' @param x A vector of incomes 
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Generalized Lorenz Curve as outputs
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export
make_genLorenz <- function(x, w = NULL){
        
        if(is.null(w)){
                w = rep(1, length(x))
        }
        
        mu = wtd.mean(x, w)
        
        t = tibble(x, w) %>%
                filter(complete.cases(.)) %>%
                group_by(x) %>%
                summarise(w = sum(w)) %>%
                arrange(x) %>%
                mutate(cum_x = mu*(cumsum(x*w)/sum(x*w)),
                       cum_w = cumsum(w)/sum(w)) %>%
                select(cum_x, cum_w)
        
        t = bind_rows(tibble(cum_x = 0,
                             cum_w = 0),
                      t)
        
        cum_x = t$cum_x
        cum_w = t$cum_w
        
        f = approxfun(x = cum_w, y = cum_x, method = "linear")
        
}

#' Makes a Generalized Lorenz Curve Function out of a quantile function
#'
#' @param quantile_func A quantile function
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Generalized Lorenz Curve as outputs
#' 
#' @import mvQuad
#' 
#' @export
make_genLorenz_fromQuantile <- function(quantile_func, gridIntegration = NULL){
        
        if(is.null(gridIntegration)){
                gridIntegration <- mvQuad::createNIGrid(dim=1, type="GLe", level=1000)
        }
        
        genLorenz_i = function(p){
                
                mvQuad::rescale(gridIntegration, domain = c(0, p))
                mvQuad::quadrature(function(p) quantile_func(p), gridIntegration)
        }
        
        genLorenz = Vectorize(genLorenz_i)
        genLorenz
}
