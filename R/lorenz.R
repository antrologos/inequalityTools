#' Makes a Lorenz Curve Function out of a income vector
#'
#' @param x A vector of incomes 
#' @param w (optional) A vector of sample weights
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Lorenz Curve as outputs
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export
make_lorenz <- function(x, w = NULL){

        if(is.null(w)){
                w = rep(1, length(x))
        }

        t = tibble(x, w) %>%
                filter(complete.cases(.)) %>%
                group_by(x) %>%
                summarise(w = sum(w)) %>%
                arrange(x) %>%
                mutate(cum_x = cumsum(x*w)/sum(x*w),
                       cum_w = cumsum(w)/sum(w)) %>%
                dplyr::select(cum_x, cum_w)

        t = bind_rows(tibble(cum_x = 0,
                             cum_w = 0),
                      t)

        cum_x = t$cum_x
        cum_w = t$cum_w

        f = approxfun(x = cum_w, y = cum_x, method = "linear")

}
#' Makes a Lorenz Curve Function out of a quantile function
#'
#' @param quantile_func A quantile function
#' @param gridIntegration (optional) A grid of class 'NIGrid' for multivariate numerical integration (mvQuad package)
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Lorenz Curve as outputs
#' 
#' @import mvQuad
#' 
#' @export
make_lorenz_fromQuantile <- function(qf, gridIntegration = NULL){

        if(is.null(gridIntegration)){
                gridIntegration <- createNIGrid(dim=1, type="GLe", level=1000)
        }

        rescale(gridIntegration, domain = c(0 + .Machine$double.eps, 1 - .Machine$double.eps))
        mu = quadrature(qf, gridIntegration)

        lorenz_i = function(p){

                mvQuad::rescale(gridIntegration, domain = c(0, p))
                mvQuad::quadrature(function(p) (1/mu)*qf(p), gridIntegration)
        }

        lorenz = Vectorize(lorenz_i)
        lorenz
}
