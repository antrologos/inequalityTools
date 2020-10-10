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

#' Makes a Lorenz Curve Function using information about means and size of income strata
#'
#' @param group_size A numerical vector informing the size (i.e population of proportion) of each group
#' @param group_means A numerical vector informing the mean income of each group. Must have the same number of elements tha 'group_size'
#' 
#' @return Returns a function which takes a vector of probabilities as inputs (p) and gives points at the Lorenz Curve as outputs
#' 
#' @export
make_lorenz_fromMeansByStrata <- function(group_size, group_means,
                                          combine_function = "mean",
                                          precision = .00005){
        
        group_data = data.frame(group_size, group_means)
        group_data = group_data[order(group_data$group_means), ]
        
        group_means = group_data$group_means
        group_size   = group_data$group_size
        
        p = cumsum(c(0, group_size))/sum(group_size)
        r = cumsum(c(0, group_means*diff(p)))/sum(group_means*diff(p))
        
        betas <- matrix(NA, ncol=3, nrow = (length(p)-2))
        r2 <- NULL
        for(i in 1:(length(p)-2)){
                reg        <- lm(r[i:(i+2)] ~ p[i:(i+2)] + I((p[i:(i+2)])^2))
                r2[i]      <- summary(reg)$r.squared
                betas[i, ] <- coef(reg)
        }
        colnames(betas) = c("intercept", "p", "p^2")
        
        # Interpolando
        make_parabolafunctions <- function(i, x, beta_matrix){
                env_tmp       <- new.env(parent = emptyenv())
                env_tmp$linha <- i
                env_tmp$beta  <- beta_matrix
                
                f = function(x) env_tmp$beta[env_tmp$linha,1] + env_tmp$beta[env_tmp$linha,2]*x + env_tmp$beta[env_tmp$linha,3]*(x^2)
                f
        }
        
        functions <- list()
        for(j in 1:nrow(betas)){
                functions[[j]] <- make_parabolafunctions(i = j,
                                                         x = x,
                                                         beta_matrix = betas)
        }
        
        size = length(p)
        threasholds <- cbind(1:(size-2),3:size)
        p_threasholds <- p[threasholds]
        dim(p_threasholds) = dim(threasholds)
        
        is_between <-function(x){
                which(x >= p_threasholds[,1] & x <= p_threasholds[,2])
        }
        
        lorenz_i <- function(y){
                
                functionsNumbers <- is_between(y)
                lorenz_value_temp = NULL
                
                count <- 1
                for(functionsNumber in functionsNumbers){
                        lorenz_value_temp_test1 <- functions[[functionsNumber]](y)
                        lorenz_value_temp[count] <- lorenz_value_temp_test1
                        count = count + 1
                }
                
                expression = paste0(combine_function,"( c(",paste(lorenz_value_temp, collapse = ", "),"), na.rm = T)")
                
                lorenz_value <- eval(expr = parse(text = expression))
                
                lorenz_value
        }
        
        lorenz_parabola = Vectorize(lorenz_i)
        
        p = seq(0, 1, precision)
        l = lorenz_parabola(p)
        
        approxfun(x = p, y = l, method = "linear")
        
}
