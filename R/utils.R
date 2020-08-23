# Makes the Inverse Function out a function
inverse = function (f, lower = -100, upper = 100, extendInt = "no") {
        function (y) uniroot((function(x) f(x) - y), lower = lower, upper = upper,  extendInt = extendInt)$root
}


# P0, as a function of z, mu, and lorenz (for the Datt-Ravallion Decompostion)
#' @import pracma
p0 <- function(z, mu, lorenz){
        derv_lorenz <- function(x){ 
                x = ifelse(x < .Machine$double.eps^(1/2), .Machine$double.eps^(1/2), x)
                x = ifelse(x > 1 - .Machine$double.eps^(1/2), 1 - .Machine$double.eps^(1/2), x)
                
                pracma::grad(f = lorenz, x0 = x, heps = .Machine$double.eps^(1/2))}
        
        inv_derv_lorenz <- inequalityTools:::inverse(f = derv_lorenz, 
                                                     lower = 0 + .Machine$double.eps^(1/2), 
                                                     upper = 1 - .Machine$double.eps^(1/2),
                                                     extendInt = "yes")
        inv_derv_lorenz(z/mu)
}


# Auxiliary function for the Concentration Curve and Concentration Coefficient
#' @import tibble
#' @import dplyr 
concentration_info <- function(x, y, w = NULL){
        
        if(is.null(w)){
                w = rep(1, length(y))
        }
        
        if(!all.equal(length(x), length(y), length(w))){
                stop("x, y, and w must be the same length")
        }
        
        data <- tibble(x, y, w) %>%
                filter(complete.cases(.)) %>%
                group_by(x, y) %>%
                summarise(w = sum(w)) %>%
                ungroup() %>%
                arrange(y) %>%
                mutate(wx     = w*x,
                       cum_w  = cumsum(w),
                       cum_x  = cumsum(wx),
                       pcum_w = cum_w/last(cum_w),
                       pcum_x = cum_x/last(cum_x))
        
        bind_rows(tibble(pcum_w = 0,
                         pcum_x = 0),
                  data %>% select(pcum_w, pcum_x))
}

# Transpose a data.frame
#' @import tibble
#' @import dplyr
transpose <- function(t, hasColNames = T){
        df <- t(t)
        statistics <- rownames(df)
        
        if(hasColNames == T){
                cols      <- paste(rownames(df)[1], df[1,], sep="_")  
                df        <- as_tibble(df[-1,])
                names(df) <- cols
                df        <- bind_cols(tibble(statistics = statistics[-1]), df)
        }else{
                cols      <- paste0("col", 1:ncol(df))  
                df        <- as_tibble(df)
                names(df) <- cols
                df        <- bind_cols(tibble(statistics), df)
        }
        
        df
}

# Replaces all NA values in a vector for given specified replacement (default to 0)
replaceNAvalues = function(x, replacement = 0){
        x[is.na(x)] <- replacement
        x
}
