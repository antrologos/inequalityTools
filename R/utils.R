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

# Makes the Concentration Coefficient to vary between -1 and 1
normalize_neg1_pos1 = function(x){
        
        2* ( (exp(x)/(exp(x)+1)) -.5)
        
}

# splinebins with extra parameters
splinebins2 <- function(bEdges, 
                        bCounts, 
                        m = NULL, 
                        numIterations = 16, 
                        openBracketInitialMultiplier = 2,
                        tailEnd_multiplier = 1.05,
                        length_cdf = 2000,
                        monoMethod = c("hyman", "monoH.FC")){
        
        monoMethod <- match.arg(monoMethod)
        L <- length(bCounts)
        if (!(is.na(bEdges[L]) | is.infinite(bEdges[L]))) 
                warning("Top bin is bounded. Expect inaccurate results.\n")
        tot <- sum(bCounts)
        if (is.null(m)) {
                warning("No mean provided: expect inaccurate results.\n")
                m <- sum(0.5 * (c(bEdges[1:(L - 1)], openBracketInitialMultiplier * bEdges[L - 1]) + 
                                        c(0, bEdges[1:(L - 1)])) * bCounts/tot)
        }
        sumbAreas <- vapply(1:L, function(i) {
                sum(bCounts[1:i])/tot
        }, numeric(1))
        
        tailEnd <- tailEnd_multiplier * bEdges[L - 1]
        
        x <- c(0, bEdges[1:(L - 1)], tailEnd, tailEnd * 1.01)
        y <- c(0, sumbAreas, 1)
        f <- splinefun(x, y, method = monoMethod)
        est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
        
        shrinkFactor <- 1
        shrinkMultiplier <- 0.995
        while (est_mean > m) {
                x <- x * shrinkMultiplier
                tailEnd <- tailEnd * shrinkMultiplier
                f <- splinefun(x, y, method = monoMethod)
                est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
                shrinkFactor <- shrinkFactor * shrinkMultiplier
        }
        
        if (bCounts[L] > 0) {
                while (est_mean < m) {
                        tailEnd <- tailEnd * 2
                        x[L + 1] <- tailEnd
                        x[L + 2] <- tailEnd * 1.01
                        f <- splinefun(x, y, method = monoMethod)
                        est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
                }
                l <- bEdges[L - 1]
                r <- tailEnd
                for (i in 1:numIterations) {
                        tailEnd <- (l + r)/2
                        x[L + 1] <- tailEnd
                        x[L + 2] <- tailEnd * 1.01
                        f <- splinefun(x, y, method = monoMethod)
                        est_mean <- tailEnd - pracma::integral(f, 0, tailEnd)
                        if (est_mean > m) 
                                r <- tailEnd
                        else l <- tailEnd
                }
        }
        
        xfix <- seq(0, tailEnd, length.out = length_cdf)
        yfix <- f(xfix)
        finv <- splinefun(yfix, xfix, method = monoMethod)
        
        splineCDF <- function(x) {
                ifelse(x < 0, 0, ifelse(x > tailEnd, 1, f(x)))
        }
        
        splinePDF <- function(x) {
                ifelse(x < 0 | x > tailEnd, 0, f(x, deriv = 1))
        }
        
        splineInvCDF <- function(x) {
                ifelse(x < 0, 0, ifelse(x > 1, tailEnd, finv(x)))
        }
        
        median_bin <- which(cumsum(bCounts) > sum(bCounts)/2)[1]
        med_fit <- splineInvCDF(0.5)
        
        if(median_bin > 1){
                fitWarn <- (med_fit < bEdges[median_bin - 1]) | (med_fit > bEdges[median_bin])
        }else{
                fitWarn <- (med_fit > bEdges[median_bin])
        }
        
        if (fitWarn) 
                warning("Inaccurate fit detected. Proceed with caution.\n")
        
        return(list(splinePDF = splinePDF, splineCDF = splineCDF, 
                    E = tailEnd, est_mean = est_mean, shrinkFactor = shrinkFactor, 
                    splineInvCDF = splineInvCDF, fitWarn = fitWarn))
}