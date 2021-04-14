#' Fits parametric distribution functions of the family Generalized Beta 2 to binned data 
#'
#' @description A wrapper function for fitting distribution functions (pdf, cdf, quantile function, and Lorenz Curve) to binned data (edged and counts) using parametric interpolations using the GB2 distributions. See the package binequality.
#'
#' @param incomeTable A data.frame with 3 columns: N (the proportion or number of observations at the bin); lower (left endpoints of each bin); upper (right endpoints of each bin)
#' @param modelsToFit A character vector with the names of the distributions to be fitted: "GB2", "GG", "BETA2", "DAGUM", "SINGMAD", "LOGNO", "WEI", "GA", "LOGLOG", "PARETO2". As a default, all the 10 distributions are adjusted to the data and then, using the BIC statistic, the best one is chosen.
#' @param obs_mean The observed mean value of the binned distribution, if available.
#' @param grid_quantile A NIGrid object created by createNIGrid (package mvQuad) indicating how precise the numerical integration of the quantile function should be. The integration of the quantile function is used to generate the Lorenz Curve
#' 
#' @return Returns a list with the following components
#' \itemize{
#'   \descriptive A tibble with four columns: distribution (the name of the distribution fitted), estMean (the distribution mean), gini (the Gini Coefficient, if \code{gini = TRUE}, NA otherwise), endpoint (the right-hand endpoint of the support of the PDF)
#'   \pdf The PDF function. Takes values as input and returns probability densities.
#'   \cdf The CDF function. Takes values as input and returns cumulative probabilities.
#'   \quantile The Quantile function. Takes cumulative probabilities as input and returns quantile values as outputs.
#'   \lorenz The Lorenz Curve function. Takes cumulative probabilities as input and returns the proportion of income accumulated at the quantile.
#'}
#' 
#' @import tibble
#' @import dplyr
#' @import binequality 
#' @import mvQuad
#' @import pracma
#' 
#' @export
tableFit_gb2 <- function(incomeTable, 
                    modelsToFit = c("GB2", "GG", 
                                    "BETA2", "DAGUM", "SINGMAD", "LOGNO", 
                                    "WEI", "GA", "LOGLOG", "PARETO2"),
                    obs_mean = NULL,
                    grid_quantile = NULL){
        
        if(is.null(obs_mean)){
                obs_mean <- as.numeric(NA)
        }
        
        if(is.null(grid_quantile)){
                grid_quantile <- createNIGrid(dim=1, type="GLe", level=1000)        
        }
        
        fit <- with(incomeTable %>% mutate(ID = "1"), 
                    run_GB_family(ID = ID,
                                  hb = N,
                                  bin_min = lower,
                                  bin_max = upper,
                                  modelsToFit = modelsToFit,
                                  obs_mean = obs_mean,
                                  ID_name = 'ID'))
        
        bestAIC_gb2 <- fit$best_model$aic[, c("distribution",
                                              "estMean",
                                              "gini",
                                              "didConverge")]
        
        bestBIC_gb2 <- fit$best_model$bic[, c("distribution",
                                              "estMean",
                                              "gini",
                                              "didConverge")]
        
        best_gb2 <- bind_cols(tibble(criteria = c("aic", "bic")),
                              bind_rows(bestAIC_gb2, bestBIC_gb2))
        
        best_gb2 <- best_gb2 %>%
                filter(didConverge == TRUE)
        
        if(nrow(best_gb2) == 0){
                stop("No model fitted the data")
        }
        
        if(nrow(best_gb2) == 2){
                best_gb2 <- best_gb2 %>%
                        filter(criteria == "bic")
        }
        
        distribution_name <- best_gb2$distribution
        parameters        <- fit$params[[best_gb2$distribution]]
        parameters        <- as.numeric(parameters)
        
        #cf. von Hippel, Scarpino, Holas (2015). "ROBUST ESTIMATION OF INEQUALITY 
        # FROM BINNED INCOMES", p. 11
        if(best_gb2$distribution  == "DAGUM"){
                distribution_name <- "GB2"
                
                parameter_tmp = parameters
                
                parameters[1]  <- parameter_tmp[1]
                parameters[2]  <- parameter_tmp[2]
                parameters[3]  <- parameter_tmp[3]
                parameters[4]  <- 1
        }
        
        if(best_gb2$distribution  == "SINGMAD"){
                distribution_name <- "GB2"
                
                parameter_tmp = parameters
                
                parameters[1]  <- parameter_tmp[1]
                parameters[2]  <- parameter_tmp[2]
                parameters[3]  <- 1
                parameters[4]  <- parameter_tmp[3]
        }
        
        if(best_gb2$distribution  == "BETA2"){
                distribution_name <- "GB2"
                
                parameter_tmp = parameters
                
                parameters[1]  <- parameter_tmp[1]
                parameters[2]  <- 1
                parameters[3]  <- parameter_tmp[2]
                parameters[4]  <- parameter_tmp[3]
                
        }
        
        
        if(best_gb2$distribution  == "LOGLOG"){
                distribution_name <- "GB2"
                
                parameter_tmp = parameters
                
                parameters[1]  <- parameter_tmp[1]
                parameters[2]  <- parameter_tmp[2]
                parameters[3]  <- 1
                parameters[4]  <- 1
                
        }
        
        if(best_gb2$distribution  == "PARETO2"){
                distribution_name <- "GB2"
                
                parameter_tmp = parameters
                
                parameters[1]  <- parameter_tmp[1]
                parameters[2]  <- 1
                parameters[3]  <- 1
                parameters[4]  <- parameter_tmp[2]
        }
        
        
        names(parameters) <- c("mu", "sigma", "nu", "tau")
        
        par_names  <- names(parameters)[!is.na(parameters)]
        parameters <- parameters[par_names]
        
        
        args   = paste(paste(par_names, parameters, sep = " = "), collapse = ", ")
        d_name = paste0("d", distribution_name)
        p_name = paste0("p", distribution_name)
        q_name = paste0("q", distribution_name)
        
        pdf = function(x) eval(parse(text = paste0(d_name, "(x, ", args, ")")))
        
        cdf = function(x) eval(parse(text = paste0(p_name, "(x, ", args, ")")))
        
        quantile = function(x) eval(parse(text = paste0(q_name, "(x, ", args, ")")))
        
        lorenz  = makeLorenz_fromQuantile(quantile, grid_quantile)
        
        
        list(descriptive = best_gb2,
             pdf = pdf,
             cdf = cdf,
             quantile = quantile,
             lorenz   = lorenz)
        
}


#' Fits non-parametric distribution functions to binned data 
#'
#' @description A wrapper function for fitting distribution functions (pdf, cdf, quantile function, and Lorenz Curve) to binned data (edged and counts) using non-parametric interpolation strategies from th binsmooth package.
#'
#' @param incomeTable A data.frame with 3 columns: N (the proportion or number of observations at the bin); lower (left endpoints of each bin); upper (right endpoints of each bin)
#' @param method A character informing the method of non-parametric interpolation: rsub (recursive subdivision of bins), spline (cublic spline interpolation of the CDF), step  (linear interpolation of the CDF)
#' @param grid_quantile A NIGrid object created by createNIGrid (package mvQuad) indicating how precise the numerical integration of the quantile function should be. The integration of the quantile function is used to generate the Lorenz Curve
#' @param gini A logical indicating if the Gini Coefficient must be returned. Default to FALSE.
#' @param lowerBound_integration A numerical value between 0 and 1 indicating the lower bound of integration of the Quantile Function. It if used to calculate the mean and the Gini Coefficient. Integration using the full interval from 0 to 1 bears imprecise results. Default to 0.004.
#' @param upperBound_integration A numerical value between 0 and 1 indicating the upper bound of integration of the Quantile Function. It if used to calculate the mean and the Gini Coefficient. Integration using the full interval from 0 to 1 bears imprecise results. Default to 0.997.
#' @param ... Other arguments passed to 'stepbins', 'rsubbins', or 'splinebins'. 
#' 
#' @return Returns a list with the following components
#' \itemize{
#'   \descriptive A tibble with four columns: distribution (the name of the interpolation method), estMean (the distribution mean), gini (the Gini Coefficient, if \code{gini = TRUE}, NA otherwise), endpoint (the right-hand endpoint of the support of the PDF)
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
tableFit_nonParametric = function(incomeTable, 
                             method = c("rsub", "spline", "step"),
                             ..., grid_quantile = NULL, gini = F,
                             lowerBound_integration = .004, 
                             upperBound_integration = .997) {
        
        method = method[1]
        
        if(is.null(grid_quantile)){
                grid_quantile <- mvQuad::createNIGrid(dim=1, type="GLe", level=1000)
        }
        
        if(method != "spline"){
                f = get(paste0(method,"bins"), envir = as.environment("package:binsmooth"))
        }else{
                f = inequalityTools:::splinebins2
        }
        
        cat("Fitting the non-parametric function...\n")
        results <- f(bEdges  = incomeTable$upper, bCounts = incomeTable$N, ...)
        
        endPoint <- results$E
        
        CDF_order      = grep(x = names(results), pattern = "[^Inv]CDF", ignore.case = T)
        quantile_order = grep(x = names(results), pattern = "InvCDF", ignore.case = T)
        PDF_order      = grep(x = names(results), pattern = "PDF", ignore.case = T)
        
        cdf = results[[CDF_order]]
        pdf = results[[PDF_order]]
        
        cat("Getting the Quantile Function...\n")
        if(length(quantile_order) == 0){
                quantile <- inequalityTools:::inverse(results[[CDF_order]], 
                                                      lower = 0, 
                                                      upper = endPoint,
                                                      extendInt = "yes") %>% Vectorize()       
        }else{
                quantile <- results[[quantile_order]]
        }
        
        
        cat("Getting the Lorenz Function...\n")
        lorenz <- inequalityTools::make_lorenz_fromQuantile(quantile, grid_quantile)
        
        cat("Estimating the mean...\n")
        mu     <- pracma::integral(quantile, lowerBound_integration, upperBound_integration)
        
        cat("Estimating the Gini Coefficient...\n")
        if(gini == T){
                gini   <- pracma::integral(function(x) { (1/mu)*(cdf(x)*(1 - cdf(x))) }, 
                                   quantile(lowerBound_integration), 
                                   quantile(upperBound_integration))
        }else{
                gini = NA_real_
        }
        
        cat("Gathering results...\n")
        list(descriptive = tibble(distribution = paste0(method, "bins"),
                                  estMean      = mu,
                                  gini         = gini,
                                  endPoint     = endPoint),
             
             pdf      = pdf,
             cdf      = cdf,
             quantile = quantile,
             lorenz   = lorenz)
}

