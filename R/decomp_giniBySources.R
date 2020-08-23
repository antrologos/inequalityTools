#' Static decomposition the Gini Coefficient by Income Sources
#'
#' @param data A data.frame with income sources and the total income (and possibliy sample weights) as variables
#' @param incomeSources_vars A character vector with names of the income sources (they must be columns of 'data')
#' @param totalIncome_var A character value naming the variable wich represents the total income (it must be a column of 'data' and for each observation in the data, the income sources must add up to the total income).
#' @param weight (optional) A character values naming the variable wich represents the sample weights
#' 
#' @return Returns a tibble with the following coluns:
#'     - Income Source (name of the income source)
#'     - Concentration Coefficient (with respect to the total income)
#'     - Share (with respect to the total income)
#'     - Concentration Coefficient*Share (sums up to the Gini Coefficient)
#'     - Contribution to Gini Coefficient
#' 
#' @import tibble
#' @import dplyr
#' @import purrr
#' 
#' @export
decomp_giniBySources <- function(data, incomeSources_vars, totalIncome_var, weight = NULL){
        
        if(is.null(weight)){
                vars <- c(incomeSources_vars, totalIncome_var)
        }else{
                vars <- c(incomeSources_vars, totalIncome_var, weight)
        }
        
        data <- data %>%
                select_at(vars) %>%
                mutate_at(.vars = incomeSources_vars, inequalityTools:::replaceNAvalues) %>%
                filter(complete.cases(.))
        
        X <- data %>%
                select_at(incomeSources_vars)
        
        y <- data[[totalIncome_var]]
        
        
        if(is.null(weight)){
                w = rep(1, nrow(data))
        }else{
                w <- data[[weight]]		
        }
        
        concentration <- purrr::map_dfr(X, function(x_i){
                calc_concCoef(x = x_i, y = y, w = w)
        }) %>%
                as.matrix() %>%
                t()
        
        share <- purrr::map_dfr(X, function(x_i) sum(x_i*w)/sum(y*w)) %>%
                as.matrix() %>%
                t()
        
        incomeSource   <- incomeSources_vars
        concentr_share <- concentration*share
        sumComponents  <- sum(concentr_share)
        contribution   <- concentr_share/sumComponents
        
        components <- tibble(incomeSources_vars, concentration, share, concentr_share, contribution)
        names(components) <- c("Income Source",
                               "Concentration Coefficient", 
                               "Share", 
                               "Concentration Coefficient*Share", 
                               "Relative Contribution to Gini Coefficient")
        
        gini <- calc_giniCoef(x = y, w = w)
        
        if(round(sumComponents, 5) != round(gini, 5)){
                warning("Components do not add up to the actual Gini")
        }
        
        components
}

#' Dynamic decomposition of the Gini Coefficient by Income Sources
#'
#' @param decomp_giniBySources0 A data.frame with with results of a static decomposition the Gini (produced by decomp_giniBySources) for group/time 0
#' @param decomp_giniBySources1 A data.frame with with results of a static decomposition the Gini (produced by decomp_giniBySources) for group/time 1
#' 
#' @return Returns a tibble with the following coluns:
#'     - Income Source (name of the income source)
#'     - Concentration Coefficient (with respect to the total income)
#'     - Share (with respect to the total income)
#'     - Concentration Coefficient*Share (sums up to the Gini Coefficient)
#'     - Contribution to Gini Coefficient
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export
dynamicDecomp_giniBySources <- function(decomp_giniBySources0, decomp_giniBySources1){
        
        components_0  <- as.matrix(decomp_giniBySources0[,2:3])
        components_1  <- as.matrix(decomp_giniBySources1[,2:3])
        
        incomeSources <- decomp_giniBySources1[[1]]
                
        diff <- components_1 - components_0
        mean <- (components_1 + components_0)/2
        
        concentration_effect <- diff[,1]*mean[,2]
        share_effect         <- diff[,2]*mean[,1]
        
        GiniDiff <- sum(concentration_effect) + sum(share_effect)
        
        effects <- cbind(concentration_effect, share_effect)
        effects <- cbind(effects,  concentration_effect + share_effect)
        effects <- rbind(effects,  colSums(effects))
        
        colnames(effects)[3]             <- "TOTAL"
        rownames(effects)[nrow(effects)] <- "TOTAL"
        
        col_names       <- colnames(effects)
        `Income Source` <- c(incomeSources, "TOTAL")
        
        effects = as_tibble(effects) %>% setNames(col_names)
        
        effects = bind_cols(tibble(`Income Source`), 
                            effects)
        
        effects
}