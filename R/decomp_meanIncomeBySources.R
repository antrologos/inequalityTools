#' Static decomposition of the mean income by sources
#'
#' @param data A data.frame with income sources and the total income (and possibliy sample weights) as variables
#' @param incomeSources_vars A character vector with names of the income sources (they must be columns of 'data')
#' @param totalIncome_var A character value naming the variable wich represents the total income (it must be a column of 'data' and for each observation in the data, the income sources must add up to the total income).
#' @param weight (optional) A character values naming the variable wich represents the sample weights
#' 
#' @return Returns a tibble with the following coluns:
#'     - Income Source: name of the income source
#'     - Mean: Mean of the income source
#'     - Prop. Pop. Recipients: Proportion of the population receiving non-zero income for each source
#'     - Mean*Prop. Pop. Recipients: mean times the proportion of recipients. This column must sum up to the global mean of the total income variable.
#' 
#' @import tibble
#' @import dplyr
#' @import purrr
#' 
#' @export
decomp_meanIncomeBySources <- function(data, incomeSources_vars, totalIncome_var, weight = NULL){
        
        if(is.null(weight)){
                vars <- c(incomeSources_vars, totalIncome_var)
        }else{
                vars <- c(incomeSources_vars, totalIncome_var, weight)
        }
        
        data <- data %>%
                select_at(vars) %>%
                filter_at(.vars = c(totalIncome_var, weight), 
                          .vars_predicate = all_vars(!is.na(.)))
        
        X <- data %>%
                select_at(incomeSources_vars) %>%
                mutate_all(function(x) { 
                        x[x %in% 0] <- NA
                        x
                        })
        
        y <- data[[totalIncome_var]]
        
        
        if(is.null(weight)){
                w = rep(1, nrow(data))
        }else{
                w <- data[[weight]]		
        }
        
        data <- bind_cols(X, tibble(w))
        
        data <- left_join(data %>%
                                  summarise_at(.vars = c(incomeSources_vars),
                                               function(x) wtd.mean(x, w)) %>%
                                  inequalityTools:::transpose(hasColNames = F) %>% setNames(c("source","mean")),
                          
                          data %>%
                                  summarise_at(.vars = c(incomeSources_vars),
                                               function(x) wtd.mean(!is.na(x), w)) %>%
                                  inequalityTools:::transpose(hasColNames = F) %>% setNames(c("source","p")),
                          
                          by = "source") %>%
                mutate(m_p = mean*p)
        
        test <- identical(round(sum(y*w)/sum(w), 6), round(sum(data$m_p), 6))
        
        if(!test){
                stop("Components do not sum up to total income")
        }
        
        data %>%
                rename(`Income Source`              = source,
                       Mean                         = mean,
                       `Prop. Pop. Recipients`      = p,
                       `Mean*Prop. Pop. Recipients` = m_p)
}


#' Dynamic decomposition of the mean income by sources
#'
#' @param decomp_meanIncomeBySources0 A data.frame with with results of a static decomposition the mean income (produced by decomp_meanIncomeBySources) for group/time 0
#' @param decomp_meanIncomeBySources1 A data.frame with with results of a static decomposition the mean income (produced by decomp_meanIncomeBySources) for group/time 1
#' 
#' @return Returns a tibble with the following coluns:
#'     - Income Source: name of the income source
#'     - mean_effect: change in mean total income due to the change in the mean of each income source
#'     - recipients_effect: change in mean total income due to the change in proportion of recipients of each income source
#'     - TOTAL: total change in mean total income due to each income source (mean_effect + recipients_effect)
#' 
#' @import tibble
#' @import dplyr
#' 
#' @export
dynamicDecomp_meanIncomeBySources <- function(decomp_meanIncomeBySources0, 
                                              decomp_meanIncomeBySources1){
        
        components_0 = as.matrix(decomp_meanIncomeBySources0[, -c(1,4)])
        components_1 = as.matrix(decomp_meanIncomeBySources1[, -c(1,4)])
        
        incomeSources <- decomp_meanIncomeBySources0[[1]]
        
        diff <- components_1 - components_0
        mean <- (components_1 + components_0)/2
        
        mean_effect        <- diff[,1]*mean[,2]
        recipients_effect  <- diff[,2]*mean[,1]
        
        effects <- cbind(mean_effect, recipients_effect)
        effects <- cbind(effects,     mean_effect + recipients_effect)
        effects <- rbind(effects,     colSums(effects))
        
        colnames(effects)[3] <- "TOTAL"
        rownames(effects)    <- c(incomeSources, "TOTAL")
        
        col_names       <- colnames(effects)
        `Income Source` <- rownames(effects)
        
        effects = as_tibble(effects) %>% 
                setNames(col_names)
        
        effects = bind_cols(tibble(`Income Source`), 
                            effects)
        
        effects
}



