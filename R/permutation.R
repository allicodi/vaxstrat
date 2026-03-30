
#' Permutation test for bounds in the naturally infected stratum
#'
#' Performs a permutation-based hypothesis test for bounds on the treatment effect
#' in the naturally infected principal stratum. Treatment assignment is permuted
#' to generate the null distribution of the bound-based estimator.
#'
#' @param data A \code{data.frame} containing the dataset.
#' @param Y_name Character; name of the outcome variable. Default is \code{"Y"}.
#' @param Z_name Character; name of the treatment variable. Default is \code{"Z"}.
#' @param S_name Character; name of the infection variable. Default is \code{"S"}.
#' @param n_permutations Integer; number of permutations to perform. Default is \code{1000}.
#' @param family Character; outcome type. Either \code{"gaussian"} (continuous)
#' or \code{"binomial"} (binary). Default is \code{"gaussian"}.
#' @param effect_dir Character; direction of beneficial effect. Either
#' \code{"positive"} (default) or \code{"negative"}. Determines which bound is used
#' for the one-sided hypothesis test.
#'
#' @details
#' The permutation test constructs a null distribution by randomly permuting the
#' treatment assignment variable (\code{Z}). The observed bound is then compared
#' to this null distribution using a one-sided test:
#' \itemize{
#'   \item If \code{effect_dir = "positive"}, the lower bound is tested against 0.
#'   \item If \code{effect_dir = "negative"}, the upper bound is tested against 0.
#' }
#'
#' Any permutation replicates that produce undefined bounds (e.g., due to violations
#' of identifying assumptions) are treated as not supporting the alternative hypothesis.
#'
#' @returns A list with class \code{"permutation_bound_nat_inf"} containing:
#' \describe{
#'   \item{original_est}{Named vector of bound estimates from the observed data.}
#'   \item{null_est}{Matrix of bound estimates from permuted datasets (rows correspond to permutations).}
#'   \item{pval_bound}{One-sided permutation p-value based on the specified direction.}
#' }
permutation_bound_nat_inf <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    n_permutations = 1e3, 
    family = "gaussian",
    effect_dir = "positive"
){
  
  original_est <- get_bound_nat_inf(
    data = data, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
  )
  
  ## Permutation approach
  null_est <- vector("list", length = n_permutations)
  for(i in 1:n_permutations){
    data_shuffle <- data
    data_shuffle[[Z_name]] <- sample(data_shuffle[[Z_name]])
    
    null_est[[i]] <- get_bound_nat_inf(
      data = data_shuffle, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
    )
  }
  
  null_est <- do.call(rbind, null_est)
  
  pval_vec <- ifelse(rep(effect_dir, nrow(null_est)) == "negative",
                     null_est[,'additive_effect_upper'] < original_est['additive_effect_upper'], # negative effect means we are interested in the upper bound < 0
                     null_est[,'additive_effect_lower'] > original_est['additive_effect_lower'])  # positive effect means we are interested in the lower bound > 0 
  
  # Temp fill in NA with FALSE? for the ones that failed to meet condition rhobar_0_n > rhobar_1_n
  pval_vec[is.na(pval_vec)] <- FALSE
                     
  out <- list(
    original_est = original_est, 
    null_est = null_est,
    pval_bound = mean(as.numeric(pval_vec))
  )
    
    class(out) <- "permutation_bound_nat_inf"
    return(out)
    
}

#' Permutation test for bounds in the doomed principal stratum
#'
#' Performs a permutation-based hypothesis test for bounds on the treatment effect
#' in the doomed principal stratum. Treatment assignment is permuted to generate
#' the null distribution of the bound-based estimator.
#'
#' @param data A \code{data.frame} containing the dataset.
#' @param Y_name Character; name of the outcome variable. Default is \code{"Y"}.
#' @param Z_name Character; name of the treatment variable. Default is \code{"Z"}.
#' @param S_name Character; name of the infection variable. Default is \code{"S"}.
#' @param n_permutations Integer; number of permutations to perform. Default is \code{1000}.
#' @param family Character; outcome type. Either \code{"gaussian"} (continuous)
#' or \code{"binomial"} (binary). Default is \code{"gaussian"}.
#' @param effect_dir Character; direction of beneficial effect. Either
#' \code{"positive"} (default) or \code{"negative"}. Determines which bound is used
#' for the one-sided hypothesis test.
#'
#' @details
#' The permutation test constructs a null distribution by randomly permuting the
#' treatment assignment variable (\code{Z}). The observed bound is compared to this
#' null distribution using a one-sided test:
#' \itemize{
#'   \item If \code{effect_dir = "positive"}, the lower bound is tested against 0.
#'   \item If \code{effect_dir = "negative"}, the upper bound is tested against 0.
#' }
#'
#' Permutation replicates that yield undefined bounds (e.g., when identifying
#' assumptions are not satisfied) are treated as not supporting the alternative.
#'
#' @returns A list with class \code{"permutation_bound_doomed"} containing:
#' \describe{
#'   \item{original_est}{Named vector of bound estimates from the observed data.}
#'   \item{null_est}{Matrix of bound estimates from permuted datasets.}
#'   \item{pval_bound}{One-sided permutation p-value.}
#' }
permutation_bound_doomed <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    n_permutations = 1e3, 
    family = "gaussian",
    effect_dir = "positive"
){
  
  original_est <- get_bound_doomed(
    data = data, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
  )
  
  ## Permutation approach
  null_est <- vector("list", length = n_permutations)
  for(i in 1:n_permutations){
    data_shuffle <- data
    data_shuffle[[Z_name]] <- sample(data_shuffle[[Z_name]])
    
    null_est[[i]] <- get_bound_doomed(
      data = data_shuffle, S_name = S_name, Y_name = Y_name, Z_name = Z_name, family = family
    )
  }
  
  null_est <- do.call(rbind, null_est)
  
  pval_vec <- ifelse(rep(effect_dir, nrow(null_est)) == "negative",
                     null_est[,'additive_effect_upper'] < original_est['additive_effect_upper'], # negative effect means we are interested in the upper bound < 0
                     null_est[,'additive_effect_lower'] > original_est['additive_effect_lower'])  # positive effect means we are interested in the lower bound > 0 
  
  # Temp fill in NA with FALSE? for the ones that failed to meet condition rhobar_0_n > rhobar_1_n
  pval_vec[is.na(pval_vec)] <- FALSE
  
  out <- list(
    original_est = original_est, 
    null_est = null_est,
    pval_bound = mean(as.numeric(pval_vec))
  )
    
    class(out) <- "permutation_bound_doomed"
    return(out)
    
}