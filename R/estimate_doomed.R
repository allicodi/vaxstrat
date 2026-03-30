#' G-computation for the doomed principal stratum
#'
#' Estimates the post-infection outcomes among individuals who would be
#' infected regardless of treatment assignment (the "doomed" principal stratum)
#' using g-computation.
#'
#' @param data A \code{data.frame} containing the dataset to predict on.
#' @param models A list of pre-fitted models as returned by \code{fit_models()}
#' or \code{fit_ml_models()}.
#'
#' @returns A named numeric vector containing:
#' \describe{
#'   \item{additive_effect}{Estimated additive effect (difference in means).}
#'   \item{log_multiplicative_effect}{Estimated log multiplicative effect.}
#'   \item{psi_1}{Estimated mean counterfactual outcome under treatment (Z = 1).}
#'   \item{psi_0}{Estimated mean counterfactual outcome under control (Z = 0).}
#' }
#' 
#' @export
do_gcomp_doomed <- function(data, models){
  
  # Psi_1 = E[P(S=1 | Z = 0, X) / P(S= 1 | Z = 0) * E[Y | Z=1, X] ]
  mu_11_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
  mu_10_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
  rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
  rho_bar_1 <- mean(rho_1_X)
  
  psi_1 <- mean(
    rho_1_X / rho_bar_1 * mu_11_X
  )
  psi_0 <- mean(
    rho_1_X / rho_bar_1 * mu_10_X
  )
  
  growth_effect <- psi_1 - psi_0
  growth_effect_log_mult <- log(psi_1 / psi_0)
  
  out <- c(growth_effect, growth_effect_log_mult, psi_1, psi_0)
  names(out) <- c("additive_effect","log_multiplicative_effect", "psi_1", "psi_0")
  
  return(out)
}

#' IPW estimator for counterfactual outcomes in the doomed principal stratum
#'
#' Estimates the average counterfactual outcomes among individuals in the doomed
#' principal stratum using inverse probability weighting (IPW).
#'
#' @param data A \code{data.frame} containing the dataset to predict on.
#' @param models A list of pre-fitted models as returned by \code{fit_models()}
#' or \code{fit_ml_models()}.
#' @param S_name Character; name of the infection indicator variable.
#' @param Y_name Character; name of the outcome variable.
#' @param Z_name Character; name of the treatment variable.
#'
#' @returns A named numeric vector containing:
#' \describe{
#'   \item{additive_effect}{Estimated additive effect.}
#'   \item{log_multiplicative_effect}{Estimated log multiplicative effect.}
#'   \item{psi_1}{Estimated mean counterfactual outcome under treatment.}
#'   \item{psi_0}{Estimated mean counterfactual outcome under control.}
#' }
#' 
#' @export
do_ipw_doomed <- function(
    data, models,
    S_name, Y_name, Z_name
){
  
  # Psi_1 = E[P(S=1 | Z = 0, X) / P(S= 1 | Z = 0) * E[Y | Z=1, X] ]
  
  rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
  rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
  pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
  pi_0_X <- 1 - pi_1_X
  rho_bar_1 <- mean(rho_1_X)
  
  S <- data[[S_name]]
  Y <- data[[Y_name]]
  Z <- data[[Z_name]]
  
  eta_1 <- mean(
    ( Z / pi_1_X ) * ( S / rho_bar_1 ) * Y
  )
  
  eta_0 <- mean(
    ( (1 - Z) / pi_0_X ) * ( S / rho_0_X ) * ( rho_1_X / rho_bar_1 ) * Y
  )
  
  growth_effect <- eta_1 - eta_0
  growth_effect_log_mult <- log(eta_1 / eta_0)
  
  out <- c(growth_effect, growth_effect_log_mult, eta_1, eta_0)
  names(out) <- c("additive_effect","log_multiplicative_effect", "psi_1", "psi_0")
  
  return(out)
}

#' Augmented IPW estimator for doomed principal stratum
#'
#' Computes an efficient augmented inverse probability weighted (AIPW) estimator
#' for the vaccine effect in the doomed principal stratum.
#' This estimator is doubly robust: it is consistent if either the outcome model
#' or the treatment/infection models are correctly specified.
#'
#' @param data A \code{data.frame} containing the dataset to predict on.
#' @param models A list of pre-fitted models as returned by \code{fit_models()}
#' or \code{fit_ml_models()}.
#' @param Y_name Character; name of the outcome variable. Default is \code{"Y"}.
#' @param Z_name Character; name of the treatment variable. Default is \code{"Z"}.
#' @param S_name Character; name of the infection variable. Default is \code{"S"}.
#' @param return_se Logical; whether to return closed form standard errors. Default is \code{FALSE}.
#'
#' @returns If \code{return_se = FALSE}, a named numeric vector containing:
#' \describe{
#'   \item{additive_effect}{Estimated additive effect.}
#'   \item{log_multiplicative_effect}{Estimated log multiplicative effect.}
#' }
#'
#' If \code{return_se = TRUE}, returns:
#' \describe{
#'   \item{additive_effect}{Estimated additive effect.}
#'   \item{additive_se}{Standard error for additive effect.}
#'   \item{log_multiplicative_effect}{Estimated log multiplicative effect.}
#'   \item{log_multiplicative_se}{Standard error for log multiplicative effect.}
#'   \item{psi_1}{Estimated mean under treatment.}
#'   \item{se_psi_1}{Standard error for \code{psi_1}.}
#'   \item{psi_0}{Estimated mean under control.}
#'   \item{se_psi_0}{Standard error for \code{psi_0}.}
#' }
#' 
#' @export
do_aipw_doomed <- function(data, 
                           models,
                           Y_name = "Y",
                           Z_name = "Z",
                           S_name = "S",
                           return_se = FALSE){
  
  
  mu_11_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
  mu_01_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
  
  rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
  rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
  
  pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
  pi_0_X <- 1 - pi_1_X
  rho_bar_1 <- mean(rho_1_X)
  
  eta_tilde_0 <- rho_1_X * mu_01_X / rho_bar_1
  eta_tilde_1 <- rho_1_X * mu_11_X / rho_bar_1
  
  eta_0_plugin <- mean( eta_tilde_0 )
  eta_1_plugin <- mean( eta_tilde_1 )
  
  Z <- data[[Z_name]]
  S <- data[[S_name]]
  Y <- data[[Y_name]]
  
  augmentation_0 <- (
    ( (1 - Z) / pi_0_X ) * ( rho_1_X / rho_bar_1 ) * ( S / rho_0_X ) * ( Y - mu_01_X ) + 
      (mu_01_X - eta_0_plugin) / rho_bar_1 * ( Z / pi_1_X ) * ( S - rho_1_X ) -
      eta_0_plugin / rho_bar_1 * ( rho_1_X - rho_bar_1 ) + eta_tilde_0 - eta_0_plugin
  )
  augmentation_1 <- (
    ( Z / pi_1_X ) * ( S / rho_bar_1 ) * ( Y - mu_11_X ) + 
      ( mu_11_X - eta_1_plugin ) / rho_bar_1 * ( Z / pi_1_X ) * ( S - rho_1_X ) - 
      eta_1_plugin / rho_bar_1 * ( rho_1_X - rho_bar_1 ) + eta_tilde_1 - eta_1_plugin
    
  )
  
  eta_0_aipw <- eta_0_plugin + mean(augmentation_0)
  eta_1_aipw <- eta_1_plugin + mean(augmentation_1)
  
  
  # Additive effect
  growth_effect <- eta_1_aipw - eta_0_aipw
  se <- sqrt(stats::var(augmentation_1 - augmentation_0) / dim(data)[1])
  
  se_psi_1 <- sqrt(stats::var(augmentation_1) / dim(data)[1])
  se_psi_0 <- sqrt(stats::var(augmentation_0) / dim(data)[1])
  
  # Multiplicative effect (log scale)
  growth_effect_log_mult <- log(eta_1_aipw / eta_0_aipw)
  
  # Yet SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  cov_matrix <- stats::cov(if_matrix) / dim(data)[1]
  
  gradient <- matrix(c(1 / eta_1_aipw, -1 / eta_0_aipw), ncol = 1)
  
  se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
  
  if(return_se){
    out <- c(growth_effect, se, growth_effect_log_mult, se_log_mult_eff, eta_1_aipw, se_psi_1, eta_0_aipw, se_psi_0)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se", "psi_1", "se_psi_1", "psi_0", "se_psi_0")
    return(out)
  }else{
    out <- c(growth_effect, growth_effect_log_mult)
    names(out) <- c("additive_effect", "log_multiplicative_effect")
    return(out)
  }
}

#' Bounds for doomed principal stratum without cross-world assumptions
#'
#' Computes nonparametric bounds for individuals in the doomed principal stratum, 
#' without invoking cross-world assumptions. 
#'
#' @param data A \code{data.frame} containing the dataset.
#' @param Y_name Character; name of the outcome variable. Default is \code{"Y"}.
#' @param Z_name Character; name of the treatment variable. Default is \code{"Z"}.
#' @param S_name Character; name of the infection variable. Default is \code{"S"}.
#' @param family Character; outcome type. Either \code{"gaussian"} for continuous
#' outcomes or \code{"binomial"} for binary outcomes.
#'
#' @details
#' The method assumes that the infection rate is lower under treatment than control.
#' If this condition is not met, the function returns \code{NA} values.
#'
#' @returns A named numeric vector containing:
#' \describe{
#'   \item{E_Y1__S0_1}{Observed mean outcome among treated and infected.}
#'   \item{E_Y0__S0_1_lower}{Lower bound.}
#'   \item{E_Y0__S0_1_upper}{Upper bound.}
#'   \item{additive_effect_lower}{Lower bound on additive effect.}
#'   \item{additive_effect_upper}{Upper bound on additive effect.}
#'   \item{mult_effect_lower}{Lower bound on multiplicative effect.}
#'   \item{mult_effect_upper}{Upper bound on multiplicative effect.}
#' }
#'
#' @export
get_bound_doomed <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    family = "gaussian"
){
  
  # Step 1: rhobar_z_n
  
  # 1.1 rhobar_0_n (or mean in subset) -- size of infected in unvaccinated (naturally infected, doomed + protected)
  rhobar_0_n <- mean(data[[S_name]][data[[Z_name]] == 0])
  
  # 1.2 rhobar_1_n -- size of infected in vaccinated (doomed)
  rhobar_1_n <- mean(data[[S_name]][data[[Z_name]] == 1])
  
  if(rhobar_0_n > rhobar_1_n){
    
    # Step 2: q_n (relative size of protected in naturally infected (doomed + protected) in unvax)
    q_n = 1 - rhobar_1_n / rhobar_0_n
    
    # Step 3: q_n^th quintiles of Y__Z0_S1 
    Y__Z0_S1 <- data[[Y_name]][which(data[[Z_name]] == 0 & data[[S_name]] == 1)]
    q_nth_quintile <- quantile(Y__Z0_S1, probs = q_n)
    one_minus_q_nth_quintile <- quantile(Y__Z0_S1, probs = 1 - q_n)
    
    # Step 4: mubar_10_l,u_n 
    if(family == "gaussian"){
      # get people < quintile lower (protected people)
      mubar_10_u_n <- sum(data[[Y_name]] * as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] > q_nth_quintile )) / 
        sum(as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] > q_nth_quintile ))
      
      # get people > quintile upper (protected people)
      mubar_10_l_n <- sum(data[[Y_name]] * as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] < one_minus_q_nth_quintile )) / 
        sum(as.numeric(data[[S_name]] == 1 & data[[Z_name]] == 0 & data[[Y_name]] < one_minus_q_nth_quintile ))
    } else{
      # Binary outcome
      
      # Unvaccinated, infected (naturally infected in placebo arm)
      data__Z0_S1 <- data[which(data[[Z_name]] == 0 & data[[S_name]] == 1),]
      
      # q_n = proportion of people in naturally infected, placebo arm that are protected
      # 1 - q_n = proportion of people in naturally infected, placebo arm that are doomed
      # target_num = proportion of doomed people * number of people in naturally infected = target number of doomed people
      target_num <- ceiling((1 - q_n)*nrow(data__Z0_S1)) # target number of doomed in placebo arm
      num_0s <- length(which(data__Z0_S1[[Y_name]] == 0))  # observed number of 0s in Doomed + Protected (naturally infected) box in the placebo arm
      num_1s <- length(which(data__Z0_S1[[Y_name]] == 1))  # observed number of 1s in Doomed + Protected (naturally infected) box in the placebo arm
      
      prop_0s__Z0_S1 <- 1 - mean(data__Z0_S1[[Y_name]])
      # check equivalence
      (num_0s >= target_num) == (prop_0s__Z0_S1 > rhobar_1_n / rhobar_0_n)
      ## Lower Bound:
      
      # Check if at least target_num 0s in the vax uninfected 
      if(num_0s >= target_num){
        # If so, muhat_10_l = 0
        mubar_10_l_n <- 0
        mubar_10_l_n_2 <- 0
      } else{
        # Else, determine fraction of 1s that need to be kept
        
        # (number of 0s to be size of doomed - num 0s observed) / number of 0s to be size of doomed = proportion of 1s to be added
        mubar_10_l_n <- ( (target_num - num_0s) / target_num )
      }
      
      ## Upper Bound:
      
      # Check if at least q_n * 100 % 1s in the vax uninfected 
      if(num_1s >= target_num){
        # If so, mubar_10_u = 1
        mubar_10_u_n <- 1
        mubar_10_u_n_2 <- 1
      } else{
        # Else, proportion of 1s in doomed box
        
        mubar_10_u_n <- num_1s / target_num
        
      }
      
    }
    
    l_n <- mubar_10_l_n
    u_n <- mubar_10_u_n
    
    #mean in vaccinated infecteds for comparison
    E_Y1__S0_1 <- mean(data[[Y_name]][data[[S_name]] == 1 & data[[Z_name]] == 1])
    
    out <- c(
      E_Y1__S0_1,
      l_n,
      u_n,
      # flip to get bounds on the effect print properly in table 
      E_Y1__S0_1 - u_n,
      E_Y1__S0_1 - l_n,
      E_Y1__S0_1 / u_n,
      E_Y1__S0_1 / l_n
    )
    
    names(out) <- c("E_Y1__S0_1",
                    "E_Y0__S0_1_lower",
                    "E_Y0__S0_1_upper",
                    "additive_effect_lower",
                    "additive_effect_upper",
                    "mult_effect_lower",
                    "mult_effect_upper")
    
  } else{
    # stop("Method not applicable unless evidence of vaccine protection.")
    
    # Get rid of this condition because permutation test
    # stop("Method not applicable unless evidence of vaccine protection.")
    out <- rep(NA, 7)
    
    names(out) <- c("E_Y1__S0_1",
                    "E_Y0__S0_1_lower",
                    "E_Y0__S0_1_upper",
                    "additive_effect_lower",
                    "additive_effect_upper",
                    "mult_effect_lower",
                    "mult_effect_upper")
    
  }

  return(out)
  
}