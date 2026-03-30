#' G-computation estimator for population-average vaccine effect
#'
#' Computes the marginal average treatment effect of vaccination on the
#' post-infection outcome using g-computation:
#' \deqn{E[Y(1) - Y(0)] = E\{E[Y \mid Z = 1, X]\} - E\{E[Y \mid Z = 0, X]\}}
#'
#' Optionally supports a two-part model that decomposes the outcome into
#' components conditional on an intermediate variable.
#'
#' @param data A data.frame containing the observed data.
#' @param models A named list of fitted models required for estimation. Must include:
#' \describe{
#'   \item{fit_Y_Z_X}{Outcome regression model \eqn{E[Y \mid Z, X]} (if `two_part_model = FALSE`)}
#'   \item{fit_Y_Z0_S0_X, fit_Y_Z0_S1_X, fit_Y_Z1_S0_X, fit_Y_Z1_S1_X}{Outcome models for each (Z, S) combination (if `two_part_model = TRUE`)}
#'   \item{fit_S_Z0_X, fit_S_Z1_X}{Models for \eqn{P(S = 1 \mid Z, X)} (if `two_part_model = TRUE`)}
#' }
#' @param two_part_model Logical; whether to use a two-part outcome model. Default is FALSE.
#' @param Z_name Character; name of the treatment variable. Default is `"Z"`.
#' @param X_name Character vector; names of covariates. Default is `"X"`.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{additive_effect}{Estimated additive effect \eqn{E[Y(1) - Y(0)]}}
#'   \item{log_multiplicative_effect}{Log multiplicative effect \eqn{\log(E[Y(1)] / E[Y(0)])}}
#'   \item{psi_1}{Estimated mean outcome under treatment}
#'   \item{psi_0}{Estimated mean outcome under control}
#' }
#'
#' @export
do_gcomp_pop <- function(data, 
                          models,
                          two_part_model = FALSE,
                          Z_name = "Z",
                          X_name = c("X")){
  # E[Y(1) - Y(0)] = E[Y(1)] - E[Y(0)] = E[E[Y | Z = 1, X = x]] - E[E[Y | Z = 0, X = x]]
  
  df_Z1 <- data.frame(Z = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_Z1) <- c(Z_name, X_name)
  df_Z0 <- data.frame(Z = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
  names(df_Z0) <- c(Z_name, X_name)
  
  if(!two_part_model){
    E_Y_Z1_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z1)
    E_Y_Z0_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z0)
  }else{
    E_Y_Z0_S0_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
    E_Y_Z0_S1_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)

    E_Y_Z1_S0_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    E_Y_Z1_S1_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)

    E_Y_Z1_X <- E_Y_Z1_S1_X * rho_1_X + E_Y_Z1_S0_X * (1 - rho_1_X)
    E_Y_Z0_X <- E_Y_Z0_S1_X * rho_0_X + E_Y_Z0_S0_X * (1 - rho_0_X)
  }

  psi_1 <- mean(E_Y_Z1_X)
  psi_0 <- mean(E_Y_Z0_X)

  pop_growth_effect <- psi_1 - psi_0
  
  pop_growth_effect_log_mult <- log(psi_1 / psi_0)
  
  out <- c(pop_growth_effect, pop_growth_effect_log_mult, psi_1, psi_0)
  names(out) <-  c("additive_effect","log_multiplicative_effect", "psi_1", "psi_0")
  return(out)
}

#' Inverse probability weighting (IPW) estimator for population-average vaccine effect
#'
#' Estimates the marginal average treatment effect using inverse probability
#' weighting:
#' \deqn{E[Y(1)] = E\left[\frac{Z Y}{P(Z = 1 \mid X)}\right], \quad
#'       E[Y(0)] = E\left[\frac{(1 - Z) Y}{P(Z = 0 \mid X)}\right]}
#'
#' @param data A data.frame containing the observed data.
#' @param models A named list containing the propensity score model:
#' \describe{
#'   \item{fit_Z_X}{Model for \eqn{P(Z = 1 \mid X)}}
#' }
#' @param Z_name Character; name of the treatment variable. Default is `"Z"`.
#' @param Y_name Character; name of the outcome variable. Default is `"Y"`.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{additive_effect}{Estimated additive effect}
#'   \item{log_multiplicative_effect}{Log multiplicative effect}
#'   \item{psi_1}{Estimated mean outcome under treatment}
#'   \item{psi_0}{Estimated mean outcome under control}
#' }
do_ipw_pop <- function(data, models, Z_name = "Z", Y_name = "Y"){
  
  pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
  pi_0_X <- 1 - pi_1_X
  
  Y <- data[[Y_name]]
  Z <- data[[Z_name]]
  
  psi_1_ipw <- mean(
    Z / pi_1_X * Y
  ) 
  psi_0_ipw <- mean(
    (1 - Z) / pi_0_X * Y
  ) 
  
  pop_growth_effect <- psi_1_ipw - psi_0_ipw
  pop_growth_effect_log_mult <- log(psi_1_ipw / psi_0_ipw)
  
  out <- c(pop_growth_effect, pop_growth_effect_log_mult, psi_1_ipw, psi_0_ipw)
  names(out) <-  c("additive_effect","log_multiplicative_effect", "psi_1", "psi_0")
  return(out)
}

#' Augmented inverse probability weighting (AIPW) estimator for population-average growth effect
#'
#' Computes a doubly robust estimate of the marginal treatment effect by combining
#' outcome regression and inverse probability weighting. The estimator remains
#' consistent if either the outcome model or the treatment model is correctly specified.
#'
#' Supports both standard and two-part outcome models.
#'
#' @param data A data.frame containing the observed data.
#' @param models A named list of fitted models required for estimation. Must include:
#' \describe{
#'   \item{fit_Y_Z_X}{Outcome model \eqn{E[Y \mid Z, X]} (if `two_part_model = FALSE`)}
#'   \item{fit_Z_X}{Treatment model \eqn{P(Z = 1 \mid X)}}
#'   \item{fit_Y_Z0_S0_X, fit_Y_Z0_S1_X, fit_Y_Z1_S0_X, fit_Y_Z1_S1_X}{Outcome models (if `two_part_model = TRUE`)}
#'   \item{fit_S_Z0_X, fit_S_Z1_X}{Models for \eqn{P(S = 1 \mid Z, X)} (if `two_part_model = TRUE`)}
#' }
#' @param Z_name Character; name of the treatment variable. Default is `"Z"`.
#' @param X_name Character vector; names of covariates. Default is `"X"`.
#' @param Y_name Character; name of the outcome variable. Default is `"Y"`.
#' @param return_se Logical; whether to return standard errors. Default is TRUE.
#' @param two_part_model Logical; whether to use a two-part outcome model. Default is FALSE.
#'
#' @return If `return_se = TRUE`, a named numeric vector containing:
#' \describe{
#'   \item{additive_effect}{Estimated additive effect}
#'   \item{additive_se}{Standard error of additive effect}
#'   \item{log_multiplicative_effect}{Log multiplicative effect}
#'   \item{log_multiplicative_se}{Standard error on log scale}
#'   \item{psi_1}{Estimated mean outcome under treatment}
#'   \item{se_psi_1}{Standard error of \eqn{\psi_1}}
#'   \item{psi_0}{Estimated mean outcome under control}
#'   \item{se_psi_0}{Standard error of \eqn{\psi_0}}
#' }
#'
#' If `return_se = FALSE`, returns only the additive and log multiplicative effects.
#'
#' The returned object also includes an `"if_matrix"` attribute containing the
#' estimated influence function contributions for \eqn{\psi_1} and \eqn{\psi_0}.
#'
#' @export
do_aipw_pop <- function(
    data, 
    models,
    Z_name = "Z",
    X_name = c("X"),
    Y_name = "Y",
    return_se = TRUE,
    two_part_model = FALSE
){
  
  
  if(!two_part_model){
    df_Z1 <- data.frame(Z = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
    names(df_Z1) <- c(Z_name, X_name)
    
    Qbar_Z1 <- simple_predict(models$fit_Y_Z_X, newdata = df_Z1)
    
    df_Z0 <- data.frame(Z = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
    names(df_Z0) <- c(Z_name, X_name)
    
    Qbar_Z0 <- simple_predict(models$fit_Y_Z_X, newdata = df_Z0)
    
  } else{
    
    E_Y_Z0_S0_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
    E_Y_Z0_S1_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    
    E_Y_Z1_S0_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    E_Y_Z1_S1_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
    
    Qbar_Z1 <- E_Y_Z1_S1_X * rho_1_X + E_Y_Z1_S0_X * (1 - rho_1_X)
    Qbar_Z0 <- E_Y_Z0_S1_X * rho_0_X + E_Y_Z0_S0_X * (1 - rho_0_X)
  }
  
  pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
  pi_0_X <- 1 - pi_1_X
  
  Y <- data[[Y_name]]
  Z <- data[[Z_name]]
  
  psi_1_plugin <- mean(Qbar_Z1)
  augmentation_1 <- Z / pi_1_X * ( Y - Qbar_Z1 ) + Qbar_Z1 - psi_1_plugin
  psi_1_aipw <- psi_1_plugin + mean(augmentation_1)
  
  psi_0_plugin <- mean(Qbar_Z0)
  augmentation_0 <- (1 - Z) / pi_0_X * ( Y - Qbar_Z0 ) + Qbar_Z0 - psi_0_plugin
  psi_0_aipw <- psi_0_plugin + mean(augmentation_0)
  
  # Additive effect
  n <- dim(data)[1]
  growth_effect <- psi_1_aipw - psi_0_aipw
  se <- sqrt(var(augmentation_1 - augmentation_0) / n)
  
  se_psi_1 <- sqrt(var(augmentation_1) / n)
  se_psi_0 <- sqrt(var(augmentation_0) / n)
  
  # Multiplicative effect (log scale)
  growth_effect_log_mult <- log(psi_1_aipw / psi_0_aipw)
  
  # Yet SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  cov_matrix <- stats::cov(if_matrix) / n
  
  gradient <- matrix(c(1 / psi_1_aipw, -1 / psi_0_aipw), ncol = 1)
  
  se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
  
  if(return_se){
    out <- c(growth_effect, se, growth_effect_log_mult, se_log_mult_eff, psi_1_aipw, se_psi_1, psi_0_aipw, se_psi_0)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se", "psi_1", "se_psi_1", "psi_0", "se_psi_0")
    
    # added to return influence fn matrix for vaccine trial project without having to restructure whole package rn. eventually change return type to list
    attr(out, "if_matrix") <- if_matrix
    
    return(out)
  }else{
    out <- c(growth_effect, growth_effect_log_mult)
    names(out) <- c("additive_effect", "log_multiplicative_effect")
    return(out)
  }
  
}
