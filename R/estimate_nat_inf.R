#' Unadjusted estimator for naturally infected principal stratum
#'
#' Computes covariate-unadjusted estimates of counterfactual post-infection
#' outcomes among the naturally infected principal stratum.
#'
#' @param data A data.frame containing observed data.
#' @param Z_name Character string specifying the treatment (vaccination) variable.
#' @param Y_name Character string specifying the outcome (e.g., growth) variable.
#' @param S_name Character string specifying the infection indicator.
#'
#' @return A named numeric vector containing:
#' \describe{
#'   \item{additive_effect}{Estimated difference \eqn{\psi_1 - \psi_0}.}
#'   \item{log_multiplicative_effect}{Estimated log ratio \eqn{\log(\psi_1 / \psi_0)}.}
#'   \item{psi_1}{Estimated mean outcome under treatment in the naturally infected stratum.}
#'   \item{psi_0}{Estimated mean outcome under control in the naturally infected stratum.}
#' }
#'
#' @export
do_unadj_nat_inf <- function(
  data, Z_name, Y_name, S_name
){
  EY1 <- mean(data[[Y_name]][data[[Z_name]] == 1])
  EY0 <- mean(data[[Y_name]][data[[Z_name]] == 0])
  mu_bar_01 <- mean(data[[Y_name]][data[[Z_name]] == 0 & data[[S_name]] == 1])
  rho_bar_0 <- mean(data[[S_name]][data[[Z_name]] == 0])

  psi_1 <- (EY1 - EY0) / rho_bar_0 + mu_bar_01
  psi_0 <- mu_bar_01

  growth_effect <- psi_1 - psi_0
  growth_effect_log_mult <- log(psi_1 / psi_0)

  out <- c(growth_effect, growth_effect_log_mult, psi_1, psi_0)
  names(out) <- c("additive_effect", "log_multiplicative_effect", "psi_1", "psi_0")
  return(out)
}

#' G-computation estimator for naturally infected principal stratum
#'
#' Computes covariate-adjusted estimates of counterfactual post-infection outcomes
#' using g-computation under specified identification assumptions.
#'
#' @param data A data.frame containing observed data.
#' @param models A named list of fitted nuisance models returned by \code{fit_models}.
#' @param Z_name Character string specifying the treatment variable.
#' @param X_name Character vector of covariate names.
#' @param exclusion_restriction Logical; whether the exclusion restriction assumption is imposed.
#' @param cross_world Logical; whether a cross-world independence assumption is imposed.
#' @param two_part_model Logical; whether outcome regression is specified via a two-part model.
#'
#' @details
#' At least one of \code{exclusion_restriction} or \code{cross_world} must be TRUE.
#' Different identification strategies are used depending on these assumptions.
#'
#' @return A named numeric vector with elements:
#' \describe{
#'   \item{additive_effect}{Estimated difference \eqn{\psi_1 - \psi_0}.}
#'   \item{log_multiplicative_effect}{Estimated log ratio \eqn{\log(\psi_1 / \psi_0)}.}
#'   \item{psi_1}{Estimated mean outcome under treatment.}
#'   \item{psi_0}{Estimated mean outcome under control.}
#' }
#'
#' @export
do_gcomp_nat_inf <- function(
  data, models, Z_name = NULL, X_name = NULL,
  exclusion_restriction = FALSE,
  cross_world = TRUE,
  two_part_model = FALSE){
  
  if(!exclusion_restriction & cross_world){
    # Psi_1 = E[P(S=1 | Z = 0, X) / P(Y = 1 | Z = 0) * E[Y | Z=1, X] ]
    E_Y_Z1_S1_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    E_Y_Z1_S0_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    P_S1_Z1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
    P_S1_Z0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
        
    P_S1_Z0 <- mean(P_S1_Z0_X)
    VE_X <- 1 - ( P_S1_Z1_X / P_S1_Z0_X )
    E_Y1_S01_X <- E_Y_Z1_S1_X * (1 - VE_X) + E_Y_Z1_S0_X * VE_X
    
    psi_1 <- mean(
      ( P_S1_Z0_X / P_S1_Z0 ) * E_Y1_S01_X
    )
    
    # Psi_0 = E[P(S=1 | Z = 0, X) / P(Y = 1 | Z = 0) * E[Y | Z=0, Y = 1, X] ]
    
    # Option 1 for estimation:
    # psi_0 <- mean(sub_Z0_S1$Y) 
    
    # Option 2 for estimation:
    E_Y_Z0_S1_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    
    psi_0 <- mean(
      ( P_S1_Z0_X / P_S1_Z0 ) * E_Y_Z0_S1_X
    )
  }else if(exclusion_restriction & !cross_world){
    df_Z1 <- data.frame(Z = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
    names(df_Z1) <- c(Z_name, X_name)

    df_Z0 <- data.frame(Z = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
    names(df_Z0) <- c(Z_name, X_name)
    
    E_Y_Z0_S1_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    
    if(!two_part_model){
      E_Y_Z1_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z1)
      E_Y_Z0_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z0)
    }else{
      E_Y_Z0_S0_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
      E_Y_Z1_S0_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
      E_Y_Z1_S1_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
      rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
      E_Y_Z1_X <- E_Y_Z1_S1_X * rho_1_X + E_Y_Z1_S0_X * (1 - rho_1_X)
      E_Y_Z0_X <- E_Y_Z0_S1_X * rho_0_X + E_Y_Z0_S0_X * (1 - rho_0_X)
    }

    rho_bar_0 <- mean(rho_0_X)
    psi_1 <- mean(E_Y_Z1_X - E_Y_Z0_X) / rho_bar_0 + mean(rho_0_X / rho_bar_0 * E_Y_Z0_S1_X)

    psi_0 <- mean(
      ( rho_0_X / rho_bar_0 ) * E_Y_Z0_S1_X
    )

  }else if(exclusion_restriction & cross_world){
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
    rho_bar_0 <- mean(rho_0_X)
    
    mu_01_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    mu_11_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    mu_10_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    mu_00_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
    
    pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
    pi_0_X <- 1 - pi_1_X
    
    # psi_0 = Weight * E[E[Y | Z = 0, Y = 1, X]]
    psi_tilde_0_X <- rho_0_X / rho_bar_0 * mu_01_X
    psi_0 <- mean( psi_tilde_0_X )
    
    # original
    # mu_dot0_X <- pi_1_X * mu_10_X + pi_0_X * mu_00_X
    
    # more generalizable
    # wt_1
    pi_1_S0_X <- (1 - rho_1_X) * pi_1_X /
      (((1 - rho_0_X)*pi_0_X) +
      (1 - rho_1_X) * pi_1_X)
    
    # wt_0
    pi_0_S0_X <- (1 - rho_0_X) * pi_0_X /
      (((1 - rho_0_X)*pi_0_X) +
      (1 - rho_1_X) * pi_1_X)
    
    mu_dot0_X <- pi_1_S0_X * mu_10_X + pi_0_S0_X * mu_00_X

    rho_bar_dot <- pi_1_X * rho_1_X + pi_0_X * rho_0_X

    psi_tilde_1_X <- rho_1_X / rho_bar_0 * mu_11_X + ( rho_0_X - rho_1_X ) / rho_bar_0 * mu_dot0_X
    psi_1 <- mean( psi_tilde_1_X )
    
  }else{
    stop("Must assume exclusion_restriction, cross_world, or both.")
  }
  growth_effect <- psi_1 - psi_0
  growth_effect_log_mult <- log(psi_1 / psi_0)
  
  out <- c(growth_effect, growth_effect_log_mult,psi_1,psi_0)
  names(out) <- c("additive_effect","log_multiplicative_effect","psi_1","psi_0")
  
  return(out)
}

#' Inverse probability weighting estimator for naturally infected principal stratum
#'
#' Computes estimates of counterfactual post-infection outcomes using inverse
#' probability weighting (IPW).
#'
#' @param data A data.frame containing observed data.
#' @param models A named list of fitted nuisance models.
#' @param exclusion_restriction Logical; whether the exclusion restriction assumption is imposed.
#' @param S_name Character string specifying the infection indicator.
#' @param Y_name Character string specifying the outcome variable.
#' @param Z_name Character string specifying the treatment variable.
#'
#' @details
#' If \code{exclusion_restriction = FALSE}, estimation relies on cross-world
#' assumptions. Otherwise, an alternative IPW representation is used.
#'
#' @return A named numeric vector with:
#' \describe{
#'   \item{additive_effect}{Estimated difference \eqn{\psi_1 - \psi_0}.}
#'   \item{log_multiplicative_effect}{Estimated log ratio \eqn{\log(\psi_1 / \psi_0)}.}
#'   \item{psi_1}{Estimated mean outcome under treatment.}
#'   \item{psi_0}{Estimated mean outcome under control.}
#' }
do_ipw_nat_inf <- function(
    data, models,
    exclusion_restriction = FALSE,
    S_name, Y_name, Z_name
){
  
  if(!exclusion_restriction){
    # Psi_1 = E[P(S=1 | Z = 0, X) / P(Y = 1 | Z = 0) * E[Y | Z=1, X] ]  
    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
    pi_0_X <- 1 - pi_1_X
    rho_bar_0 <- mean(rho_0_X)
    S <- data[[S_name]]
    Y <- data[[Y_name]]
    Z <- data[[Z_name]]
    
    psi_1 <- mean(
      ( 1 / rho_bar_0 ) * ( Z / pi_1_X ) * 
        ( S + ( rho_0_X - rho_1_X ) * ( 1 - S ) / ( 1 - rho_1_X ) ) * Y
    )
    
    psi_0 <- mean(
      ( S / rho_bar_0 ) * ( ( 1 - Z ) / pi_0_X ) * Y
    )
  }else{
    
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
    pi_0_X <- 1 - pi_1_X
    
    rho_bar_0 <- mean(rho_0_X)
    
    S <- data[[S_name]]
    Y <- data[[Y_name]]
    Z <- data[[Z_name]]
    
    rho_bar_0 <- mean(S[Z == 0])
    pi_Z_X <- ifelse(Z, pi_1_X, pi_0_X)

    psi_1 <- mean( (2*Z - 1) / pi_Z_X * Y ) / mean( ( (1 - Z) / pi_0_X) * S ) + 
      mean(( S / rho_bar_0 ) * ( (1 - Z) / pi_0_X ) * Y)
    
    psi_0 <- mean(
      ( S / rho_bar_0 ) * ( (1 - Z) / pi_0_X ) * Y
    )

  }
  
  growth_effect <- psi_1 - psi_0
  growth_effect_log_mult <- log(psi_1 / psi_0)
  
  out <- c(growth_effect, growth_effect_log_mult, psi_1, psi_0)
  names(out) <- c("additive_effect","log_multiplicative_effect","psi_1","psi_0")
  
  return(out)
}

#' Augmented inverse probability weighted (AIPW) estimator
#'
#' Computes an efficient, doubly robust estimator of the growth effect in the
#' naturally infected principal stratum.
#'
#' @param data A data.frame containing observed data.
#' @param models A named list of fitted nuisance models.
#' @param Y_name Character string specifying the outcome variable. Default is "Y".
#' @param Z_name Character string specifying the treatment variable. Default is "Z".
#' @param S_name Character string specifying the infection indicator. Default is "S".
#' @param X_name Character vector of covariate names.
#' @param exclusion_restriction Logical; whether the exclusion restriction assumption is imposed.
#' @param cross_world Logical; whether a cross-world independence assumption is imposed.
#' @param two_part_model Logical; whether outcome regression uses a two-part model.
#' @param return_se Logical; if TRUE, returns standard errors and influence function outputs.
#'
#' @details
#' This estimator is doubly robust: it is consistent if either the outcome model
#' or the treatment/infection models are correctly specified.
#'
#' @return
#' If \code{return_se = FALSE}, returns:
#' \describe{
#'   \item{additive_effect}{Estimated difference \eqn{\psi_1 - \psi_0}.}
#'   \item{log_multiplicative_effect}{Estimated log ratio.}
#' }
#'
#' If \code{return_se = TRUE}, additionally returns standard errors and estimates
#' of \eqn{\psi_1} and \eqn{\psi_0}. The influence function matrix is stored as an attribute.
#'
#' @export
do_aipw_nat_inf <- function(
  data, models,
  exclusion_restriction = FALSE,
  cross_world = TRUE,
  Y_name = "Y",
  Z_name = "Z",
  S_name = "S",
  X_name = "X",
  return_se = FALSE,
  two_part_model = FALSE
){
  
  if(cross_world & !exclusion_restriction){
    rho_0 <- simple_predict(models$fit_S_Z0_X, newdata = data)
    mu_01 <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    pi_1 <- simple_predict(models$fit_Z_X, newdata = data)
    pi_0 <- 1 - pi_1
    rho_bar_0 <- mean(rho_0)
    
    
    # psi_0 = Weight * E[E[Y | Z = 0, Y = 1, X]]
    
    psi_tilde_0 <- rho_0 / rho_bar_0 * mu_01
    
    psi_0 <- mean( psi_tilde_0 )
    
    augmentation_0 <- (
      (1 - data[[Z_name]]) / pi_0 * ( data[[S_name]] / rho_bar_0 ) * (data[[Y_name]] - mu_01) + 
        (1 - data[[Z_name]]) / pi_0 * ( mu_01 - psi_0 ) / rho_bar_0 * ( data[[S_name]] - rho_0 ) + 
        ( psi_0 / rho_bar_0 ) * ( rho_0 - rho_bar_0 ) + 
        psi_tilde_0 - psi_0
    )
    
    psi_0_aipw <- psi_0 + mean(augmentation_0)
    
    # psi_1 = Weight * E[E[Y | Z = 1, X]]
    mu_11 <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    mu_10 <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    rho_1 <- simple_predict(models$fit_S_Z1_X, newdata = data)
    
    psi_tilde_1 <- rho_1 / rho_bar_0 * mu_11 + ( rho_0 - rho_1 ) / rho_bar_0 * mu_10
    psi_1 <- mean( psi_tilde_1 )
    
    augmentation_1 <- (
      (data[[Z_name]] / pi_1) * (data[[S_name]] / rho_bar_0) * (data[[Y_name]] - mu_11) +
        (data[[Z_name]] / pi_1) * ((1 - data[[S_name]]) / (1 - rho_1)) * (rho_0 - rho_1) / rho_bar_0 * (data[[Y_name]] - mu_10) + 
        (data[[Z_name]] / pi_1) * (mu_11 - mu_10) / rho_bar_0 * (data[[S_name]] - rho_1) + 
        ((1 - data[[Z_name]]) / pi_0) * (mu_10 - psi_1) / rho_bar_0 * (data[[S_name]] - rho_0) - 
        psi_1 / rho_bar_0 * (rho_0 - rho_bar_0) + psi_tilde_1 - psi_1
    )
    
    psi_1_aipw <- psi_1 + mean(augmentation_1)
  }else if(exclusion_restriction & !cross_world){
    
    E_Y_Z0_S1_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
    pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
    pi_0_X <- 1 - pi_1_X
    
    if(!two_part_model){
      df_Z1 <- data.frame(Z = 1, X = data[,colnames(data) %in% X_name, drop = FALSE])
      names(df_Z1) <- c(Z_name, X_name)
      
      df_Z0 <- data.frame(Z = 0, X = data[,colnames(data) %in% X_name, drop = FALSE])
      names(df_Z0) <- c(Z_name, X_name)
      
      E_Y_Z1_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z1)
      E_Y_Z0_X <- simple_predict(models$fit_Y_Z_X, newdata = df_Z0)
    } else{
      # same logic as gcomp
      E_Y_Z0_S0_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
      E_Y_Z1_S0_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
      E_Y_Z1_S1_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
      
      E_Y_Z1_X <- E_Y_Z1_S1_X * rho_1_X + E_Y_Z1_S0_X * (1 - rho_1_X)
      E_Y_Z0_X <- E_Y_Z0_S1_X * rho_0_X + E_Y_Z0_S0_X * (1 - rho_0_X)
      
    }
    
    Z <- data[[Z_name]]
    S <- data[[S_name]]
    Y <- data[[Y_name]]
    E_Y_Z_X <- ifelse(Z, E_Y_Z1_X, E_Y_Z0_X)
    pi_Z_X <- ifelse(Z, pi_1_X, 1 - pi_1_X)
    
    ate <- mean(E_Y_Z1_X - E_Y_Z0_X)
    ate_augmentation <- (2*Z - 1) / pi_Z_X * ( Y - E_Y_Z_X ) + E_Y_Z1_X - E_Y_Z0_X - ate
    ate_aipw <- ate + mean(ate_augmentation)

    rho_bar_0 <- mean(rho_0_X)
    rho_bar_0_augmentation <- (1 - Z) / (1 - pi_1_X) * ( S - rho_0_X ) + rho_0_X - rho_bar_0
    rho_bar_0_aipw <- rho_bar_0 + mean(rho_bar_0_augmentation)

    psi_tilde_0 <- rho_0_X / rho_bar_0 * E_Y_Z0_S1_X
    psi_0 <- mean( psi_tilde_0 )
    psi_0_augmentation <- (
      (1 - Z) / pi_0_X * ( S / rho_bar_0 ) * (Y - E_Y_Z0_S1_X) + 
        (1 - Z) / pi_0_X * ( E_Y_Z0_S1_X - psi_0 ) / rho_bar_0 * ( S - rho_0_X ) + 
        ( psi_0 / rho_bar_0 ) * ( rho_0_X - rho_bar_0 ) + 
        psi_tilde_0 - psi_0
    )
    psi_0_aipw <- psi_0 + mean(psi_0_augmentation)

    psi_1_aipw <- ate_aipw / rho_bar_0_aipw + psi_0_aipw
    
    if_matrix <- cbind(ate_augmentation, rho_bar_0_augmentation, psi_0_augmentation)
    psi_1_gradient <- matrix(c(
      1 / rho_bar_0_aipw, - ate_aipw / rho_bar_0_aipw^2, 1
    ), ncol = 1)
    # augmentation_1 <- c( t(psi_1_gradient) %*% if_matrix )
    # fixed so dims line up??
    augmentation_1 <- if_matrix %*% psi_1_gradient
    augmentation_0 <- psi_0_augmentation
  }else if(exclusion_restriction & cross_world){
    rho_0_X <- simple_predict(models$fit_S_Z0_X, newdata = data)
    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
    
    mu_01_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
    pi_0_X <- 1 - pi_1_X
    rho_bar_0 <- mean(rho_0_X)
    
    # psi_0 = Weight * E[E[Y | Z = 0, Y = 1, X]]
    psi_tilde_0_X <- rho_0_X / rho_bar_0 * mu_01_X
    
    psi_0 <- mean( psi_tilde_0_X )
    
    augmentation_0 <- (
      (1 - data[[Z_name]]) / pi_0_X * ( data[[S_name]] / rho_bar_0 ) * (data[[Y_name]] - mu_01_X) + 
        (1 - data[[Z_name]]) / pi_0_X * ( mu_01_X - psi_0 ) / rho_bar_0 * ( data[[S_name]] - rho_0_X ) + 
        ( psi_0 / rho_bar_0 ) * ( rho_0_X - rho_bar_0 ) + 
        psi_tilde_0_X - psi_0
    )
    
    psi_0_aipw <- psi_0 + mean(augmentation_0)
    
    mu_11_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    mu_10_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    mu_00_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
    
    # original
    # mu_dot0_X <- pi_1_X * mu_10_X + pi_0_X * mu_00_X
    
    # more generalizable
    # wt_1
    pi_1_S0_X <- (1 - rho_1_X) * pi_1_X /
      (((1 - rho_0_X)*pi_0_X) +
      (1 - rho_1_X) * pi_1_X)
    
    # wt_0
    pi_0_S0_X <- (1 - rho_0_X) * pi_0_X /
      (((1 - rho_0_X)*pi_0_X) +
      (1 - rho_1_X) * pi_1_X)
    
    mu_dot0_X <- pi_1_S0_X * mu_10_X + pi_0_S0_X * mu_00_X

    rho_bar_dot <- pi_1_X * rho_1_X + pi_0_X * rho_0_X

    psi_tilde_1_X <- rho_1_X / rho_bar_0 * mu_11_X + ( rho_0_X - rho_1_X ) / rho_bar_0 * mu_dot0_X
    psi_1 <- mean( psi_tilde_1_X )
    
    augmentation_1 <- (
      (data[[Z_name]] / pi_1_X) * (data[[S_name]] / rho_bar_0) * (data[[Y_name]] - mu_11_X) +
        ((1 - data[[S_name]]) / (1 - rho_bar_dot)) * (rho_0_X - rho_1_X) / rho_bar_0 * (data[[Y_name]] - mu_dot0_X) + 
        (data[[Z_name]] / pi_1_X) * (mu_11_X - mu_dot0_X) / rho_bar_0 * (data[[S_name]] - rho_1_X) + 
        ((1 - data[[Z_name]]) / pi_0_X) * (mu_dot0_X - psi_1) / rho_bar_0 * (data[[S_name]] - rho_0_X) - 
        psi_1 / rho_bar_0 * (rho_0_X - rho_bar_0) + psi_tilde_1_X - psi_1
    )
    
    psi_1_aipw <- psi_1 + mean(augmentation_1)
  }else{
    stop("Must assume exclusion_restriction, cross_world, or both.")
  }
  
  # Additive effect
  efficient_growth_effect <- psi_1_aipw - psi_0_aipw
  se <- sqrt(stats::var(augmentation_1 - augmentation_0) / dim(data)[1])
  
  se_psi_1 <- sqrt(stats::var(augmentation_1) / dim(data)[1])
  se_psi_0 <- sqrt(stats::var(augmentation_0) / dim(data)[1])
  
  # Multiplicative effect (log scale)
  efficient_growth_effect_log_mult <- log(psi_1_aipw / psi_0_aipw)
  
  # Get SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  colnames(if_matrix) <- c("augmentation_1", "augmentation_2")
  
  cov_matrix <- stats::cov(if_matrix) / dim(data)[1]

  gradient <- matrix(c(1 / psi_1_aipw, -1 / psi_0_aipw), ncol = 1)
  
  se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
  
  if(return_se){
    out <- c(efficient_growth_effect, se, efficient_growth_effect_log_mult, se_log_mult_eff, psi_1_aipw, se_psi_1, psi_0_aipw, se_psi_0)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se", "psi_1", "se_psi_1", "psi_0", "se_psi_0")
    
    # added to return influence fn matrix for vaccine trial project without having to restructure whole package rn. eventually change return type to list
    attr(out, "if_matrix") <- if_matrix
    
    return(out)
  }else{
    out <- c(efficient_growth_effect, efficient_growth_effect_log_mult)
    names(out) <- c("additive_effect", "log_multiplicative_effect")
    return(out)
  }
}

#' Targeted maximum likelihood estimator (TMLE)
#'
#' Computes a targeted maximum likelihood estimate of the growth effect in the
#' naturally infected principal stratum.
#'
#' @param data A data.frame containing observed data.
#' @param models A named list of fitted nuisance models.
#' @param exclusion_restriction Logical; currently not supported.
#' @param Y_name Character string specifying the outcome variable.
#' @param Z_name Character string specifying the treatment variable.
#' @param S_name Character string specifying the infection indicator.
#' @param return_se Logical; whether to return standard errors.
#' @param max_iter Maximum number of targeting iterations.
#' @param tol Convergence tolerance for the efficient influence function.
#'
#' @details
#' TMLE updates initial nuisance estimates via iterative targeting steps until
#' the empirical mean of the efficient influence function is approximately zero.
#'
#' @return
#' A named numeric vector containing the additive and log multiplicative effects.
#' If \code{return_se = TRUE}, standard errors and component estimates are included.
do_tmle_nat_inf <- function(
    data, models, 
    exclusion_restriction = FALSE,
    Y_name = "Y", Z_name = "Z", S_name = "S",
    return_se = FALSE, max_iter = 10,
    tol = 1 / (sqrt(dim(data)[1]) * log(dim(data)[1]))
){
  
  if(exclusion_restriction){
    stop("TMLE with exclusion restriction not implemented yet")
    }else{
    idx_Z0 <- which(data[[Z_name]] == 0)
    idx_Z1 <- which(data[[Z_name]] == 1)
    idx_Z0_S1 <- which(data[[Z_name]] == 0 & data[[S_name]] == 1)
    l <- min(data[[Y_name]])
    u <- max(data[[Y_name]])
    
    
    pi_1 <- simple_predict(models$fit_Z_X, newdata = data)
    
    rho_0 <- simple_predict(models$fit_S_Z0_X, newdata = data)
    rho_1 <- simple_predict(models$fit_S_Z1_X, newdata = data)
    
    mu_11 <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    mu_10 <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    mu_01 <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    
    pi_0 <- 1 - pi_1
    rho_bar_0 <- mean(rho_0)
    
    psi_tilde_0 <- rho_0 / rho_bar_0 * mu_01
    psi_0 <- mean( psi_tilde_0 )
    
    psi_tilde_1 <- rho_1 / rho_bar_0 * mu_11 + ( rho_0 - rho_1 ) / rho_bar_0 * mu_10
    psi_1 <- mean( psi_tilde_1 )
    
    phi_0 <- function(data, Z_name, S_name, Y_name, pi_0, rho_0, rho_bar_0, mu_01, psi_tilde_0, psi_0) {
      (
        (1 - data[[Z_name]]) / pi_0 * (data[[S_name]] / rho_bar_0) * (data[[Y_name]] - mu_01) +
          (1 - data[[Z_name]]) / pi_0 * (mu_01 - psi_0) / rho_bar_0 * (data[[S_name]] - rho_0) +
          (psi_0 / rho_bar_0) * (rho_0 - rho_bar_0) +
          psi_tilde_0 - psi_0
      )
    }
    
    phi_1 <- function(data, Z_name, S_name, Y_name, pi_1, pi_0, rho_0, rho_bar_0, rho_1, mu_11, mu_10, psi_tilde_1, psi_1) {
      (
        (data[[Z_name]] / pi_1) * (data[[S_name]] / rho_bar_0) * (data[[Y_name]] - mu_11) +
          (data[[Z_name]] / pi_1) * ((1 - data[[S_name]]) / (1 - rho_1)) * (rho_0 - rho_1) / rho_bar_0 * (data[[Y_name]] - mu_10) +
          (data[[Z_name]] / pi_1) * (mu_11 - mu_10) / rho_bar_0 * (data[[S_name]] - rho_1) +
          ((1 - data[[Z_name]]) / pi_0) * (mu_10 - psi_1) / rho_bar_0 * (data[[S_name]] - rho_0) -
          psi_1 / rho_bar_0 * (rho_0 - rho_bar_0) +
          psi_tilde_1 - psi_1
      )
    }
    
    trim_logit <- function(p, tol = 1e-3){ 
      p[p < tol] <- tol
      p[p > 1 - tol] <- 1 - tol
      return(qlogis(p))
    }
    
    scale_01 <- function(x, l, u){
      ( x - l ) / ( u - l )
    }
    
    rescale_01 <- function(x, l, u){
      x * (u - l) + l
    }
    
    phi_0_data <- phi_0(data, Z_name, S_name, Y_name, pi_0, rho_0, rho_bar_0, mu_01, psi_tilde_0, psi_0)
    phi_1_data <- phi_1(data, Z_name, S_name, Y_name, pi_1, pi_0, rho_0, rho_bar_0, rho_1, mu_11, mu_10, psi_tilde_1, psi_1)
    phi_ge_data <- phi_1_data - phi_0_data
    
    mean_phi_0 <- mean(phi_0_data)
    mean_phi_1 <- mean(phi_1_data)
    mean_phi_ge_data <- mean(phi_ge_data)
    
    iter <- 0
    
    mu_11_star <- mu_11
    mu_10_star <- mu_10
    mu_01_star <- mu_01
    rho_1_star <- rho_1
    rho_0_star <- rho_0
    rho_bar_0_star <- rho_bar_0
    psi_tilde_0_star <- psi_tilde_0
    psi_tilde_1_star <- psi_tilde_1
    psi_0_star <- psi_0
    psi_1_star <- psi_1
    
    while((mean_phi_0^2 + mean_phi_1^2 + mean_phi_ge_data^2) > tol & iter <= max_iter){
      # cat("iter", iter, "\n")
      # cat("mean_eif", mean_phi_ge_data, "\n")
      
      # target mu's
      Y_scale <- scale_01(data[[Y_name]], l, u)
      
      # target mu_11
      mu_11_star_scale <- scale_01(mu_11_star, l, u)
      logit_mu_11_star_scale <- trim_logit(mu_11_star_scale)
      target_wt <- (
        (data[[Z_name]] / pi_1) * (data[[S_name]] / rho_bar_0_star)
      )
      
      target_data <- data.frame(
        Y_scale = Y_scale,
        target_wt = target_wt,
        logit_mu_11_star_scale = logit_mu_11_star_scale
      )
      target_fit <- suppressWarnings(stats::glm(
        Y_scale ~ offset(logit_mu_11_star_scale), 
        weight = target_wt,
        family = binomial(),
        data = target_data,
        start = c(0)
      ))
      mu_11_star <- rescale_01(target_fit$fitted.values, l, u)
      
      # target mu_01
      mu_01_star_scale <- scale_01(mu_01_star, l, u)
      logit_mu_01_star_scale <- trim_logit(mu_01_star_scale)
      target_wt <- (
        ((1 - data[[Z_name]]) / (1 - pi_1)) * (data[[S_name]] / rho_bar_0_star)
      )
      
      target_data <- data.frame(
        Y_scale = Y_scale,
        target_wt = target_wt,
        logit_mu_01_star_scale = logit_mu_01_star_scale
      )
      target_fit <- suppressWarnings(stats::glm(
        Y_scale ~ offset(logit_mu_01_star_scale), 
        weight = target_wt,
        family = binomial(),
        data = target_data,
        start = c(0)
      ))
      mu_01_star <- rescale_01(target_fit$fitted.values, l, u)
      
      # target mu_10
      mu_10_star_scale <- scale_01(mu_10_star, l, u)
      logit_mu_10_star_scale <- trim_logit(mu_10_star_scale)
      target_wt <- with(data, 
                        ( data[[Z_name]] / pi_1 ) * ( (1 - data[[S_name]]) / rho_bar_0_star ) 
      )
      H1 <- ( rho_0_star - rho_1_star ) / ( 1 - rho_1_star )
      target_data <- data.frame(
        Y_scale = Y_scale,
        target_wt = target_wt,
        H1 = H1,
        logit_mu_10_star_scale = logit_mu_10_star_scale
      )
      target_fit <- suppressWarnings(stats::glm(
        Y_scale ~ -1 + offset(logit_mu_10_star_scale) + H1, 
        weight = target_wt,
        family = binomial(),
        data = target_data,
        start = c(0)
      ))
      mu_10_star <- rescale_01(target_fit$fitted.values, l, u)
      
      psi_tilde_1_star <- rho_1_star / rho_bar_0_star * mu_11_star + ( rho_0_star - rho_1_star ) / rho_bar_0_star * mu_10_star
      psi_1_star <- mean( psi_tilde_1_star )
      
      psi_tilde_0_star <- rho_0_star / rho_bar_0_star * mu_01_star
      psi_0_star <- mean( psi_tilde_0_star )
      
      # target rho_0
      H1 <- mu_10_star - psi_1_star
      H0 <- mu_01_star - psi_0_star
      logit_rho_0_star <- trim_logit(rho_0_star)
      target_wt <- (1 - data[[Z_name]]) / pi_0
      
      # with linear models, these may be perfectly correlated, but numerically
      # R thinks they are not and tries to fit a glm, which blows up. setting 
      # H0 to a constant in these cases will remove the term from the model because
      # the model also includes an intercept
      if(stats::sd(H1) > 0 & stats::sd(H0) > 0){
        if(stats::cor(H1, H0) > 0.99999){
          H0 <- 1
        }
      
        target_data <- data.frame(
          S_inf = data[[S_name]],
          target_wt = target_wt,
          H1 = H1,
          H0 = H0,
          logit_rho_0_star = logit_rho_0_star
        )
        target_data <- stats::setNames(target_data, c(S_name, names(target_data[-1])))
        
        # include intercept so rho_bar_0_star is still mean(Y[Z == 0])
        target_fit <- stats::glm(
          as.formula(paste0(S_name," ~ offset(logit_rho_0_star) + H1 + H0")), 
          family = binomial(),
          weight = target_wt,
          data = target_data,
          start = c(0, 0, 0)
        )
        rho_0_star <- target_fit$fitted.values
      }
      
      # shouldn't change because of intercept, but just in case
      rho_bar_0_star <- mean(rho_0_star)
      
      ## sanity check
      # tmp <- with(data, 
      # ( (1 - Z) / pi_0 ) * ( mu_01_star - psi_0_star ) / rho_bar_0_star * ( S_inf - rho_0_star )
      # )
      # mean(tmp) # should be small
      # tmp <- with(data, 
      #     ( (1 - Z) / pi_0 ) * ( (mu_10_star - psi_1_star) / rho_bar_0_star ) * (S_inf - rho_0_star) 
      # )
      # mean(tmp) # should be small
      
      # target rho_1
      H1 <- mu_11_star - mu_10_star
      logit_rho_1_star <- trim_logit(rho_1_star)
      target_wt <- data[[Z_name]] / pi_1
      target_data <- data.frame(
        S_name = data[[S_name]],
        target_wt = target_wt,
        H1 = H1,
        logit_rho_1_star = logit_rho_1_star
      )
      target_data <- stats::setNames(target_data, c(S_name, names(target_data[-1])))
      
      # include intercept so rho_bar_0_star is still mean(Y[Z == 0])
      target_fit <- stats::glm(
        as.formula(paste0(S_name, " ~ -1 + offset(logit_rho_1_star) + H1")), 
        family = binomial(),
        weight = target_wt,
        data = target_data,
        start = c(0)
      )
      rho_1_star <- target_fit$fitted.values
      
      psi_tilde_1_star <- rho_1_star / rho_bar_0_star * mu_11_star + ( rho_0_star - rho_1_star ) / rho_bar_0_star * mu_10_star
      psi_1_star <- mean( psi_tilde_1_star )
      
      psi_tilde_0_star <- rho_0_star / rho_bar_0_star * mu_01_star
      psi_0_star <- mean( psi_tilde_0_star )
      
      phi_0_data <- phi_0(data, Z_name, S_name, Y_name, pi_0, rho_0_star, rho_bar_0_star, mu_01_star, psi_tilde_0_star, psi_0_star)
      phi_1_data <- phi_1(data, Z_name, S_name, Y_name, pi_1, pi_0, rho_0_star, rho_bar_0_star, rho_1_star, mu_11_star, mu_10_star, psi_tilde_1_star, psi_1_star)
      phi_ge_data <- phi_1_data - phi_0_data
      
      mean_phi_0 <- mean(phi_0_data)
      mean_phi_1 <- mean(phi_1_data)
      mean_phi_ge_data <- mean(phi_ge_data)
      
      iter <- iter + 1
    }
  }
  
  # Additive growth effect
  tmle_ge <- psi_1_star - psi_0_star
  
  # Log multiplicative growth effect
  tmle_ge_log_mult <- log(psi_1_star / psi_0_star)
  
  if(return_se){
    se <- sqrt(stats::var(phi_ge_data) / dim(data)[1])
    
    if_matrix <- cbind(phi_1_data, phi_0_data)
    cov_matrix <- stats::cov(if_matrix) / dim(data)[1]
    # 1/psi_1, -1/psi_0
    gradient <- matrix(c(1 / psi_1_star, -1 / psi_0_star), ncol = 1)
    se_log_mult_eff <- sqrt(t(gradient) %*% cov_matrix %*% gradient)
    
    se_psi_1 <- sqrt(diag(cov_matrix))[1]
    se_psi_0 <- sqrt(diag(cov_matrix))[2]
    
    out <- c(tmle_ge, se, tmle_ge_log_mult, se_log_mult_eff, psi_1_star, se_psi_1, psi_0_star, se_psi_0)
    names(out) <- c("additive_effect", "additive_se", "log_multiplicative_effect", "log_multiplicative_se", "psi_1", "se_psi_1", "psi_0", "se_psi_0")
    return(out)
  }else{
    out <- c(tmle_ge, tmle_ge_log_mult)
    names(out) <- c("additive_effect", "log_multiplicative_effect")
    return(out)
  }
  
}

#' Sensitivity analysis using AIPW estimator
#'
#' Performs sensitivity analysis for violations of cross-world assumptions using
#' a parameter \eqn{\epsilon}.
#'
#' @param data A data.frame containing observed data.
#' @param models A named list of fitted nuisance models.
#' @param Y_name Character string specifying the outcome variable.
#' @param Z_name Character string specifying the treatment variable.
#' @param S_name Character string specifying the infection indicator.
#' @param epsilon Numeric vector of sensitivity parameters.
#' @param return_se Logical; whether to return standard errors.
#'
#' @return A list with class \code{"sens_cw"} containing:
#' \describe{
#'   \item{epsilon}{Grid of sensitivity parameter values.}
#'   \item{psi_1_epsilon}{Estimated \eqn{\psi_1} for each epsilon.}
#'   \item{psi_0_aipw}{Baseline estimate of \eqn{\psi_0}.}
#'   \item{additive_effect}{Additive effects across epsilon values.}
#'   \item{log_multiplicative_effect}{Log multiplicative effects across epsilon values.}
#' }
#' 
#' @export
do_sens_cw_aipw_nat_inf <- function(data,
                                   models,
                                   Y_name = "Y",
                                   Z_name = "Z",
                                   S_name = "S",
                                   epsilon = exp(seq(log(0.5), log(2), length = 49)),
                                   return_se = FALSE){
  
  
  
  pi_1 <- simple_predict(models$fit_Z_X, newdata = data)
  
  # vaccine probabilities
  # pi_1 <- models$fit_Z_X$fitted.values
  pi_0 <- 1 - pi_1
  
  # Get weight
  sub_Z0 <- data[data[[Z_name]] == 0,]
  
  # rho_bar_0 <- mean(sub_Z0[[S_name]])
  
  data_0 <- data; data_0[[Z_name]] <- 0
  data_1 <- data; data_1[[Z_name]] <- 1
  
  rho_0_X <- simple_predict(models$fit_S_Z_X, newdata = data_0)
  mu_01_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
  
  rho_bar_0 <- mean(rho_0_X)
  
  psi_tilde_0_X <- rho_0_X / rho_bar_0 * mu_01_X
  
  psi_0 <- mean( psi_tilde_0_X )
  
  Z_i <- data[[Z_name]]
  S_i <- data[[S_name]]
  Y_i <- data[[Y_name]]
  
  augmentation_0 <- (
    (1 - Z_i) / pi_0 * ( S_i / rho_bar_0 ) * (Y_i - mu_01_X) + 
      (1 - Z_i) / pi_0 * ( mu_01_X - psi_0 ) / rho_bar_0 * ( S_i - rho_0_X ) + 
      ( psi_0 / rho_bar_0 ) * ( rho_0_X - rho_bar_0 ) + 
      psi_tilde_0_X - psi_0
  )
  
  psi_0_aipw <- psi_0 + mean(augmentation_0)
  
  if(inherits(models$fit_Y_Z1_S1_X, "SuperLearner")){
    mu_11_X <- predict(models$fit_Y_Z1_S1_X, newdata = data, type = "response")$pred
    mu_10_X <- predict(models$fit_Y_Z1_S0_X, newdata = data, type = "response")$pred
    rho_1_X <- predict(models$fit_S_Z_X, newdata = data_1, type = "response")$pred
  } else{
    mu_11_X <- stats::predict(models$fit_Y_Z1_S1_X, newdata = data, type = "response")
    mu_10_X <- stats::predict(models$fit_Y_Z1_S0_X, newdata = data, type = "response")
    rho_1_X <- stats::predict(models$fit_S_Z_X, newdata = data_1, type = "response")
  }
  
  psi_11_epsilon_X <- rho_1_X / rho_bar_0 * mu_11_X 
  psi_10_epsilon_X <- sapply(epsilon, function(eps){
    (rho_0_X - rho_1_X) / rho_bar_0 * (1 - rho_1_X) / ((1 - eps) * rho_0_X - rho_1_X + eps) * mu_10_X
  }, simplify = FALSE)
  psi_11_epsilon <- mean(psi_11_epsilon_X)
  psi_10_epsilon <- lapply(psi_10_epsilon_X, mean)
  
  psi_1_epsilon <- lapply(psi_10_epsilon, function(psi_10_eps){
    psi_11_epsilon + psi_10_eps
  })
  
  augmentation_1_epsilon <- mapply(
  eps = epsilon, psi_10_eps_X = psi_10_epsilon_X, psi_10_eps = psi_10_epsilon,
  function(eps, psi_10_eps_X, psi_10_eps){
    ( Z_i / pi_1) * ( S_i / rho_bar_0 ) * ( Y_i - mu_11_X ) + 
      Z_i / pi_1 * ( mu_11_X / rho_bar_0 ) * ( S_i - rho_1_X ) - 
      psi_11_epsilon / rho_bar_0 * ( 1 - Z_i ) / pi_0 * (S_i - rho_0_X) -
      psi_11_epsilon / rho_bar_0 * (rho_0_X - rho_bar_0) +
      psi_11_epsilon_X - psi_11_epsilon + 
      Z_i / pi_1 * (1 - S_i) / (rho_bar_0) * (rho_0_X - rho_1_X) / ((1 - eps) * rho_0_X - rho_1_X + eps) * ( Y_i - mu_10_X ) + 
      ( 1 - Z_i ) / pi_0 * (1 - rho_1_X) / ((1 - eps) * rho_0_X - rho_1_X + eps) * mu_10_X / rho_bar_0 * ( S_i - rho_0_X ) -
      Z_i / pi_1 * (1 - rho_1_X) / ((1 - eps) * rho_0_X - rho_1_X + eps) * mu_10_X / rho_bar_0 * ( S_i - rho_1_X ) -
      psi_10_eps / rho_bar_0 * (1 - Z_i) / pi_0 * ( S_i - rho_0_X ) -
      psi_10_eps / rho_bar_0 * (rho_0_X - rho_bar_0) -
      Z_i / pi_1 * ( rho_0_X - rho_1_X ) / rho_bar_0 * mu_10_X / ((1 - eps) * rho_0_X - rho_1_X + eps) * ( S_i - rho_1_X ) - 
      (1 - eps) * (1 - Z_i) / (pi_0) * (rho_0_X - rho_1_X) / rho_bar_0 * (1 - rho_1_X) / ((1 - eps) * rho_0_X - rho_1_X + eps)^2 * mu_10_X * (S_i - rho_0_X) + 
      Z_i / pi_1 * (rho_0_X - rho_1_X) / rho_bar_0 * (1 - rho_1_X) / ((1 - eps) * rho_0_X - rho_1_X + eps)^2 * mu_10_X * ( S_i - rho_1_X ) + 
      psi_10_eps_X - psi_10_eps
  }, SIMPLIFY = FALSE
  )
  
  psi_1_epsilon_aipw <- mapply(
    psi_1_eps = psi_1_epsilon, augmentation_1_eps = augmentation_1_epsilon, 
    function(psi_1_eps, augmentation_1_eps){
      psi_1_eps + mean(augmentation_1_eps)
    },
    SIMPLIFY = FALSE
  )
  
  # Additive effect
  efficient_growth_effect_epsilon <- lapply(psi_1_epsilon_aipw, function(psi_1_eps){
    psi_1_eps - psi_0_aipw
  })
  se_epsilon <- lapply(augmentation_1_epsilon, function(augmentation_1_eps){
    sqrt(stats::var(augmentation_1_eps - augmentation_0) / dim(data)[1])
  })
  
  # Multiplicative effect (log scale)
  efficient_growth_effect_log_mult_epsilon <- lapply(psi_1_epsilon_aipw, function(psi_1_eps){
    log(psi_1_eps / psi_0_aipw)
  })
  
  cov_matrices <- mapply(
    augmentation_1_eps = augmentation_1_epsilon, 
    function(augmentation_1_eps){
      if_matrix <- cbind(augmentation_1_eps, augmentation_0)
      cov_matrix <- stats::cov(if_matrix) / dim(data)[1]
      return(cov_matrix)
    }, SIMPLIFY =  FALSE
  )
  
  # Get SE using IF matrix same way as TMLE
  se_log_mult_eff <- mapply(
    augmentation_1_eps = augmentation_1_epsilon, 
    psi_1_eps_aipw = psi_1_epsilon_aipw,
    function(augmentation_1_eps, psi_1_eps_aipw){
      if_matrix <- cbind(augmentation_1_eps, augmentation_0)
      cov_matrix <- stats::cov(if_matrix) / dim(data)[1]
      gradient <- matrix(c(1 / psi_1_eps_aipw, -1 / psi_0_aipw), ncol = 1)
      return(sqrt(t(gradient) %*% cov_matrix %*% gradient))
    })

  if(return_se){
    out <- list(
      epsilon = epsilon,
      psi_1_epsilon = unlist(psi_1_epsilon_aipw, use.names = FALSE),
      psi_0_aipw = psi_0_aipw,
      additive_effect = unlist(efficient_growth_effect_epsilon, use.names = FALSE), 
      additive_se = unlist(se_epsilon, use.names = FALSE), 
      log_multiplicative_effect = unlist(efficient_growth_effect_log_mult_epsilon, use.names = FALSE), 
      log_multiplicative_se = unlist(se_log_mult_eff, use.names = FALSE),
      cov_matrices = cov_matrices
    )
  }else{
    out <- list(
      epsilon = epsilon,
      psi_1_epsilon = unlist(psi_1_epsilon_aipw, use.names = FALSE),
      psi_0_aipw = psi_0_aipw,
      additive_effect = unlist(efficient_growth_effect_epsilon, use.names = FALSE), 
      log_multiplicative_effect = unlist(efficient_growth_effect_log_mult_epsilon, use.names = FALSE)
    )
  }
  
  class(out) <- "sens_cw"
  
  return(out)
}

#' Sensitivity analysis to monotonicity using AIPW estimator
#'
#' Performs sensitivity analysis for violations of monotonicity using
#' a parameter \eqn{\epsilon}.
#'
#' @param data A data.frame containing observed data.
#' @param models A named list of fitted nuisance models.
#' @param Y_name Character string specifying the outcome variable.
#' @param Z_name Character string specifying the treatment variable.
#' @param S_name Character string specifying the infection indicator.
#' @param epsilon Numeric vector of sensitivity parameters.
#'
#' @return A list with class \code{"monotonicity_sens"} containing:
#' \describe{
#'   
#' }
#' 
#' @export
do_sens_mono_nat_inf <- function(data,
                                 models,
                                 Y_name = "Y",
                                 Z_name = "Z",
                                 S_name = "S",
                                 epsilon = exp(seq(log(0.5), log(2), length = 49))
                                 ){
  
  # ---------------------------------------
  # Get bounds on p
  # ---------------------------------------
  
  # 1. fit your logistic regression of S ~ Z + X
  # 2. get predictions for all X setting Z = 1 => these are estimates \rho_{1,n}(X_i) of \rho_1(X_i), i = 1, \dots, n
  # 3. get predictions for all X setting Z = 0 => these are estimates \rho_{1,0}(X_i) of \rho_0(X_i), i = 1, \dots, n
  # 4. calculate max(0, max_{i=1,...n} (\rho_{1,n}(X_i) - \rho_{0,n}(X_i))) -> call this p_min
  # 5. calculate min( min_{i = 1, ..., n} \rho_{1,n}(X_i), min_{i=1,...,n} (1 - \rho_{0,n}(X_i)) ) -> call this p_max
  # 6. let p vary between p_min and p_max
  
  data_0 <- data; data_0[[Z_name]] <- 0
  data_1 <- data; data_1[[Z_name]] <- 1
  
  # 2. get predictions for all X setting Z = 1 => these are estimates \rho_{1,n}(X_i) of \rho_1(X_i), i = 1, \dots, n
  rho_0_X <- simple_predict(models$fit_S_Z_X, newdata = data_0)
  
  # 3. get predictions for all X setting Z = 0 => these are estimates \rho_{1,0}(X_i) of \rho_0(X_i), i = 1, \dots, n
  rho_1_X <- simple_predict(models$fit_S_Z_X, newdata = data_1)
  
  # 4. calculate max(0, max_{i=1,...n} (\rho_{1,n}(X_i) - \rho_{0,n}(X_i))) -> call this p_min
  rho_dif <- rho_1_X - rho_0_X
  p_min <- max(0, rho_dif)
  
  # 5. calculate min( min_{i = 1, ..., n} \rho_{1,n}(X_i), min_{i=1,...,n} (1 - \rho_{0,n}(X_i)) ) -> call this p_max
  p_max <- min(min(rho_1_X), min(1 - rho_0_X))
  
  # 6. let p vary between p_min and p_max, equally spaced on log scale but if starting at 0 hard code
  if(p_min > 0){
    p_range <- exp(seq(log(p_min), log(p_max), length = length(epsilon)))
  } else{
    p_range <- c(
      0,
      exp(seq(log(1e-6), log(p_max), length.out = length(epsilon)))[-1]
    )
  }

  # grid of epsilons and ps
  p_epsilon_grid <- expand.grid(p = p_range, epsilon = epsilon)
  
  # ---------------------------------------
  # Plug-in estimate
  # ---------------------------------------
  
  # addditional components
  mu_11_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
  mu_10_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
  mu_01_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
  mu_00_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
  
  pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
  pi_0_X <- 1 - pi_1_X
  
  # wt_1
  pi_1_S0_X <- (1 - rho_1_X) * pi_1_X /
    (((1 - rho_0_X)*pi_0_X) +
    (1 - rho_1_X) * pi_1_X)
  
  # wt_0
  pi_0_S0_X <- (1 - rho_0_X) * pi_0_X /
    (((1 - rho_0_X)*pi_0_X) +
    (1 - rho_1_X) * pi_1_X)
  
  mu_dot0_X <- pi_1_S0_X * mu_10_X + pi_0_S0_X * mu_00_X
  
  rho_bar_0 <- mean(rho_0_X)
  
  # put together to estimate \psi_1 = E[Y(1) | S(0) = 1] ---------------------
  # outer weight
  wt_X <- rho_0_X / rho_bar_0
  
  psi_1_epsilon_p_X <- lapply(seq_len(nrow(p_epsilon_grid)), function(i) {
    
    p_i <- p_epsilon_grid$p[i]
    epsilon_i <- p_epsilon_grid$epsilon[i]
    
    term1_X <- mu_11_X /
      (((rho_1_X - p_i) / rho_1_X) + (p_i / (epsilon_i * rho_1_X))) *
      ((rho_1_X - p_i) / rho_0_X)
    
    term2_X <- mu_dot0_X *
      ((rho_0_X - rho_1_X + p_i) / rho_0_X)
    
    wt_X * (term1_X + term2_X)
  })
  
  psi_1_epsilon_p <- sapply(psi_1_epsilon_p_X, mean)
  
  
  # psi_0
  psi_tilde_0_X <- wt_X * mu_01_X
  psi_0 <- mean(psi_tilde_0_X)
  
  # effects
  additive_growth_effect <- psi_1_epsilon_p - psi_0
  log_mult_growth_effect <- log(psi_1_epsilon_p / psi_0)
  
  out <- list(
    epsilon = p_epsilon_grid$epsilon,
    p = p_epsilon_grid$p,
    grid = p_epsilon_grid,
    p_min = p_min,
    p_max = p_max,
    psi_1_epsilon_p = psi_1_epsilon_p,
    psi_0 = psi_0,
    additive_effect = additive_growth_effect,
    log_multiplicative_effect = log_mult_growth_effect
  )
  
  class(out) <- "sens_mono"
  
  return(out)
  
}

#' Nonparametric bounds for naturally infected principal stratum
#' 
#' @param data A data.frame containing observed data.
#' @param Y_name Outcome variable name.
#' @param Z_name Treatment variable name.
#' @param S_name Infection indicator name.
#' @param family Outcome type: "gaussian" (continuous) or "binomial" (binary).
#' 
#' @return A named numeric vector with:
#' \describe{
#'   \item{E_Y0__S0_1}{Observed mean outcome among unvaccinated infected.}
#'   \item{E_Y1__S0_1_lower}{Lower bound.}
#'   \item{E_Y1__S0_1_upper}{Upper bound.}
#'   \item{additive_effect_lower}{Lower bound on additive effect.}
#'   \item{additive_effect_upper}{Upper bound on additive effect.}
#'   \item{mult_effect_lower}{Lower bound on multiplicative effect.}
#'   \item{mult_effect_upper}{Upper bound on multiplicative effect.}
#' }
#'
#' @export
get_bound_nat_inf <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    family = "gaussian"
){
  
  # Step 1: rhobar_z_n
  
  # 1.1 rhobar_0_n (or mean in subset)
  rhobar_0_n <- mean(data[[S_name]][data[[Z_name]] == 0])
  
  # 1.2 rhobar_1_n
  rhobar_1_n <- mean(data[[S_name]][data[[Z_name]] == 1])
  
  if(rhobar_0_n > rhobar_1_n){
    # Step 2: mubar_11_n 
    mubar_11_n_num <- sum(data[[Y_name]]*data[[S_name]]*data[[Z_name]])
    mubar_11_n_denom <- sum(data[[S_name]]*data[[Z_name]])
    
    if(mubar_11_n_denom == 0){ # denominator NA in some bootstrap replicates, change to 0
      mubar_11_n <- 0
    } else{
      mubar_11_n <- mubar_11_n_num / mubar_11_n_denom
    }
    
    # Step 3: q_n (relative size of protected in (immune + protected) in vax)
    q_n = 1 - (1 - rhobar_0_n) / (1 - rhobar_1_n)
    
    # Step 4: q_n^th quintiles of S__Z1_S0 (aka Y__Z1_S0, need to rename everything at some point)
    Y__Z1_S0 <- data[[Y_name]][which(data[[Z_name]] == 1 & data[[S_name]] == 0)]
    q_nth_quintile <- stats::quantile(Y__Z1_S0, probs = q_n) # NOTE failing here if condition on line 709 not met
    one_minus_q_nth_quintile <- stats::quantile(Y__Z1_S0, probs = 1 - q_n)
    
    # Step 5: mubar_10_l,u_n 
    if(family == "gaussian"){
      mubar_10_l_n <- sum(data[[Y_name]] * as.numeric(data[[S_name]] == 0 & data[[Z_name]] == 1 & data[[Y_name]] < q_nth_quintile )) / 
        sum(as.numeric(data[[S_name]] == 0 & data[[Z_name]] == 1 & data[[Y_name]] < q_nth_quintile ))
      
      mubar_10_u_n <- sum(data[[Y_name]] * as.numeric(data[[S_name]] == 0 & data[[Z_name]] == 1 & data[[Y_name]] > one_minus_q_nth_quintile )) / 
        sum(as.numeric(data[[S_name]] == 0 & data[[Z_name]] == 1 & data[[Y_name]] > one_minus_q_nth_quintile ))
    } else{
      # Binary outcome
      
      # Vaccinated, Uninfected
      data__Z1_S0 <- data[which(data[[Z_name]] == 1 & data[[S_name]] == 0),]
      
      target_num <- ceiling(q_n * nrow(data__Z1_S0))
      num_0s <- length(which(data__Z1_S0[[Y_name]] == 0))
      num_1s <- length(which(data__Z1_S0[[Y_name]] == 1))
      
      ## Lower Bound:
      
      # Check if at least q_n * 100 % 0s in the vax uninfected 
      if(num_0s >= target_num){
        # If so, muhat_10_l = 0
        mubar_10_l_n <- 0
      } else{
        # Else, mubar_10_l = (q_n - prop zeros in vax uninf) / q_n
        mubar_10_l_n <- (q_n - (num_0s / nrow(data__Z1_S0)) ) /
          q_n
      }
      
      ## Upper Bound:
      
      # Check if at least q_n * 100 % 1s in the vax uninfected 
      if(num_1s >= target_num){
        # If so, mubar_10_u = 1
        mubar_10_u_n <- 1
      } else{
        # Else, mubar_10_u = (q_n - prop ones in vax uninf) / q_n
        
        mubar_10_u_n <- 1 - ((q_n - (num_1s / nrow(data__Z1_S0)) ) /
                               q_n)
      }
      
    }
    
    # Step 6: final estimates of the bounds (just doing both each time for now)
    
    l_n <- mubar_11_n * (rhobar_1_n / rhobar_0_n) + mubar_10_l_n * (1 - (rhobar_1_n / rhobar_0_n))
    u_n <- mubar_11_n * (rhobar_1_n / rhobar_0_n) + mubar_10_u_n * (1 - (rhobar_1_n / rhobar_0_n))
    
    #mean in unvaccinated infecteds for comparison
    E_Y0__S0_1 <- mean(data[[Y_name]][data[[S_name]] == 1 & data[[Z_name]] == 0])
    
    out <- c(
      E_Y0__S0_1,
      l_n,
      u_n,
      l_n - E_Y0__S0_1,
      u_n - E_Y0__S0_1,
      l_n / E_Y0__S0_1,
      u_n / E_Y0__S0_1
    )
    
    names(out) <- c("E_Y0__S0_1",
                    "E_Y1__S0_1_lower",
                    "E_Y1__S0_1_upper",
                    "additive_effect_lower",
                    "additive_effect_upper",
                    "mult_effect_lower",
                    "mult_effect_upper")
    
  } else{
    # Get rid of this condition ?? because permutation test
    # stop("Method not applicable unless evidence of vaccine protection.")
    out <- rep(NA, 7)
    
    names(out) <- c("E_Y0__S0_1",
                    "E_Y1__S0_1_lower",
                    "E_Y1__S0_1_upper",
                    "additive_effect_lower",
                    "additive_effect_upper",
                    "mult_effect_lower",
                    "mult_effect_upper")
  }
  
  return(out)
  
}


#' Covariate-adjusted bounds for naturally infected principal stratum
#'
#' Computes bounds on counterfactual outcomes within strata of a discrete
#' covariate and aggregates across strata.
#'
#' @param data A data.frame containing observed data.
#' @param X_name Name of discrete covariate used for stratification.
#' @param Y_name Outcome variable name.
#' @param Z_name Treatment variable name.
#' @param S_name Infection indicator name.
#' @param family Outcome type: "gaussian" or "binomial".
#'
#' @details
#' Levels of \code{X_name} not observed in both treatment arms are pooled to
#' ensure identifiability.
#'
#' @return A named numeric vector containing covariate-adjusted bounds for:
#' \describe{
#'   \item{E_Y0__S0_1}{Observed mean outcome among unvaccinated infected.}
#'   \item{E_Y1__S0_1_lower}{Lower bound.}
#'   \item{E_Y1__S0_1_upper}{Upper bound.}
#'   \item{additive_effect_lower}{Lower bound on additive effect.}
#'   \item{additive_effect_upper}{Upper bound on additive effect.}
#'   \item{mult_effect_lower}{Lower bound on multiplicative effect.}
#'   \item{mult_effect_upper}{Upper bound on multiplicative effect.}
#' }
#'
#' @export
get_cov_adj_bound_nat_inf <- function(
    data, 
    X_name = "X",
    Y_name = "Y",
    Z_name = "Z",
    S_name = "S",
    family = "gaussian"
){
  n <- dim(data)[1]
  
  # regroup by levels being in both vaccine and placebo arm
  x_levels_original <- sort(unique(data[[X_name]]))
  x_levels_z1 <- sort(unique(data[[X_name]][data[[Z_name]] == 1]))
  x_levels_z0 <- sort(unique(data[[X_name]][data[[Z_name]] == 0]))
  x_levels_not_in_both_z <- x_levels_original[
    ( !(x_levels_original %in% x_levels_z1) ) | ( !(x_levels_original %in% x_levels_z0 ) )
  ]
  if(length(x_levels_not_in_both_z) > 0){
    x_levels_in_both_z <- setdiff(x_levels_original, x_levels_not_in_both_z)
    if(length(x_levels_not_in_both_z) > 0){
      if(length(x_levels_in_both_z) > 0){
        x_set_level <- x_levels_in_both_z[1]
      }else{
        x_set_level <- 1
      }
      for(x_val in x_levels_not_in_both_z){
        data[[X_name]][data[[X_name]] == x_val] <- x_set_level
      }
    }
    x_levels <- sort(unique(data[[X_name]]))
  }else{
    x_levels <- x_levels_original
  }
  
  n_x_levels <- length(x_levels)
  P_Xisx_level <- rep(NA, n_x_levels)
  l_x_level <- rep(NA, n_x_levels)
  u_x_level <- rep(NA, n_x_levels)
  P_Sis1__Zis0_Xisx_level <- rep(NA, n_x_levels)
  E_Y__Zis0_Xisx_level <- rep(NA, n_x_levels)

  ct <- 0
  for(x_level in x_levels){
    ct <- ct + 1
    x_level_idx <- which(data[[X_name]] == x_level)
    n_x_level <- length(x_level_idx)
    data_x_level <- data[x_level_idx, , drop = FALSE]
    bound_nat_inf_x_level <- get_bound_nat_inf(
      data = data_x_level, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family
    )
    l_x_level[ct] <- bound_nat_inf_x_level["E_Y1__S0_1_lower"]
    u_x_level[ct] <- bound_nat_inf_x_level["E_Y1__S0_1_upper"]
    P_Xisx_level[ct] <- n_x_level / n
    P_Sis1__Zis0_Xisx_level[ct] <- mean(data_x_level[[S_name]][data_x_level[[Z_name]] == 0])
    E_Y__Zis0_Xisx_level[ct] <- mean(data_x_level[[Y_name]][data_x_level[[Z_name]] == 0 & data_x_level[[S_name]] == 1])
  }
  
  P_Sis1__Zis0 <- mean(data[[S_name]][data[[Z_name]] == 0])
  E_Y0__S0_1 <- sum(P_Sis1__Zis0_Xisx_level / P_Sis1__Zis0 * E_Y__Zis0_Xisx_level * P_Xisx_level)
  E_Y1__S0_1_lower <- sum(l_x_level * P_Sis1__Zis0_Xisx_level / P_Sis1__Zis0 * P_Xisx_level)
  E_Y1__S0_1_upper <- sum(u_x_level * P_Sis1__Zis0_Xisx_level / P_Sis1__Zis0 * P_Xisx_level)
  
  out <- c(
    "E_Y0__S0_1" = E_Y0__S0_1,
    "E_Y1__S0_1_lower" = E_Y1__S0_1_lower,
    "E_Y1__S0_1_upper" = E_Y1__S0_1_upper,
    "additive_effect_lower" = E_Y1__S0_1_lower - E_Y0__S0_1,
    "additive_effect_upper" = E_Y1__S0_1_upper - E_Y0__S0_1,
    "mult_effect_lower" = E_Y1__S0_1_lower / E_Y0__S0_1,
    "mult_effect_upper" = E_Y1__S0_1_upper / E_Y0__S0_1
  )

  return(out)
}
