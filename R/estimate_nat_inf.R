#' Function for covariate-unadjusted estimation of counterfactual post-infection outcomes in the 
#' naturally infected principal strata
#' @inheritParams vegrowth
#' 
#' @export
#' 
#' @return unadjusted estimates of growth effect in the naturally infected strata
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

#' Function for g-computation of counterfactual post-infection outcomes in the 
#' naturally infected principal strata
#' 
#' @inheritParams vegrowth
#' @param models A list of fitted models returned from \code{fit_models}
#'  
#' @export
#' 
#' @return g-comp estimate of growth effect in the naturally infected strata
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
    mu_01_X <- simple_predict(models$fit_Y_Z0_S1_X, newdata = data)
    pi_1_X <- simple_predict(models$fit_Z_X, newdata = data)
    pi_0_X <- 1 - pi_1_X
    rho_bar_0 <- mean(rho_0_X)
    
    # psi_0 = Weight * E[E[Y | Z = 0, Y = 1, X]]
    psi_tilde_0_X <- rho_0_X / rho_bar_0 * mu_01_X
    
    psi_0 <- mean( psi_tilde_0_X )
    
    mu_11_X <- simple_predict(models$fit_Y_Z1_S1_X, newdata = data)
    mu_10_X <- simple_predict(models$fit_Y_Z1_S0_X, newdata = data)
    mu_00_X <- simple_predict(models$fit_Y_Z0_S0_X, newdata = data)
    mu_dot0_X <- pi_1_X * mu_10_X + pi_0_X * mu_00_X

    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
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

#' Function for IPW of counterfactual post-infection outcomes in the 
#' naturally infected principal strata
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param S_name TODO
#' @param Y_name TODO
#' @param Z_name TODO
#' 
#' @return IPW estimate of growth effect in the naturally infected principal stratum
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

#' Function for efficient AIPW estimator
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param Y_name name of growth outcome variable, default Y
#' @param Z_name name of vaccine treatment variable, default Z
#' @param S_name name of infection variable, default Y
#' @param return_se flag to return standard error, defualt FALSE
#' @param exclusion_restriction boolean indicating whether an exclusion restriction is assumed
#' @param cross_world boolean indicating whether a cross-world independence assumption is assumed
#' 
#' @export
#' 
#' @return AIPW estimate of growth effect in naturally infected strata (+ standard error if return_se = TRUE)
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
    mu_dot0_X <- pi_1_X * mu_10_X + pi_0_X * mu_00_X

    rho_1_X <- simple_predict(models$fit_S_Z1_X, newdata = data)
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
  se <- sqrt(var(augmentation_1 - augmentation_0) / dim(data)[1])
  
  se_psi_1 <- sqrt(var(augmentation_1) / dim(data)[1])
  se_psi_0 <- sqrt(var(augmentation_0) / dim(data)[1])
  
  # Multiplicative effect (log scale)
  efficient_growth_effect_log_mult <- log(psi_1_aipw / psi_0_aipw)
  
  # Get SE using IF matrix same way as TMLE
  if_matrix <- cbind(augmentation_1, augmentation_0)
  colnames(if_matrix) <- c("augmentation_1", "augmentation_2")
  
  cov_matrix <- cov(if_matrix) / dim(data)[1]

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

#' Function for efficient TMLE estimator
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param Y_name name of growth outcome variable, default Y
#' @param Z_name name of vaccine treatment variable, default Z
#' @param S_name name of infection variable, default Y
#' @param return_se flag to return standard error, defualt FALSE
#' @param max_iter TODO 
#' @param tol TOOD 
#' 
#' @return TMLE estimate of growth effect (+ standard error if return_se = TRUE)
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
      target_fit <- suppressWarnings(glm(
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
      target_fit <- suppressWarnings(glm(
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
      target_fit <- suppressWarnings(glm(
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
      if(sd(H1) > 0 & sd(H0) > 0){
        if(cor(H1, H0) > 0.99999){
          H0 <- 1
        }
      
        target_data <- data.frame(
          S_inf = data[[S_name]],
          target_wt = target_wt,
          H1 = H1,
          H0 = H0,
          logit_rho_0_star = logit_rho_0_star
        )
        target_data <- setNames(target_data, c(S_name, names(target_data[-1])))
        
        # include intercept so rho_bar_0_star is still mean(Y[Z == 0])
        target_fit <- glm(
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
      target_data <- setNames(target_data, c(S_name, names(target_data[-1])))
      
      # include intercept so rho_bar_0_star is still mean(Y[Z == 0])
      target_fit <- glm(
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
    se <- sqrt(var(phi_ge_data) / dim(data)[1])
    
    if_matrix <- cbind(phi_1_data, phi_0_data)
    cov_matrix <- cov(if_matrix) / dim(data)[1]
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

#' Function for efficient AIPW estimator for sensitivity analysis
#' 
#' @param data dataset to predict on
#' @param models list of pre-fit models needed for estimation
#' @param Y_name name of growth outcome variable, default Y
#' @param Z_name name of vaccine treatment variable, default Z
#' @param S_name name of infection variable, default Y
#' @param epislon a vector of values for the sensitivity parameter
#' @param return_se flag to return standard error, defualt FALSE
#' 
#' 
#' @return AIPW estimate of growth effect (+ standard error if return_se = TRUE)
#' 
do_sens_aipw_nat_inf <- function(data,
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
  
  data_0 <- data; data[[Z_name]] <- 0
  data_1 <- data; data[[Z_name]] <- 1
  
  
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
    mu_11_X <- predict(models$fit_Y_Z1_S1_X, newdata = data, type = "response")
    mu_10_X <- predict(models$fit_Y_Z1_S0_X, newdata = data, type = "response")
    rho_1_X <- predict(models$fit_S_Z_X, newdata = data_1, type = "response")
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
    sqrt(var(augmentation_1_eps - augmentation_0) / dim(data)[1])
  })
  
  # Multiplicative effect (log scale)
  efficient_growth_effect_log_mult_epsilon <- lapply(psi_1_epsilon_aipw, function(psi_1_eps){
    log(psi_1_eps / psi_0_aipw)
  })
  
  cov_matrices <- mapply(
    augmentation_1_eps = augmentation_1_epsilon, 
    function(augmentation_1_eps){
      if_matrix <- cbind(augmentation_1_eps, augmentation_0)
      cov_matrix <- cov(if_matrix) / dim(data)[1]
      return(cov_matrix)
    }, SIMPLIFY =  FALSE
  )
  
  # Get SE using IF matrix same way as TMLE
  se_log_mult_eff <- mapply(
    augmentation_1_eps = augmentation_1_epsilon, 
    psi_1_eps_aipw = psi_1_epsilon_aipw,
    function(augmentation_1_eps, psi_1_eps_aipw){
      if_matrix <- cbind(augmentation_1_eps, augmentation_0)
      cov_matrix <- cov(if_matrix) / dim(data)[1]
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
  
  class(out) <- "sens"
  
  return(out)
}

#' Function for bounds on naturally infected estimate without use of cross-world assumption
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param S_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary outcome
#' 
#' @export
#' 
#' @return list containing estimate of E[Y(0) | Y(0) = 1], bounds on E[Y(1) | Y(0) = 1], bounds on additive effect, bounds on multiplicative effect
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
  
  # get rid of this condition bc permutation test?
  if(rhobar_0_n > rhobar_1_n){
    # Step 2: mubar_11_n 
    mubar_11_n_num <- sum(data[[Y_name]]*data[[S_name]]*data[[Z_name]])
    mubar_11_n_denom <- sum(data[[S_name]]*data[[Z_name]])
    
    if(mubar_11_n_denom == 0){ # denominator NA in some bootstrap replicates, change to 0
      mubar_11_n <- 0
    } else{
      mubar_11_n <- mubar_11_n_num / mubar_11_n_denom
    }
    
    # Step 3: q_n (relative size of protected? in (immune + protected) in vax)
    q_n = 1 - (1 - rhobar_0_n) / (1 - rhobar_1_n)
    
    # Step 4: q_n^th quintiles of S__Z1_S0 (aka Y__Z1_S0, need to rename everything at some point)
    Y__Z1_S0 <- data[[Y_name]][which(data[[Z_name]] == 1 & data[[S_name]] == 0)]
    q_nth_quintile <- quantile(Y__Z1_S0, probs = q_n) # NOTE failing here if condition on line 709 not met
    one_minus_q_nth_quintile <- quantile(Y__Z1_S0, probs = 1 - q_n)
    
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
    
    # Step 6: final estimates of the bounds (just doing both each time for now??)
    
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


#' Function for covariate-adjusted bounds on naturally infected estimate
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param X_name discrete-valued covariate name
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param S_name infection variable name
#' @param family gaussian for continuous outcome, binomial for binary outcome
#' 
#' @export
#' 
#' @return list containing estimate of E[Y(0) | Y(0) = 1], bounds on E[Y(1) | Y(0) = 1], bounds on additive effect, bounds on multiplicative effect
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
  E_Y1__S0_1_lower <- sum(l_x_level * P_Xisx_level)
  E_Y1__S0_1_upper <- sum(u_x_level * P_Xisx_level)
  
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

# ------------------------------------------------------------------------------
# Old or in-progress
# ------------------------------------------------------------------------------

#' Function for Hudgens-style bounds on effect that incorporate covariates
#' 
#' Currently assumes that the conditional mean of Y follows a linear model
#' with Normal errors.
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param models list of pre-fit models needed for estimation
#' @param family gaussian for continuous outcome Y, binomial for binary outcome
#' @param lower_bound A boolean. If TRUE, then adds the smallest growth measures 
#'    to the infected vaccines thereby yielding a lower
#'    bound on the effect of interest. If FALSE, then adds the largest
#'    growth measures to the infected vacccinees thereby yielding an upper
#'    bound on the effect of interest.
#' 
#' @examples
#' 
#' n <- 1e4
#' X <- sample(seq(-1,1), n, replace = TRUE)
#' p_immune <- 0.5 + 0.25 * X
#' p_doomed <- 0.1 + 0.05 * X
#' p_helped <- 1 - (p_immune + p_doomed)
#' ps <- mapply(
#'   p_i = p_immune, p_d = p_doomed, p_h = p_helped, 
#'   FUN = function(p_i, p_d, p_h){
#'     sample(
#'       c("immune", "doomed", "helped"), 
#'       size = 1, prob = c(p_i, p_d, p_h)
#'     )
#'   }
#' )
#' 
#' Z <- rbinom(n, 1, 0.5)
#' Y0 <- ifelse(ps == "immune", 0, 1)
#' Y1 <- ifelse(ps == "doomed", 1, 0)
#' Y <- ifelse(Z == 1, Y1, Y0)
#' Y1 <- 1*X - 0.5 * Y1 + rnorm(n, 0, 0.5)
#' Y0 <- 1*X - 0.5 * Y0 + rnorm(n, 0, 0.5)
#' Y <- ifelse(Z == 1, Y1, Y0)
#' 
#' marginal_effect <- mean(Y1 - Y0)
#' ps_effect <- mean(Y1[Y0 == 1] - Y0[Y0 == 1])
#' 
#' data <- data.frame(X, Z, Y, Y)
#' models <- fit_models(data)
#' 
#' get_adjusted_hudgens_stat(data, models, lower_bound = TRUE)
#' # compare to unadjusted
#' get_hudgens_stat(data, lower_bound = TRUE)
#' 
#' get_adjusted_hudgens_stat(data, models, lower_bound = FALSE)
#' # compare to unadjusted
#' get_hudgens_stat(data, lower_bound = FALSE)
#' 
#' # binary outcome
#' Y_binary <- as.numeric(Y > 1)
#' data <- data.frame(X, Z, Y, Y = Y_binary)
#' models <- fit_models(data, family = binomial())
#' get_adjusted_hudgens_stat(data, models, family = "binomial", lower_bound = TRUE)
#' get_adjusted_hudgens_stat(data, models, family = "binomial", lower_bound = FALSE)
#' 
#' 
#' @return Hudgens-style estimate of bound on effect in naturally infected

# get_adjusted_hudgens_stat <- function(
#     data, 
#     models,
#     family = "gaussian",
#     lower_bound = TRUE
# ){
#   
#   E_Y_Z0_S1_X <- predict(models$fit_Y_Z0_S1_X, newdata = data, type = "response")
#   
#   E_Y_Z1_S1_X <- predict(models$fit_Y_Z1_S1_X, newdata = data, type = "response")
#   E_Y_Z1_S0_X <- predict(models$fit_Y_Z1_S0_X, newdata = data, type = "response")
#   
#   P_S1_Z1_X <- predict(models$fit_S_Z1_X, newdata = data, type = "response")
#   P_S1_Z0_X <- predict(models$fit_S_Z0_X, newdata = data, type = "response")
#   P_S0_Z1_X <- 1 - P_S1_Z1_X
#   P_S0_Z0_X <- 1 - P_S1_Z0_X
#   
#   P_S1_Z0 <- mean(P_S1_Z0_X)
#   
#   VE_X <- 1 - ( P_S1_Z1_X / P_S1_Z0_X )
#   if(any(VE_X < 0)){
#     warning("Some condtional ZE estimates < 0 -- truncating these estimates at 0.")
#   }
#   VE_X[VE_X < 0] <- 0
#   ZE_is_zero <- (VE_X == 0)
#   ZE_is_nonzero <- (VE_X > 0)
#   
#   q_X_low <- 1 - P_S0_Z0_X / P_S0_Z1_X
#   q_X_high <- 1 - q_X_low
#   
#   if(family == "gaussian"){
#     sd_Y <- (mean(models$fit_Y_Z1_S0_X$residuals^2))^(1/2)
#   }
#   
#   E_Y_Z1_S0_truncY_X <- E_Y_Z1_S0_X
#   
#   if(lower_bound){
#     
#     if(family == "gaussian"){
#       # calculate mean of Normal given less than q_X_low
#       beta_X <- (q_X_low - E_Y_Z1_S0_X) / sd_Y
#       E_Y_Z1_S0_truncY_X[ZE_is_nonzero] <- E_Y_Z1_S0_X[ZE_is_nonzero] - sd_Y * dnorm(beta_X[ZE_is_nonzero]) / pnorm(beta_X[ZE_is_nonzero])
#     }else{
#       E_Y_Z1_S0_truncY_X[ZE_is_nonzero] <- as.numeric(E_Y_Z1_S0_X[ZE_is_nonzero] > q_X_low)
#     }
#   }else{
#     
#     if(family == "gaussian"){
#       # calculate mean of Normal given greater than q_X_high
#       alpha_X <- (q_X_high - E_Y_Z1_S0_X) / sd_Y
#       E_Y_Z1_S0_truncY_X[ZE_is_nonzero] <- E_Y_Z1_S0_X[ZE_is_nonzero] + sd_Y * dnorm(alpha_X[ZE_is_nonzero]) / pnorm(alpha_X[ZE_is_nonzero], lower.tail = FALSE)
#     }else{
#       E_Y_Z1_S0_truncY_X[ZE_is_nonzero] <- as.numeric(E_Y_Z1_S0_X[ZE_is_nonzero] > q_X_high)
#     }
#   }
#   
#   E_Y_Z1_S0_X_bound <- E_Y_Z1_S1_X * (1 - VE_X) + E_Y_Z1_S0_truncY_X * VE_X
#   
#   effect <- mean(
#     P_S1_Z0_X / P_S1_Z0 * (E_Y_Z1_S0_X_bound - E_Y_Z0_S1_X)
#   )
#   
#   return(effect)
#   
# }

#' Function to get chop-lump style test-statistic
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param S_name infection variable name
#' 
#' @return dataframe with chop-lump style test statistics for mean in vax, mean in placebo
#' 
# get_chop_lump_statistic <- function(data,
#                                     Y_name = "Y",
#                                     Z_name = "Z",
#                                     S_name = "S"){
#   
#   # 1. Yet number of people in relevant groups
#   n_no_inf_plc <- sum(data[[S_name]] == 0 & data[[Z_name]] == 0)
#   n_no_inf_vax <- sum(data[[S_name]] == 0 & data[[Z_name]] == 1)
#   n_plc <- sum(data[[Z_name]] == 0)
#   n_vax <- sum(data[[Z_name]] == 1)
#   n_inf_plc <- n_plc - n_no_inf_plc
#   n_inf_vax <- n_vax - n_no_inf_vax
#   
#   # 2. Chop based on group with more infections
#   if(n_inf_plc > n_inf_vax){
#     
#     # placebo - everyone is infected, simple mean
#     mean_Y_plc <- mean(data[[Y_name]][data[[Z_name]] == 0 & data[[S_name]] == 1])
#     
#     # vaccinated - weighted mean
#     mean_Y_noinf_vax <- mean(data[[Y_name]][data[[Z_name]] == 1 & data[[S_name]] == 0])
#     mean_Y_inf_vax <- mean(data[[Y_name]][data[[Z_name]] == 1 & data[[S_name]] == 1])
#     
#     mean_Y_vax  <- mean_Y_noinf_vax * ((n_inf_plc - n_inf_vax) / n_inf_plc) +
#       mean_Y_inf_vax * (n_inf_vax / n_inf_plc)
#     
#   } else if (n_inf_plc < n_inf_vax){
#     
#     # vaccinated - everyone is infected, simple mean
#     mean_Y_vax <- mean(data[[Y_name]][data[[Z_name]] == 1 & data[[S_name]] == 1])
#     
#     # placebo - weighted mean
#     mean_Y_noinf_plc <- mean(data[[Y_name]][data[[Z_name]] == 0 & data[[S_name]] == 0])
#     mean_Y_inf_plc <- mean(data[[Y_name]][data[[Z_name]] == 0 & data[[S_name]] == 1])
#     
#     mean_Y_plc  <- mean_Y_noinf_plc * ((n_inf_vax - n_inf_plc) / n_inf_vax) +
#       mean_Y_inf_plc * (n_inf_plc / n_inf_vax)
#     
#   } else{
#     
#     # vaccinated - everyone is infected, simple mean
#     mean_Y_vax <- mean(data[[Y_name]][data[[Z_name]] == 1 & data[[S_name]] == 1])
#     
#     # placebo - everyone is infected, simple mean
#     mean_Y_plc <- mean(data[[Y_name]][data[[Z_name]] == 0 & data[[S_name]] == 1])
#     
#   }
#   
#   return(data.frame(mean_Y_plc = mean_Y_plc,
#                     mean_Y_vax = mean_Y_vax))
# }

#' Function to do permutation test for chop-lump style test-statistics
#' 
#' @param data dataframe containing dataset to use for analysis
#' @param Y_name growth outcome variable name
#' @param Z_name vaccination variable name
#' @param S_name infection variable name
#' @param n_permutations number of permutations to complete
#' 
#' @return chop lump test statistic
# do_chop_lump_test <- function(data, 
#                               Y_name = "Y",
#                               Z_name = "Z",
#                               S_name = "S",
#                               n_permutations = 1e4){
#   
#   original_means <- get_chop_lump_statistic(data, 
#                                             Y_name = Y_name,
#                                             Z_name = Z_name,
#                                             S_name = S_name)
#   ## Permutation approach
#   null_means <- vector("list", length = n_permutations)
#   for(i in 1:n_permutations){
#     data_shuffle <- data
#     data_shuffle[[Z_name]] <- sample(data_shuffle[[Z_name]])
#     
#     null_means[[i]] <- get_chop_lump_statistic(data_shuffle,
#                                                Y_name = Y_name,
#                                                Z_name = Z_name,
#                                                S_name = S_name)
#   }
#   
#   null_df <- do.call(rbind, null_means)
#   
#   ## Hypothesis test
#   observed_diff <- original_means$mean_Y_vax - original_means$mean_Y_plc
#   null_df$mean_diff <- null_df$mean_Y_vax - null_df$mean_Y_plc
#   
#   out <- list(
#     obs_diff = observed_diff,
#     null_diffs = null_df$mean_diff,
#     pval = mean(null_df$mean_diff > observed_diff)
#   )
#   
#   class(out) <- "chop_lump_res"
#   return(out)
# }


#' Helper function to plot chop-lump test results
#' 
#' @param out chop_lump_res object to plot
#' 
#' @return A histogram plot with the distribution of test statistics under the null hypothesis,
#' and a vertical line showing the observed test statistic.

# plot.chop_lump_res <- function(out){
#   hist(
#     out$null_diffs,
#     xlab = "Test statistic",
#     xlim = range(c(out$null_diffs, out$observed_diff)),
#     main = "Distribution of test statistic under null"
#   )
#   abline(v = out$observed_diff, col = 2, lwd = 2)
# }


