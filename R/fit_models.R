#' Fit parametric models required for estimation procedures
#'
#' Fits all regression models needed for the specified estimands and estimation
#' methods using generalized linear models (`glm`). The function constructs and
#' fits only the models required based on the requested estimands, methods, and
#' assumptions (e.g., exclusion restriction, cross-world assumptions).
#'
#' @param data A data.frame containing the observed data.
#' @param Y_name Character; name of the outcome variable. Default is `"Y"`.
#' @param Z_name Character; name of the treatment variable. Default is `"Z"`.
#' @param X_name Character vector; names of baseline covariates. Default is `"X"`.
#' @param S_name Character; name of the intermediate (e.g., infection) variable. Default is `"S"`.
#' @param estimand Character vector specifying estimands of interest. Options include:
#'   `"nat_inf"` (naturally infected), `"doomed"`, and `"pop"` (population-level).
#' @param method Character vector specifying estimation methods. Options include
#'   `"gcomp"`, `"ipw"`, `"aipw"`, `"tmle"`, `"bound"`, and `"sens"`.
#' @param exclusion_restriction Logical; whether exclusion restriction is assumed.
#' @param cross_world Logical; whether cross-world assumptions are invoked.
#' @param Y_Z_X_model Optional formula for outcome regression \eqn{E[Y \mid Z, X]}.
#'   If NULL, defaults to `Y ~ Z + X`.
#' @param Y_X_S1_model Optional formula for outcome regression among \eqn{S = 1}.
#'   Defaults to `Y ~ X`.
#' @param Y_X_S0_model Optional formula for outcome regression among \eqn{S = 0}.
#'   Defaults to `Y ~ X`.
#' @param S_X_model Optional formula for \eqn{P(S = 1 \mid X)}. Defaults to `S ~ X`.
#' @param S_Z_X_model Optional formula for \eqn{P(S = 1 \mid Z, X)} (used for sensitivity analysis).
#' @param Z_X_model Formula for \eqn{P(Z = 1 \mid X)}. Default is intercept-only (`Z ~ 1`).
#' @param family Family for outcome regression models (passed to `glm`). Default is `"gaussian"`.
#'
#' @return A named list of fitted `glm` objects. The exact components depend on the
#' requested estimands and methods, and may include:
#' \describe{
#'   \item{fit_Y_Z_X}{Outcome model \eqn{E[Y \mid Z, X]}}
#'   \item{fit_Y_Z1_S1_X, fit_Y_Z1_S0_X, fit_Y_Z0_S1_X, fit_Y_Z0_S0_X}{Outcome models stratified by treatment and intermediate variable}
#'   \item{fit_S_Z1_X, fit_S_Z0_X}{Models for \eqn{P(S = 1 \mid Z, X)} fitted separately by treatment}
#'   \item{fit_S_Z_X}{Model for \eqn{P(S = 1 \mid Z, X)} (pooled; used for sensitivity analysis)}
#'   \item{fit_Z_X}{Model for \eqn{P(Z = 1 \mid X)}}
#' }
#'
#' @details
#' Models are fit only when required by the specified estimand–method combination.
#' For example:
#' \itemize{
#'   \item Outcome regression (`fit_Y_Z_X`) is used by g-computation, AIPW, and TMLE.
#'   \item Treatment models (`fit_Z_X`) are used by weighting-based estimators.
#'   \item Stratified outcome and intermediate models are used for principal
#'   stratification estimands (e.g., `"nat_inf"` and `"doomed"`).
#' }
#'
#' @export
fit_models <- function(data,
                       Y_name = "Y",
                       Z_name = "Z",
                       X_name = c("X"),
                       S_name = "S",
                       estimand = c("nat_inf", "doomed", "pop"),
                       method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
                       exclusion_restriction = c(TRUE, FALSE),
                       cross_world = c(TRUE, FALSE),
                       Y_Z_X_model = NULL,
                       Y_X_S1_model = NULL,
                       Y_X_S0_model = NULL,
                       S_X_model = NULL,
                       S_Z_X_model = NULL,
                       Z_X_model = paste0(Z_name, " ~ 1"),
                       family = "gaussian"){
  
  # Prep model formulas if not pre-specified
  if(is.null(Y_Z_X_model)){
    Y_Z_X_model <- stats::as.formula(paste0(Y_name, "~",
                                     Z_name, "+",
                                     paste0(X_name, collapse = "+")))
  }
  
  if(is.null(Y_X_S1_model)){
    Y_X_S1_model <- stats::as.formula(paste0(Y_name, "~",
                                      paste0(X_name, collapse = "+")))
  }
  
  if(is.null(Y_X_S0_model)){
    Y_X_S0_model <- stats::as.formula(paste0(Y_name, "~",
                                      paste0(X_name, collapse = "+")))
  }
  
  if(is.null(S_X_model)){
    S_X_model <- stats::as.formula(paste0(S_name, "~",
                                   paste0(X_name, collapse = "+")))
  }
  
  out <- list()
  
  if(
    "pop" %in% estimand | 
    ("nat_inf" %in% estimand & 
     TRUE %in% exclusion_restriction & 
     any(c("gcomp", "aipw", "tmle") %in% method))
  ){
    out$fit_Y_Z_X <- stats::glm(Y_Z_X_model, data = data, family = family)
  } 
  
  if(any(c("nat_inf", "doomed") %in% estimand)){
    sub_Z0 <- data[data[[Z_name]] == 0,]
    out$fit_S_Z0_X <- stats::glm(S_X_model, sub_Z0, family = "binomial")
    
    sub_Z1 <- data[data[[Z_name]] == 1,]
    out$fit_S_Z1_X <- stats::glm(S_X_model, sub_Z1, family = "binomial")
    
    sub_Z1_S1 <- data[data[[Z_name]] == 1 & data[[S_name]] == 1,]
    out$fit_Y_Z1_S1_X <- stats::glm(Y_X_S1_model, data = sub_Z1_S1, family = family)
    
    sub_Z1_S0 <- data[data[[Z_name]] == 1 & data[[S_name]] == 0,]
    out$fit_Y_Z1_S0_X <- stats::glm(Y_X_S0_model, data = sub_Z1_S0, family = family)
    
    sub_Z0_S1 <- data[data[[Z_name]] == 0 & data[[S_name]] == 1,]
    out$fit_Y_Z0_S1_X <- stats::glm(Y_X_S1_model, data = sub_Z0_S1, family = family)

    sub_Z0_S0 <- data[data[[Z_name]] == 0 & data[[S_name]] == 0,]
    out$fit_Y_Z0_S0_X <- stats::glm(Y_X_S0_model, data = sub_Z0_S0, family = family)
    
  }
  
  # only needed for AIPW sensitivity analysis
  if(any(c("sens") %in% method)){
    if(is.null(S_Z_X_model)){
      S_Z_X_model <- paste0(S_name, "~", Z_name, "+", paste0(X_name, collapse = "+"))
    }
    
    out$fit_S_Z_X <- stats::glm(
      S_Z_X_model,
      family = "binomial",
      data = data
    )
  }
  
  # not needed for gcomp
  if(any(c("aipw", "sens", "ipw", "tmle") %in% method)){
    out$fit_Z_X <- stats::glm(
      Z_X_model, family = "binomial", data = data
    )
  }
  
  return(out)
  
}


#' Fit machine learning models required for estimation procedures using SuperLearner
#'
#' Fits all models required for the specified estimands and estimation methods
#' using the \code{SuperLearner} ensemble learning framework. Only models needed
#' for the requested estimand–method combinations are fit.
#'
#' @param data A data.frame containing the observed data.
#' @param Y_name Character; name of the outcome variable. Default is `"Y"`.
#' @param Z_name Character; name of the treatment variable. Default is `"Z"`.
#' @param X_name Character vector; names of baseline covariates. Default is `"X"`.
#' @param S_name Character; name of the intermediate (e.g., infection) variable. Default is `"S"`.
#' @param estimand Character vector specifying estimands of interest:
#'   `"nat_inf"`, `"doomed"`, `"pop"`.
#' @param method Character vector specifying estimation methods:
#'   `"gcomp"`, `"ipw"`, `"aipw"`, `"tmle"`, `"bound"`, `"sens"`.
#' @param exclusion_restriction Logical; whether exclusion restriction is assumed.
#' @param cross_world Logical; whether cross-world assumptions are invoked.
#' @param Y_Z_X_library Character vector of SuperLearner libraries for outcome model
#'   \eqn{E[Y \mid Z, X]}. Default is `"SL.glm"`.
#' @param Y_X_library Character vector of SuperLearner libraries for outcome models
#'   conditional on \eqn{X}. Default is `"SL.glm"`.
#' @param S_X_library Character vector of libraries for \eqn{P(S = 1 \mid X)}.
#' @param S_Z_X_library Character vector of libraries for \eqn{P(S = 1 \mid Z, X)}.
#' @param Z_X_library Character vector of libraries for \eqn{P(Z = 1 \mid X)}.
#'   Default is `"SL.mean"`.
#' @param family Family for outcome regression (passed to SuperLearner). Default is `"gaussian"`.
#' @param v_folds Integer; number of cross-validation folds. Default is 3.
#'
#' @return A named list of fitted \code{SuperLearner} objects. Components mirror
#' those returned by \code{fit_models()}, including:
#' \describe{
#'   \item{fit_Y_Z_X}{Outcome model \eqn{E[Y \mid Z, X]}}
#'   \item{fit_Y_Z1_S1_X, fit_Y_Z1_S0_X, fit_Y_Z0_S1_X, fit_Y_Z0_S0_X}{Stratified outcome models}
#'   \item{fit_S_Z1_X, fit_S_Z0_X}{Models for \eqn{P(S = 1 \mid Z, X)}}
#'   \item{fit_S_Z_X}{Model for \eqn{P(S = 1 \mid Z, X)} (pooled)}
#'   \item{fit_Z_X}{Model for \eqn{P(Z = 1 \mid X)}}
#' }
#'
#' @details
#' This function provides a flexible alternative to \code{fit_models()} by allowing
#' nonparametric estimation via ensemble learning. Users can specify custom
#' SuperLearner libraries for each component of the estimation procedure.
#'
#' @seealso \code{\link[SuperLearner]{SuperLearner}}
#'
#' @export
fit_ml_models <- function(data,
                          Y_name = "Y",
                          Z_name = "Z",
                          X_name = c("X"),
                          S_name = "S",
                          estimand = c("nat_inf", "doomed", "pop"),
                          method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
                          exclusion_restriction = c(TRUE, FALSE),
                          cross_world = c(TRUE, FALSE),
                          Y_Z_X_library = c("SL.glm"),
                          Y_X_library = c("SL.glm"),
                          S_X_library = c("SL.glm"),
                          S_Z_X_library = c("SL.glm"),
                          Z_X_library = c("SL.mean"),
                          family = "gaussian",
                          v_folds = 3){
  
  out <- list()
  
  if(
    "pop" %in% estimand | 
    ("nat_inf" %in% estimand & 
     TRUE %in% exclusion_restriction & 
     any(c("gcomp", "aipw", "tmle") %in% method))
  ){
    out$fit_Y_Z_X <- SuperLearner::SuperLearner(Y = data[[Y_name]],
                                                X = data[, colnames(data) %in% c(Z_name, X_name), drop = FALSE],
                                                family = family,
                                                SL.library = Y_Z_X_library, 
                                                cvControl = list(V = v_folds))
    
  } 
  
  # needed for any weight-based estimator
  if(any(c("nat_inf", "doomed") %in% estimand)){
    sub_Z0 <- data[data[[Z_name]] == 0,]
    out$fit_S_Z0_X <- SuperLearner::SuperLearner(Y = sub_Z0[[S_name]],
                                                 X = sub_Z0[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = S_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_Z1 <- data[data[[Z_name]] == 1,]
    out$fit_S_Z1_X <- SuperLearner::SuperLearner(Y = sub_Z1[[S_name]],
                                                 X = sub_Z1[, X_name, drop = FALSE],
                                                 family = stats::binomial(),
                                                 SL.library = S_X_library, 
                                                 cvControl = list(V = v_folds))
    
    sub_Z1_S1 <- data[data[[Z_name]] == 1 & data[[S_name]] == 1,]
    out$fit_Y_Z1_S1_X <- SuperLearner::SuperLearner(Y = sub_Z1_S1[[Y_name]],
                                                    X = sub_Z1_S1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(V = v_folds))
    
    sub_Z1_S0 <- data[data[[Z_name]] == 1 & data[[S_name]] == 0,]
    out$fit_Y_Z1_S0_X <- SuperLearner::SuperLearner(Y = sub_Z1_S0[[Y_name]],
                                                    X = sub_Z1_S0[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(V = v_folds))
    
    sub_Z0_S1 <- data[data[[Z_name]] == 0 & data[[S_name]] == 1,]
    out$fit_Y_Z0_S1_X <- SuperLearner::SuperLearner(Y = sub_Z0_S1[[Y_name]],
                                                    X = sub_Z0_S1[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(V = v_folds))

    sub_Z0_S0 <- data[data[[Z_name]] == 0 & data[[S_name]] == 0,]
    out$fit_Y_Z0_S0_X <- SuperLearner::SuperLearner(Y = sub_Z0_S0[[Y_name]],
                                                    X = sub_Z0_S0[, X_name, drop = FALSE],
                                                    family = family,
                                                    SL.library = Y_X_library, 
                                                    cvControl = list(V = v_folds))

  }
  
  # only needed for AIPW sensitivity analysis
  if(any(c("sens") %in% method)){
    out$fit_S_Z_X <- SuperLearner::SuperLearner(
      Y = data[[S_name]],
      X = data[, c(Z_name, X_name), drop = FALSE],
      family = "binomial",
      SL.library = S_Z_X_library,
      cvControl = list(V = v_folds))
  }
  
  # needed for all but gcomp
  if(any(c("aipw", "sens", "tmle") %in% method)){
    out$fit_Z_X <- SuperLearner::SuperLearner(Y = data[[Z_name]],
                                              X = data[, X_name, drop = FALSE],
                                              family = stats::binomial(),
                                              SL.library = Z_X_library, 
                                              cvControl = list(V = v_folds))
  }
  
  return(out)
  
}
