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
#'   `"gcomp"`, `"ipw"`, `"aipw"`, `"tmle"`, `"bound"`, and `"sens_cw"`.
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
fit_one_glm <- function(data,
                         Y_name = "Y",
                         Z_name = "Z",
                         X_name = c("X"),
                         S_name = "S",
                         estimand = c("nat_inf", "doomed", "pop"),
                         method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens_cw"),
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
  if(any(c("sens_cw", "sens_mono") %in% method)){
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
  if(any(c("aipw", "sens_cw", "sens_mono", "ipw", "tmle") %in% method)){
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
#'   `"gcomp"`, `"ipw"`, `"aipw"`, `"tmle"`, `"bound"`, `"sens_cw"`.
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
fit_one_ml <- function(data,
                        Y_name = "Y",
                        Z_name = "Z",
                        X_name = c("X"),
                        S_name = "S",
                        estimand = c("nat_inf", "doomed", "pop"),
                        method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens_cw"),
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
  if(any(c("sens_cw", "sens_mono") %in% method)){
    out$fit_S_Z_X <- SuperLearner::SuperLearner(
      Y = data[[S_name]],
      X = data[, c(Z_name, X_name), drop = FALSE],
      family = stats::binomial(),
      SL.library = S_Z_X_library,
      cvControl = list(V = v_folds))
  }
  
  # needed for all but gcomp
  if(any(c("aipw", "sens_cw", "sens_mono", "tmle") %in% method)){
    out$fit_Z_X <- SuperLearner::SuperLearner(Y = data[[Z_name]],
                                              X = data[, X_name, drop = FALSE],
                                              family = stats::binomial(),
                                              SL.library = Z_X_library, 
                                              cvControl = list(V = v_folds))
  }
  
  return(out)
  
}

fit_models <- function(data,
                       Y_name = "G",
                       Z_name = "Z",
                       X_name = "X",
                       S_name = "S",
                       estimand = c("nat_inf", "doomed", "pop"),
                       method = c("gcomp", "ipw", "aipw", "tmle",
                                  "bound", "cov_adj_bound", "sens_cw", "sens_mono"),
                       exclusion_restriction = c(TRUE, FALSE),
                       cross_world = c(TRUE, FALSE),
                       ml = TRUE,
                       Y_Z_X_model = NULL,
                       Y_X_S1_model = NULL,
                       Y_X_S0_model = NULL,
                       S_X_model = NULL,
                       S_Z_X_model = NULL,
                       Z_X_model = paste0(Z_name, " ~ 1"),
                       Y_Z_X_library = c("SL.glm"),
                       Y_X_library = c("SL.glm"),
                       S_X_library = c("SL.glm"),
                       S_Z_X_library = c("SL.glm"),
                       Z_X_library = c("SL.mean"),
                       family = "gaussian",
                       v_folds = 3){
  
  # ml_models only for AIPW and TMLE
  model_list <- list(
    models = NULL,
    ml_models = NULL
  )
  
  # everything but bounds
  needs_models <- any(method %in% c(
    "gcomp", "ipw", "aipw", "tmle", "sens_cw", "sens_mono"
  ))
  
  if (!needs_models) {
    return(model_list)
  }
  
  # ML models only for estimators aipw, tmle, sensitivity
  if (ml && any(method %in% c("aipw", "tmle", "sens_cw", "sens_mono"))) {
    
    model_list$ml_models <- fit_one_ml(
      data = data,
      Y_name = Y_name,
      Z_name = Z_name,
      X_name = X_name,
      S_name = S_name,
      estimand = estimand,
      method = method,
      exclusion_restriction = exclusion_restriction,
      cross_world = cross_world,
      Y_Z_X_library = Y_Z_X_library,
      Y_X_library = Y_X_library,
      S_X_library = S_X_library,
      S_Z_X_library = S_Z_X_library,
      Z_X_library = Z_X_library,
      family = family,
      v_folds = v_folds
    )
  }
  
  # GLM models for gcomp/ipw, or for everything if ml = FALSE
  if (!ml || any(method %in% c("gcomp", "ipw"))) {
    
    model_list$models <- fit_one_glm(
      data = data,
      Y_name = Y_name,
      Z_name = Z_name,
      X_name = X_name,
      S_name = S_name,
      estimand = estimand,
      method = method,
      exclusion_restriction = exclusion_restriction,
      cross_world = cross_world,
      Y_Z_X_model = Y_Z_X_model,
      Y_X_S1_model = Y_X_S1_model,
      Y_X_S0_model = Y_X_S0_model,
      S_X_model = S_X_model,
      S_Z_X_model = S_Z_X_model,
      Z_X_model = Z_X_model,
      family = family
    )
  }
  
  model_list
}


#' Fit nuisance models using cross-fitting
#'
#' Fits the nuisance models required for the requested estimands and estimation
#' methods using sample splitting. The data are partitioned into cross-fitting
#' folds; for each fold, all required nuisance models are fit on the observations
#' outside that fold. The fitted models are then reorganized by nuisance model
#' type so that out-of-fold predictions can be obtained using
#' \code{\link{predict.cross_fit}}.
#'
#' @param data A \code{data.frame} containing the observed data.
#' @param id_name Optional character string giving the name of a unique subject
#'   identifier in \code{data}. If \code{NULL}, an internal identifier is created
#'   using the row order of \code{data}.
#' @param Y_name Character string giving the name of the outcome variable.
#'   Default is \code{"G"}.
#' @param Z_name Character string giving the name of the treatment variable.
#'   Default is \code{"Z"}.
#' @param X_name Character vector giving the names of baseline covariates.
#'   Default is \code{"X"}.
#' @param S_name Character string giving the name of the intermediate variable,
#'   such as infection status. Default is \code{"S"}.
#' @param estimand Character vector specifying estimands of interest. Options
#'   include \code{"nat_inf"} for the naturally infected estimand,
#'   \code{"doomed"} for the doomed estimand, and \code{"pop"} for the
#'   population-level estimand.
#' @param method Character vector specifying estimation methods. Options include
#'   \code{"gcomp"}, \code{"ipw"}, \code{"aipw"}, \code{"tmle"},
#'   \code{"bound"}, \code{"cov_adj_bound"}, \code{"sens_cw"}, and
#'   \code{"sens_mono"}.
#' @param exclusion_restriction Logical vector indicating whether the exclusion
#'   restriction is assumed.
#' @param cross_world Logical vector indicating whether cross-world assumptions
#'   are invoked.
#' @param seed Numeric seed used to randomly assign observations to
#'   cross-fitting folds.
#' @param ml Logical; whether to fit machine learning nuisance models using
#'   \code{SuperLearner}. If \code{FALSE}, only parametric generalized linear
#'   models are fit when required.
#' @param Y_Z_X_model Optional formula for the outcome regression
#'   \eqn{E[Y \mid Z, X]}. If \code{NULL}, the default formula is constructed
#'   from \code{Y_name}, \code{Z_name}, and \code{X_name}.
#' @param Y_X_S1_model Optional formula for outcome regressions among
#'   observations with \eqn{S = 1}.
#' @param Y_X_S0_model Optional formula for outcome regressions among
#'   observations with \eqn{S = 0}.
#' @param S_X_model Optional formula for models of \eqn{P(S = 1 \mid X)}.
#' @param S_Z_X_model Optional formula for models of
#'   \eqn{P(S = 1 \mid Z, X)}, used for sensitivity analyses.
#' @param Z_X_model Formula for the treatment model
#'   \eqn{P(Z = 1 \mid X)}. Default is an intercept-only model.
#' @param Y_Z_X_library Character vector of \code{SuperLearner} libraries for
#'   the outcome model \eqn{E[Y \mid Z, X]}.
#' @param Y_X_library Character vector of \code{SuperLearner} libraries for
#'   outcome models conditional on \eqn{X}.
#' @param S_X_library Character vector of \code{SuperLearner} libraries for
#'   models of \eqn{P(S = 1 \mid X)}.
#' @param S_Z_X_library Character vector of \code{SuperLearner} libraries for
#'   models of \eqn{P(S = 1 \mid Z, X)}.
#' @param Z_X_library Character vector of \code{SuperLearner} libraries for
#'   models of \eqn{P(Z = 1 \mid X)}.
#' @param family Family for outcome regression models. Passed to \code{glm} or
#'   \code{SuperLearner}. Default is \code{"gaussian"}.
#' @param v_folds Integer; number of internal cross-validation folds used by
#'   \code{SuperLearner}. Default is 3.
#' @param cf_folds Integer; number of cross-fitting folds. Must be at least 2.
#'   Default is 5.
#'
#' @details
#' Cross-fitting separates nuisance model training from nuisance prediction. For
#' each validation fold, models are fit using all observations not in that fold.
#' The returned object stores fold-specific models in a format that allows
#' \code{predict.cross_fit()} to return out-of-fold predictions for each
#' observation.
#'
#' The returned \code{fold_assignments} object contains:
#' \describe{
#'   \item{\code{id}}{The subject identifier used for fold matching.}
#'   \item{\code{idx}}{The original row index in \code{data}.}
#'   \item{\code{validation_fold}}{The validation fold assigned to each row.}
#' }
#'
#' @return An object of class \code{"cross_fit_models"}, which is a list with
#'   components:
#' \describe{
#'   \item{\code{models}}{Cross-fitted parametric nuisance models, or
#'   \code{NULL} if no parametric models were fit.}
#'   \item{\code{ml_models}}{Cross-fitted machine learning nuisance models, or
#'   \code{NULL} if no machine learning models were fit.}
#'   \item{\code{fold_assignments}}{A \code{data.frame} mapping observations to
#'   validation folds.}
#' }
#'
#' @seealso \code{\link{fit_models}}, \code{\link{reorg_cross_fit_models}},
#'   \code{\link{predict.cross_fit}}
#'
#' @export
cross_fit_models <- function(data,
                             id_name = NULL,
                             Y_name = "G",
                             Z_name = "Z",
                             X_name = "X",
                             S_name = "S",
                             estimand = c("nat_inf", "doomed", "pop"),
                             method = c("gcomp", "ipw", "aipw", "tmle",
                                        "bound", "cov_adj_bound", "sens_cw", "sens_mono"),
                             exclusion_restriction = c(TRUE, FALSE),
                             cross_world = c(TRUE, FALSE),
                             seed = 12345,
                             ml = TRUE,
                             Y_Z_X_model = NULL,
                             Y_X_S1_model = NULL,
                             Y_X_S0_model = NULL,
                             S_X_model = NULL,
                             S_Z_X_model = NULL,
                             Z_X_model = paste0(Z_name, " ~ 1"),
                             Y_Z_X_library = c("SL.glm"),
                             Y_X_library = c("SL.glm"),
                             S_X_library = c("SL.glm"),
                             S_Z_X_library = c("SL.glm"),
                             Z_X_library = c("SL.mean"),
                             family = "gaussian",
                             v_folds = 3,
                             cf_folds = 5) {
  
  set.seed(seed)
  
  # split into cf_folds equal sized validation folds and randomly shuffle
  if(cf_folds > 1){
    folds <- cut(seq(1,nrow(data)),breaks=cf_folds,labels=FALSE)
    folds <- sample(folds)
    
    fold_assignments <- data.frame(
      id = data[[id_name]],
      idx = data$.vaxstrat_idx,
      validation_fold = folds
    )
    
  } else{
    stop("cf_folds must be at least 2 for cross-fitting.")
  }
  
  # make training data and fit models
  cross_fit_model_list <- vector("list", length = cf_folds)
  for(k in 1:cf_folds){
    validation_idx <- fold_assignments$idx[
      fold_assignments$validation_fold == k
    ]
    
    training_data <- data[
      !(data$.vaxstrat_idx %in% validation_idx),
      ,
      drop = FALSE
    ]
    
    cross_fit_model_list[[k]] <- fit_models(
      data = training_data,
      Y_name = Y_name,
      Z_name = Z_name,
      X_name = X_name,
      S_name = S_name,
      estimand = estimand,
      method = method,
      exclusion_restriction = exclusion_restriction,
      cross_world = cross_world,
      ml = ml,
      Y_Z_X_model = Y_Z_X_model,
      Y_X_S1_model = Y_X_S1_model,
      Y_X_S0_model = Y_X_S0_model,
      S_X_model = S_X_model,
      S_Z_X_model = S_Z_X_model,
      Z_X_model = Z_X_model,
      Y_Z_X_library = Y_Z_X_library,
      Y_X_library = Y_X_library,
      S_X_library = S_X_library,
      S_Z_X_library = S_Z_X_library,
      Z_X_library = Z_X_library,
      family = family,
      v_folds = v_folds
    )
  }
  
  # reorganize by model type
  reorg_cross_fit <- reorg_cross_fit_models(fold_models = cross_fit_model_list,
                                            fold_assignments = fold_assignments,
                                            id_name = id_name)

  out <- list(models = reorg_cross_fit$models,
              ml_models = reorg_cross_fit$ml_models, 
              fold_assignments = fold_assignments)
  
  class(out) <- "cross_fit_models"
  
  return(out)
}

#' Reorganize fitted cross-fit models by model type
#'
#' Converts a list organized as:
#'
#' fold -> model source -> model type
#'
#' into:
#'
#' model source -> model type -> fold-specific models
#'
#' Each model-type object is assigned class `"cross_fit"` and contains the
#' corresponding fold assignments needed for out-of-fold prediction.
#'
#' @param fold_models A list of fitted model collections, with one element
#'   per cross-fitting fold.
#' @param fold_assignments A data.frame containing columns `id` and
#'   `validation_fold`.
#' @param id_name Character; name of the ID column expected in prediction data.
#'
#' @return A list containing model collections organized by model type.
#'
reorg_cross_fit_models <- function(fold_models,
                                    fold_assignments,
                                    id_name) {
  
  model_types <- c("models", "ml_models")
  
  out <- stats::setNames(
    vector("list", length(model_types)),
    model_types
  )
  
  for (mtype in model_types) {
    
    # Determine which nuisance models were fit under this type
    model_names <- unique(unlist(
      lapply(
        fold_models,
        function(x) {
          if (is.null(x[[mtype]])) {
            return(NULL)
          }
          
          names(x[[mtype]])
        }
      ),
      use.names = FALSE
    ))
    
    if (length(model_names) == 0L) {
      out[[mtype]] <- NULL
      next
    }
    
    out[[mtype]] <- stats::setNames(
      vector("list", length(model_names)),
      model_names
    )
    
    for (model_name in model_names) {
      
      fold_specific_models <- lapply(
        fold_models,
        function(x) {
          if (is.null(x[[mtype]])) {
            return(NULL)
          }
          
          x[[mtype]][[model_name]]
        }
      )
      
      names(fold_specific_models) <- 1:length(fold_specific_models)
      
      cross_fit_model <- list(fold_specific_models = fold_specific_models, 
                              fold_assignments = fold_assignments,
                              id_name = id_name)
      
      class(cross_fit_model) <- "cross_fit"
      
      out[[mtype]][[model_name]] <- cross_fit_model
    }
  }
  
  out
}

