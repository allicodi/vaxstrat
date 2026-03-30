#' Principal stratification analysis for vaccine effects on post-infection outcomes
#'
#' \code{vaxstrat()} is the main entry point for the \pkg{vaxstrat} package. It implements
#' principal stratification methods to estimate causal effects of vaccination on outcomes
#' defined after infection (e.g., growth following infection). The function supports
#' multiple estimands (naturally infected, doomed, population-level) 
#' and estimation strategies (g-computation, inverse
#' probability weighting (IPW), augmented IPW (AIPW), targeted maximum likelihood
#' estimation (TMLE), bounds, and sensitivity analyses)
#'
#' @param data A data frame containing the analysis dataset.
#' @param Y_name Character. Name of the outcome variable (e.g., growth).
#' @param Z_name Character. Name of the treatment (vaccination) variable.
#' @param X_name Character vector. Name(s) of baseline covariates.
#' @param S_name Character. Name of the post-treatment intermediate variable (e.g., infection indicator).
#'
#' @param estimand Character vector specifying estimand(s) of interest:
#' \describe{
#'   \item{"nat_inf"}{Effect among the naturally infected principal strata.}
#'   \item{"doomed"}{Effect among the doomed principal stratum.}
#'   \item{"pop"}{Population-level (marginal) effect.}
#' }
#'
#' @param method Character vector specifying estimation method(s):
#' \describe{
#'   \item{"gcomp"}{G-computation.}
#'   \item{"ipw"}{Inverse probability weighting.}
#'   \item{"aipw"}{Augmented inverse probability weighting.}
#'   \item{"tmle"}{Targeted maximum likelihood estimation (naturally infected only).}
#'   \item{"bound"}{Nonparametric bounds without cross-world assumptions.}
#'   \item{"cov_adj_bound"}{Covariate-adjusted bounds.}
#'   \item{"sens"}{Sensitivity analysis for violations of cross-world assumptions.}
#' }
#'
#' @param exclusion_restriction Logical or vector of logicals. Whether to impose the
#' exclusion restriction assumption for naturally infected estimators.
#'
#' @param cross_world Logical or vector of logicals. Whether to impose cross-world
#' assumptions required for point identification of naturally infected effects.
#'
#' @param two_part_model Logical. If \code{TRUE}, models \eqn{E(Y \mid Z, X)} via
#' separate models for \eqn{E(Y \mid Z, X, S)} and \eqn{P(S \mid Z, X)}; otherwise uses
#' a single model. Currently implemented for select estimators.
#'
#' @param n_boot Integer. Number of bootstrap replicates for standard error estimation.
#' @param permutation Logical. Whether to perform permutation-based inference for bounds.
#' @param n_perm Integer. Number of permutations if \code{permutation = TRUE}.
#' @param seed Integer. Random seed for reproducibility.
#'
#' @param return_se Logical. If \code{TRUE}, uses closed-form standard errors when available
#' (AIPW/TMLE); otherwise uses bootstrap.
#'
#' @param ml Logical. If \code{TRUE}, uses SuperLearner for nuisance estimation in AIPW/TMLE.
#'
#' @param Y_Z_X_model Optional model formula for \eqn{E(Y \mid Z, X)}.
#' @param Y_X_S1_model Optional model formula for \eqn{E(Y \mid X, S = 1)}.
#' @param Y_X_S0_model Optional model formula for \eqn{E(Y \mid X, S = 0)}.
#' @param S_X_model Optional model formula for \eqn{P(S \mid X)}.
#' @param S_Z_X_model Optional model formula for \eqn{P(S \mid Z, X)}.
#' @param Z_X_model Model formula for \eqn{P(Z \mid X)}. Defaults to intercept-only.
#'
#' @param Y_Z_X_library Character vector of SuperLearner libraries for \eqn{E(Y \mid Z, X)}.
#' @param Y_X_library Character vector of SuperLearner libraries for \eqn{E(Y \mid X)}.
#' @param S_X_library Character vector of SuperLearner libraries for \eqn{P(S \mid X)}.
#' @param S_Z_X_library Character vector of SuperLearner libraries for \eqn{P(S \mid Z, X)}.
#' @param Z_X_library Character vector of SuperLearner libraries for \eqn{P(Z \mid X)}.
#'
#' @param null_hypothesis_value Numeric. Null value for hypothesis testing (default 0).
#' @param alpha_level Numeric. Significance level for hypothesis tests (default 0.05).
#' @param return_models Logical. If \code{TRUE}, returns fitted nuisance models.
#'
#' @param family Character. Outcome family (e.g., \code{"gaussian"}).
#' @param v_folds Integer. Number of cross-validation folds for SuperLearner.
#'
#' @param effect_dir Character. Direction of beneficial effect (\code{"positive"} or \code{"negative"}).
#' Used for one-sided inference with bounds.
#'
#' @param epsilon Numeric vector. Sensitivity parameter values (used when \code{method = "sens"}).
#'
#' @return An object of class \code{"vaxstrat"},
#' structured as a nested list with components corresponding to each estimand:
#'
#' \describe{
#'   \item{\code{nat_inf}}{Results for the naturally infected principal strata.}
#'   \item{\code{doomed}}{Results for the doomed principal stratum.}
#'   \item{\code{pop}}{Results for the population (marginal) estimand.}
#' }
#'
#' Each estimand contains one or more estimator-specific results (e.g.,
#' \code{aipw}, \code{gcomp}, \code{bound}, \code{sens}), with structure depending
#' on the method:
#'
#' \describe{
#'   \item{\code{pt_est}}{Point estimates. Typically includes:
#'     \itemize{
#'       \item \code{additive_effect}
#'       \item \code{log_multiplicative_effect}
#'       \item corresponding standard errors (if available)
#'     }
#'   }
#'
#'   \item{\code{boot_se}}{Bootstrap-based standard errors and confidence intervals
#'   (for methods using resampling). May include separate lower/upper bound SEs
#'   for partially identified estimands.}
#'
#'   \item{\code{test_stat}}{Test statistics for hypothesis tests (additive and multiplicative scales).}
#'
#'   \item{\code{p_val}}{P-values corresponding to test statistics.}
#'
#'   \item{\code{reject}}{Logical indicators for rejection of the null hypothesis
#'   at the specified \code{alpha_level}.}
#'
#'   \item{\code{permutation}}{(Bounds only, optional) Results from permutation-based
#'   inference, including null distributions and permutation p-values.}
#' }
#'
#' For sensitivity analyses (\code{method = "sens"}), results additionally include:
#' \describe{
#'   \item{\code{epsilon}}{Grid of sensitivity parameter values.}
#'   \item{\code{pt_est}}{Estimates as a function of \code{epsilon}.}
#'   \item{\code{reject}}{Data frames summarizing hypothesis tests across \code{epsilon}.}
#' }
#'
#' If \code{return_models = TRUE}, the output also includes:
#' \describe{
#'   \item{\code{models}}{List of fitted nuisance models (GLM and/or SuperLearner).}
#' }
#'
#' @details
#' The function supports estimation under varying identification assumptions:
#' \itemize{
#'   \item Exclusion restriction
#'   \item Cross-world assumptions
#' }
#'
#' Bootstrap is used for inference unless closed-form standard errors are available
#' (AIPW and TMLE).
#'
#' @examples
#' \dontrun{
#' fit <- vaxstrat(
#'   data = trial_data,
#'   Y_name = "growth",
#'   Z_name = "vaccine",
#'   S_name = "infection",
#'   X_name = c("age", "sex"),
#'   estimand = "nat_inf",
#'   method = c("aipw", "tmle")
#' )
#'
#' summary(fit)
#' }
#'
#' @export
vaxstrat <- function(data,
                     Y_name = "G",
                     Z_name = "Z",
                     X_name = "X",
                     S_name = "S",
                     estimand = c("nat_inf", "doomed", "pop"),
                     method = c("gcomp", "ipw", "aipw", "tmle", "bound", "cov_adj_bound", "sens"),
                     exclusion_restriction = c(TRUE, FALSE),
                     cross_world = c(TRUE, FALSE),
                     two_part_model = FALSE,
                     n_boot = 1000,
                     permutation = FALSE,
                     n_perm = 1000,
                     seed = 12345,
                     return_se = TRUE,
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
                     null_hypothesis_value = 0,
                     alpha_level = 0.05,
                     return_models = TRUE,
                     family = "gaussian",
                     v_folds = 3,
                     effect_dir = "positive",
                     epsilon = exp(seq(log(0.5), log(2), length = 49))){
  
  set.seed(seed)
  
  # ----------------------------------------------------------------------------
  # 1. Model Fitting -----------------------------------------------------------
  # ----------------------------------------------------------------------------
  
  # ml_models only for AIPW and TMLE
  model_list <- list(models = NULL, 
                     ml_models = NULL)
  
  # Estimation methods requiring model fitting (everything but bounds)
  if(any(method %in% c("gcomp", "ipw", "aipw", "tmle", "sens"))){
    
    # If ML specified and aipw, tmle, and/or sens are in method, fit ML models; otherwise only fit GLMs
    if(ml){
      if(any(method %in% c("aipw", "tmle", "sens"))){
        ml_models <- vaxstrat::fit_ml_models(data = data, 
                                             estimand = estimand,
                                             method = method,
                                             exclusion_restriction = exclusion_restriction,
                                             cross_world = cross_world, 
                                             Y_name = Y_name,
                                             Z_name = Z_name,
                                             S_name = S_name,
                                             X_name = X_name,
                                             Y_Z_X_library = Y_Z_X_library,
                                             Y_X_library = Y_X_library,
                                             S_X_library = S_X_library,
                                             S_Z_X_library = S_Z_X_library,
                                             Z_X_library = Z_X_library,
                                             family = family,
                                             v_folds = v_folds)
        
        model_list$ml_models <- ml_models
      } 
      
      # still force gcomp & ipw to be glm only
      if(any(method %in% c("gcomp", "ipw"))){
        models <- vaxstrat::fit_models(data = data, 
                                       estimand = estimand,
                                       method = method,
                                       exclusion_restriction = exclusion_restriction,
                                       cross_world = cross_world,
                                       Y_name = Y_name,
                                       Z_name = Z_name,
                                       S_name = S_name,
                                       X_name = X_name,
                                       Y_Z_X_model = Y_Z_X_model,
                                       Y_X_S1_model = Y_X_S1_model,
                                       Y_X_S0_model = Y_X_S0_model,
                                       S_X_model = S_X_model,
                                       S_Z_X_model = S_Z_X_model,
                                       Z_X_model = Z_X_model,
                                       family = family)
        
        model_list$models <- models
      }
    } else{
      # ML not specified; use glms for all 
      
      models <- vaxstrat::fit_models(data = data, 
                                     estimand = estimand,
                                     method = method,
                                     exclusion_restriction = exclusion_restriction,
                                     cross_world = cross_world,
                                     Y_name = Y_name,
                                     Z_name = Z_name,
                                     S_name = S_name,
                                     X_name = X_name,
                                     Y_Z_X_model = Y_Z_X_model,
                                     Y_X_S1_model = Y_X_S1_model,
                                     Y_X_S0_model = Y_X_S0_model,
                                     S_X_model = S_X_model,
                                     S_Z_X_model = S_Z_X_model,
                                     Z_X_model = Z_X_model,
                                     family = family)
      model_list$models <- models
    } 
  }
  
  # ----------------------------------------------------------------------------
  # 2. Bootstrap standard error & confidence intervals -------------------------
  # ----------------------------------------------------------------------------
  if(return_se == FALSE | any(method %in% c("gcomp", "ipw", "bound", "cov_adj_bound", "sens"))){
    out <- bootstrap_estimates(data = data, 
                               Y_name = Y_name,
                               Z_name = Z_name,
                               S_name = S_name,
                               X_name = X_name,
                               n_boot = n_boot, 
                               family = family,
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
                               v_folds = v_folds,
                               estimand = estimand, 
                               method = method, 
                               exclusion_restriction = exclusion_restriction,
                               cross_world = cross_world,
                               two_part_model = two_part_model,
                               effect_dir = effect_dir,
                               epsilon = epsilon,
                               return_se = return_se)
  } else{
    out <- vector("list", length = length(estimand))
    names(out) <- estimand
  }
 
 # ----------------------------------------------------------------------------
 # 3. Point estimates for effects of interest & tests
 # ----------------------------------------------------------------------------
  
  # Naturally infected --------------------------------------------------------
  
  if("nat_inf" %in% estimand){
    
    # For each exclusion restriction option (TRUE, FALSE)
    for(er in exclusion_restriction){
      
      # Exclustion restriction -- save as <estimator>_ER 
      er_suffix <- if (er) "_ER" else ""
      
        if("ipw" %in% method){
          
          estimator <- paste0("ipw", er_suffix)
          
          out$nat_inf[[estimator]]$pt_est <- do_ipw_nat_inf(data = data, models = models, S_name = S_name, Y_name = Y_name, Z_name = Z_name, exclusion_restriction = er)
          
          out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
            out$nat_inf[[estimator]]$boot_se$se_additive
          
          out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
            out$nat_inf[[estimator]]$boot_se$se_log_mult
          
          out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
          out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
          
          out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                out$nat_inf[[estimator]]$boot_se$se_additive) > stats::qnorm(1 - alpha_level/2)
          out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                            out$nat_inf[[estimator]]$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level/2)
          
          class(out$nat_inf[[estimator]]) <- estimator
        }
      
        
        for(cw in cross_world){
          
          cw_suffix <- if(cw) "_CW" else ""
          
          # can't have both false
          if(er == FALSE & cw == FALSE) next 
          
            if("gcomp" %in% method){
              
              estimator <- paste0("gcomp", er_suffix, cw_suffix)
              
              out$nat_inf[[estimator]]$pt_est <- do_gcomp_nat_inf(data = data, models = models, Z_name = Z_name, X_name = X_name, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
              
              out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                out$nat_inf[[estimator]]$boot_se$se_additive
              
              out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                out$nat_inf[[estimator]]$boot_se$se_log_mult
              
              out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
              out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
              
              out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                             out$nat_inf[[estimator]]$boot_se$se_additive) > stats::qnorm(1 - alpha_level/2)
              out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                                         out$nat_inf[[estimator]]$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level/2)
              
              class(out$nat_inf[[estimator]]) <- estimator
            }
          
            if("aipw" %in% method){
              
              estimator <- paste0("aipw", er_suffix, cw_suffix)
              
              if(ml){
                out$nat_inf[[estimator]]$pt_est <- do_aipw_nat_inf(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
              } else{
                out$nat_inf[[estimator]]$pt_est <- do_aipw_nat_inf(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
              }
              
              if(is.null(out$nat_inf[[estimator]]$boot_se)){
                # closed form SE
                
                out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                  out$nat_inf[[estimator]]$pt_est['additive_se']
                
                out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                  out$nat_inf[[estimator]]$pt_est['log_multiplicative_se']
                
                out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
                out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
                
                out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                               out$nat_inf[[estimator]]$pt_est['additive_se']) > stats::qnorm(1 - alpha_level / 2)
                out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                                           out$nat_inf[[estimator]]$pt_est['log_multiplicative_se']) > stats::qnorm(1 - alpha_level / 2)
              } else{
                # bootstrap se
                
                out$nat_inf[[estimator]]$test_stat$additive <- (out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                  out$nat_inf[[estimator]]$boot_se$se_additive
                
                out$nat_inf[[estimator]]$test_stat$mult <- (out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                  out$nat_inf[[estimator]]$boot_se$se_log_mult
                
                out$nat_inf[[estimator]]$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$additive)))
                out$nat_inf[[estimator]]$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$nat_inf[[estimator]]$test_stat$mult)))
                
                out$nat_inf[[estimator]]$reject$additive <- (abs(out$nat_inf[[estimator]]$pt_est['additive_effect'] - null_hypothesis_value) / 
                                                               out$nat_inf[[estimator]]$boot_se$se_additive) > stats::qnorm(1 - alpha_level / 2)
                out$nat_inf[[estimator]]$reject$mult <- (abs(out$nat_inf[[estimator]]$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                                           out$nat_inf[[estimator]]$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level / 2)
                
              }
              
              class(out$nat_inf[[estimator]]) <- estimator
              
            }
          }
          
        }
      
    if("tmle" %in% method){
      if(ml){
        out$nat_inf$tmle$pt_est <- do_tmle_nat_inf(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$nat_inf$tmle$pt_est <- do_tmle_nat_inf(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
      
      if(is.null(out$nat_inf$tmle$boot_se)){
        # closed form SE
        
        out$nat_inf$tmle$test_stat$additive <- (out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$pt_est['additive_se']
        
        out$nat_inf$tmle$test_stat$mult <- (out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$pt_est['log_multiplicative_se']
        
        out$nat_inf$tmle$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$nat_inf$tmle$test_stat$additive)))
        out$nat_inf$tmle$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$nat_inf$tmle$test_stat$mult)))
        
        out$nat_inf$tmle$reject$additive <- (abs(out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$nat_inf$tmle$pt_est['additive_se']) > stats::qnorm(1 - alpha_level / 2)
        out$nat_inf$tmle$reject$mult <- (abs(out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$nat_inf$tmle$pt_est['log_multiplicative_se']) > stats::qnorm(1 - alpha_level / 2)
      } else{
        # bootstrap se
        
        out$nat_inf$tmle$test_stat$additive <- (out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$boot_se$se_additive
        
        out$nat_inf$tmle$test_stat$mult <- (out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$nat_inf$tmle$boot_se$se_log_mult
        
        out$nat_inf$tmle$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$nat_inf$tmle$test_stat$additive)))
        out$nat_inf$tmle$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$nat_inf$tmle$test_stat$mult)))
        
        out$nat_inf$tmle$reject$additive <- (abs(out$nat_inf$tmle$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$nat_inf$tmle$boot_se$se_additive) > stats::qnorm(1 - alpha_level / 2)
        out$nat_inf$tmle$reject$mult <- (abs(out$nat_inf$tmle$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                               out$nat_inf$tmle$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level / 2)
      }
      
      class(out$nat_inf$tmle) <- "tmle"

    }
    
    if("sens" %in% method){
      if(ml){
        out$nat_inf$sens$pt_est <- do_sens_aipw_nat_inf(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      } else{
        out$nat_inf$sens$pt_est <- do_sens_aipw_nat_inf(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      }
      
      if(is.null(out$nat_inf$sens$boot_se)){
        # closed form SE
        
        # test stats
        test_stat_add <- (out$nat_inf$sens$pt_est$additive_effect - null_hypothesis_value) /
          out$nat_inf$sens$pt_est$additive_se
        
        test_stat_mult <- (out$nat_inf$sens$pt_est$log_multiplicative_effect - null_hypothesis_value) /
          out$nat_inf$sens$pt_est$log_multiplicative_se
        
        # p-values
        p_val_add <- 2 * (1 - stats::pnorm(abs(test_stat_add)))
        p_val_mult <- 2 * (1 - stats::pnorm(abs(test_stat_mult)))
        
        # results data frames
        out$nat_inf$sens$reject$additive <- data.frame(
          epsilon   = out$nat_inf$sens$pt_est$epsilon,
          pt_est    = out$nat_inf$sens$pt_est$additive_effect,
          se        = out$nat_inf$sens$pt_est$additive_se,
          test_stat = test_stat_add,
          p_val      = p_val_add,
          reject    = p_val_add < alpha_level
        )
        
        out$nat_inf$sens$reject$mult <- data.frame(
          epsilon   = out$nat_inf$sens$pt_est$epsilon,
          pt_est    = exp(out$nat_inf$sens$pt_est$log_multiplicative_effect),
          se        = out$nat_inf$sens$pt_est$log_multiplicative_se,
          test_stat = test_stat_mult,
          p_val      = p_val_mult,
          reject    = p_val_mult < alpha_level
        )
        
      } else{
        # bootstrap se
        
        # additive
        test_stat_add <- (out$nat_inf$sens$pt_est$additive_effect - null_hypothesis_value) /
          out$nat_inf$sens$boot_se$se_additive
        p_val_add <- 2 * (1 - stats::pnorm(abs(test_stat_add)))
        
        out$nat_inf$sens$reject$additive <- data.frame(
          epsilon = out$nat_inf$sens$boot_se$epsilon,
          pt_est = out$nat_inf$sens$pt_est$additive_effect,
          se = out$nat_inf$sens$boot_se$se_additive,
          test_stat = test_stat_add,
          p_val = p_val_add,
          reject = p_val_add < alpha_level
        )
        
        # multiplicative
        test_stat_mult <- (out$nat_inf$sens$pt_est$log_multiplicative_effect - null_hypothesis_value) /
          out$nat_inf$sens$boot_se$se_mult
        p_val_mult <- 2 * (1 - stats::pnorm(abs(test_stat$mult)))
        
        out$nat_inf$sens$reject$mult <- data.frame(
          epsilon = out$nat_inf$sens$boot_se$epsilon,
          pt_est = exp(out$nat_inf$sens$pt_est$log_multiplicative_effect),
          se = out$nat_inf$sens$boot_se$se_mult,
          test_stat = test_stat_mult,
          p_val = p_val_mult,
          reject = p_val_mult < alpha_level
        )
        
      }
      
      class(out$nat_inf$sens) <- "sens"

    }
    
    if("bound" %in% method){
      
      out$nat_inf$bound$pt_est <- get_bound_nat_inf(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
      
      # Bounds test - one sided
      # If effect direction < 0, test upper bound; Else, test lower bound 
      if(effect_dir == "negative"){
        out$nat_inf$bound$test_stat$additive <- out$nat_inf$bound$pt_est['additive_effect_upper'] / out$nat_inf$bound$boot_se$se_additive_upper
        out$nat_inf$bound$p_val$additive <- stats::pnorm(out$nat_inf$bound$test_stat$additive, lower.tail = TRUE) 
        out$nat_inf$bound$reject$additive <- out$nat_inf$bound$test_stat$additive < stats::qnorm(alpha_level)
      
        out$nat_inf$bound$test_stat$mult <- log(out$nat_inf$bound$pt_est['mult_effect_upper']) / out$nat_inf$bound$boot_se$se_log_mult_upper
        out$nat_inf$bound$p_val$mult <- stats::pnorm(out$nat_inf$bound$test_stat$mult, lower.tail = TRUE) # had note about * 2 but i think we want one-sided right? 
        out$nat_inf$bound$reject$mult <- out$nat_inf$bound$test_stat$mult < stats::qnorm(alpha_level)
      } else{
        out$nat_inf$bound$test_stat$additive <- out$nat_inf$bound$pt_est['additive_effect_lower'] / out$nat_inf$bound$boot_se$se_additive_lower
        out$nat_inf$bound$p_val$additive <- stats::pnorm(out$nat_inf$bound$test_stat$additive, lower.tail = FALSE)
        out$nat_inf$bound$reject$additive <- out$nat_inf$bound$test_stat$additive > stats::qnorm(1-alpha_level)
        
        out$nat_inf$bound$test_stat$mult <- log(out$nat_inf$bound$pt_est['mult_effect_lower']) / out$nat_inf$bound$boot_se$se_log_mult_lower
        out$nat_inf$bound$p_val$mult <- stats::pnorm(out$nat_inf$bound$test_stat$mult, lower.tail = FALSE) # had note about * 2 but i think we want one-sided right? 
        out$nat_inf$bound$reject$mult <- out$nat_inf$bound$test_stat$mult > stats::qnorm(1-alpha_level)
        
      }
      
      # Permutation test
      if(permutation){
        out$nat_inf$bound$permutation <- permutation_bound_nat_inf(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, n_permutations = n_perm, family = family, effect_dir = effect_dir)
      }
      
      class(out$nat_inf$bound) <- "bound"
      
    }
    
    
    if("cov_adj_bound" %in% method){
      
      out$nat_inf$cov_adj_bound$pt_est <- get_cov_adj_bound_nat_inf(data = data, X_name = X_name, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
    
      # Bounds test - one sided
      # If effect direction < 0, test upper bound; Else, test lower bound 
      if(effect_dir == "negative"){
        out$nat_inf$cov_adj_bound$test_stat$additive <- out$nat_inf$cov_adj_bound$pt_est['additive_effect_upper'] / out$nat_inf$cov_adj_bound$boot_se$se_additive_upper
        out$nat_inf$cov_adj_bound$p_val$additive <- stats::pnorm(out$nat_inf$cov_adj_bound$test_stat$additive, lower.tail = TRUE) 
        out$nat_inf$cov_adj_bound$reject$additive <- out$nat_inf$cov_adj_bound$test_stat$additive < stats::qnorm(alpha_level)
        
        out$nat_inf$cov_adj_bound$test_stat$mult <- log(out$nat_inf$cov_adj_bound$pt_est['mult_effect_upper']) / out$nat_inf$bound$boot_se$se_log_mult_upper
        out$nat_inf$cov_adj_bound$p_val$mult <- stats::pnorm(out$nat_inf$cov_adj_bound$test_stat$mult, lower.tail = TRUE) # had note about * 2 but i think we want one-sided right? 
        out$nat_inf$cov_adj_bound$reject$mult <- out$nat_inf$cov_adj_bound$test_stat$mult < stats::qnorm(alpha_level)
      } else{
        out$nat_inf$cov_adj_bound$test_stat$additive <- out$nat_inf$cov_adj_bound$pt_est['additive_effect_lower'] / out$nat_inf$cov_adj_bound$boot_se$se_additive_lower
        out$nat_inf$cov_adj_bound$p_val$additive <- stats::pnorm(out$nat_inf$cov_adj_bound$test_stat$additive, lower.tail = FALSE)
        out$nat_inf$cov_adj_bound$reject$additive <- out$nat_inf$cov_adj_bound$test_stat$additive > stats::qnorm(1-alpha_level)
        
        out$nat_inf$cov_adj_bound$test_stat$mult <- log(out$nat_inf$cov_adj_bound$pt_est['mult_effect_lower']) / out$nat_inf$cov_adj_bound$boot_se$se_log_mult_lower
        out$nat_inf$cov_adj_bound$p_val$mult <- stats::pnorm(out$nat_inf$cov_adj_bound$test_stat$mult, lower.tail = FALSE) # had note about * 2 but i think we want one-sided right? 
        out$nat_inf$cov_adj_bound$reject$mult <- out$nat_inf$cov_adj_bound$test_stat$mult > stats::qnorm(1-alpha_level)
        
      }
   
      # Permutation test not included for covariate adjusted bound
      
      class(out$nat_inf$cov_adj_bound) <- "cov_adj_bound"
      
    }
    class(out$nat_inf) <- "nat_inf"
    
  }
  
  # Doomed --------------------------------------------------------------------
  
  if("doomed" %in% estimand){
    
    if("gcomp" %in% method){
      out$doomed$gcomp$pt_est <- do_gcomp_doomed(data = data, models = models)
      
      out$doomed$gcomp$test_stat$additive <- (out$doomed$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$doomed$gcomp$boot_se$se_additive
      
      out$doomed$gcomp$test_stat$mult <- (out$doomed$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$doomed$gcomp$boot_se$se_log_mult
      
      out$doomed$gcomp$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$doomed$gcomp$test_stat$additive)))
      out$doomed$gcomp$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$doomed$gcomp$test_stat$mult)))
      
      out$doomed$gcomp$reject$additive <- (abs(out$doomed$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
                                              out$doomed$gcomp$boot_se$se_additive) > stats::qnorm(1 - alpha_level/2)
      out$doomed$gcomp$reject$mult <- (abs(out$doomed$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                          out$doomed$gcomp$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level/2)
      
      class(out$doomed$gcomp) <- "gcomp"
    }
    
    if("ipw" %in% method){
      out$doomed$ipw$pt_est <- do_ipw_doomed(data = data, models = models, S_name = S_name, Y_name = Y_name, Z_name = Z_name)
      
      out$doomed$ipw$test_stat$additive <- (out$doomed$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$doomed$ipw$boot_se$se_additive
      
      out$doomed$ipw$test_stat$mult <- (out$doomed$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$doomed$ipw$boot_se$se_log_mult
      
      out$doomed$ipw$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$doomed$ipw$test_stat$additive)))
      out$doomed$ipw$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$doomed$ipw$test_stat$mult)))
      
      out$doomed$ipw$reject$additive <- (abs(out$doomed$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                            out$doomed$ipw$boot_se$se_additive) > stats::qnorm(1 - alpha_level/2)
      out$doomed$ipw$reject$mult <- (abs(out$doomed$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                        out$doomed$ipw$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level/2)
      
      class(out$doomed$ipw) <- "ipw"
    }
    
    if("aipw" %in% method){
      if(ml){
        out$doomed$aipw$pt_est <- do_aipw_doomed(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$doomed$aipw$pt_est <- do_aipw_doomed(data = data, models = models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
      
      if(is.null(out$doomed$aipw$boot_se)){
        # closed form SE
        
        out$doomed$aipw$test_stat$additive <- (out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$pt_est['additive_se']
        
        out$doomed$aipw$test_stat$mult <- (out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$pt_est['log_multiplicative_se']
        
        out$doomed$aipw$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$doomed$aipw$test_stat$additive)))
        out$doomed$aipw$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$doomed$aipw$test_stat$mult)))
        
        out$doomed$aipw$reject$additive <- (abs(out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$doomed$aipw$pt_est['additive_se']) > stats::qnorm(1 - alpha_level / 2)
        out$doomed$aipw$reject$mult <- (abs(out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$doomed$aipw$pt_est['log_multiplicative_se']) > stats::qnorm(1 - alpha_level / 2)
      } else{
        # bootstrap se
        
        out$doomed$aipw$test_stat$additive <- (out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$boot_se$se_additive
        
        out$doomed$aipw$test_stat$mult <- (out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$doomed$aipw$boot_se$se_log_mult
        
        out$doomed$aipw$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$doomed$aipw$test_stat$additive)))
        out$doomed$aipw$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$doomed$aipw$test_stat$mult)))
        
        out$doomed$aipw$reject$additive <- (abs(out$doomed$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$doomed$aipw$boot_se$se_additive) > stats::qnorm(1 - alpha_level / 2)
        out$doomed$aipw$reject$mult <- (abs(out$doomed$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$doomed$aipw$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level / 2)
        
      }
      
      class(out$doomed$aipw) <- "aipw"
      
    }
    
    if("bound" %in% method){
      
      out$doomed$bound$pt_est <- get_bound_doomed(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
      
      # Bounds test - one sided
      # If effect direction < 0, test upper bound; Else, test lower bound 
      if(effect_dir == "negative"){
        out$doomed$bound$test_stat$additive <- out$doomed$bound$pt_est['additive_effect_upper'] / out$doomed$bound$boot_se$se_additive_upper
        out$doomed$bound$p_val$additive <- stats::pnorm(out$doomed$bound$test_stat$additive, lower.tail = TRUE) 
        out$doomed$bound$reject$additive <- out$doomed$bound$test_stat$additive < stats::qnorm(alpha_level)
        
        out$doomed$bound$test_stat$mult <- log(out$doomed$bound$pt_est['mult_effect_upper']) / out$doomed$bound$boot_se$se_log_mult_upper
        out$doomed$bound$p_val$mult <- stats::pnorm(out$doomed$bound$test_stat$mult, lower.tail = TRUE) 
        out$doomed$bound$reject$mult <- out$doomed$bound$test_stat$mult < stats::qnorm(alpha_level)
        
      } else{
        out$doomed$bound$test_stat$additive <- out$doomed$bound$pt_est['additive_effect_lower'] / out$doomed$bound$boot_se$se_additive_lower
        out$doomed$bound$p_val$additive <- stats::pnorm(out$doomed$bound$test_stat$additive, lower.tail = FALSE)
        out$doomed$bound$reject$additive <- out$doomed$bound$test_stat$additive > stats::qnorm(1-alpha_level)
        
        out$doomed$bound$test_stat$mult <- log(out$doomed$bound$pt_est['mult_effect_lower']) / out$doomed$bound$boot_se$se_log_mult_lower
        out$doomed$bound$p_val$mult <- stats::pnorm(out$doomed$bound$test_stat$mult, lower.tail = FALSE) 
        out$doomed$bound$reject$mult <- out$doomed$bound$test_stat$mult > stats::qnorm(1-alpha_level)
      }
      
      # Permutation test
      if(permutation){
        out$doomed$bound$permutation <- permutation_bound_doomed(data = data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, n_permutations = n_perm, family = family, effect_dir = effect_dir)
      }
      
      class(out$doomed$bound) <- "bound"
      
    }
    
    if(!is.null(out$doomed)){
      class(out$doomed) <- "doomed"
    }

  }
  
  # Population ----------------------------------------------------------------
  
  if("pop" %in% estimand){
    
    if("gcomp" %in% method){
      out$pop$gcomp$pt_est <- do_gcomp_pop(data = data, models = models,  Z_name = Z_name, X_name = X_name, two_part_model = two_part_model)
      
      out$pop$gcomp$test_stat$additive <- (out$pop$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$pop$gcomp$boot_se$se_additive
      
      out$pop$gcomp$test_stat$mult <- (out$pop$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$pop$gcomp$boot_se$se_log_mult
      
      out$pop$gcomp$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$pop$gcomp$test_stat$additive)))
      out$pop$gcomp$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$pop$gcomp$test_stat$mult)))
      
      out$pop$gcomp$reject$additive <- (abs(out$pop$gcomp$pt_est['additive_effect'] - null_hypothesis_value) / 
                                              out$pop$gcomp$boot_se$se_additive) > stats::qnorm(1 - alpha_level/2)
      out$pop$gcomp$reject$mult <- (abs(out$pop$gcomp$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                          out$pop$gcomp$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level/2)
      
      class(out$pop$gcomp) <- "gcomp"
    }
    
    if("ipw" %in% method){
      out$pop$ipw$pt_est <- do_ipw_pop(data = data, models = models, Y_name = Y_name, Z_name = Z_name)
      
      out$pop$ipw$test_stat$additive <- (out$pop$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
        out$pop$ipw$boot_se$se_additive
      
      out$pop$ipw$test_stat$mult <- (out$pop$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
        out$pop$ipw$boot_se$se_log_mult
      
      out$pop$ipw$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$pop$ipw$test_stat$additive)))
      out$pop$ipw$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$pop$ipw$test_stat$mult)))
      
      out$pop$ipw$reject$additive <- (abs(out$pop$ipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                            out$pop$ipw$boot_se$se_additive) > stats::qnorm(1 - alpha_level/2)
      out$pop$ipw$reject$mult <- (abs(out$pop$ipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                        out$pop$ipw$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level/2)
      
      class(out$pop$ipw) <- "ipw"
    }
    
    if("aipw" %in% method){
      if(ml){
        out$pop$aipw$pt_est <- do_aipw_pop(data = data, models = ml_models, Y_name = Y_name, Z_name = Z_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      } else{
        out$pop$aipw$pt_est <- do_aipw_pop(data = data, models = models, Y_name = Y_name, Z_name = Z_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      }
      
      if(is.null(out$pop$aipw$boot_se)){
        # closed form SE
        
        out$pop$aipw$test_stat$additive <- (out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$pop$aipw$pt_est['additive_se']
        
        out$pop$aipw$test_stat$mult <- (out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$pop$aipw$pt_est['log_multiplicative_se']
        
        out$pop$aipw$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$pop$aipw$test_stat$additive)))
        out$pop$aipw$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$pop$aipw$test_stat$mult)))
        
        out$pop$aipw$reject$additive <- (abs(out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$pop$aipw$pt_est['additive_se']) > stats::qnorm(1 - alpha_level / 2)
        out$pop$aipw$reject$mult <- (abs(out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$pop$aipw$pt_est['log_multiplicative_se']) > stats::qnorm(1 - alpha_level / 2)
      } else{
        # bootstrap se
        
        out$pop$aipw$test_stat$additive <- (out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
          out$pop$aipw$boot_se$se_additive
        
        out$pop$aipw$test_stat$mult <- (out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
          out$pop$aipw$boot_se$se_log_mult
        
        out$pop$aipw$p_val$additive <- 2 * (1 - stats::pnorm(abs(out$pop$aipw$test_stat$additive)))
        out$pop$aipw$p_val$mult <- 2 * (1 - stats::pnorm(abs(out$pop$aipw$test_stat$mult)))
        
        out$pop$aipw$reject$additive <- (abs(out$pop$aipw$pt_est['additive_effect'] - null_hypothesis_value) / 
                                               out$pop$aipw$boot_se$se_additive) > stats::qnorm(1 - alpha_level / 2)
        out$pop$aipw$reject$mult <- (abs(out$pop$aipw$pt_est['log_multiplicative_effect'] - null_hypothesis_value) / 
                                           out$pop$aipw$boot_se$se_log_mult) > stats::qnorm(1 - alpha_level / 2)
        
      }
      
      class(out$pop$aipw) <- "aipw"
      
    }
    
    if(!is.null(out$pop)){
      class(out$pop) <- "pop"
    }
    
  }
  
  if(return_models){
    out$models <- model_list
  }
 
  class(out) <- "vaxstrat"
  
  return(out)

}
