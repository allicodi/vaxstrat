#' Single bootstrap replicate for vaxstrat estimators
#'
#' Generates one bootstrap sample from the observed data and computes
#' specified estimators of growth effects within principal strata
#' (e.g., naturally infected, doomed, and population).
#'
#' @param data A data.frame containing the observed data.
#' @param Y_name Character string specifying the outcome (growth) variable.
#' @param Z_name Character string specifying the treatment (vaccination) variable.
#' @param X_name Character vector of covariate names.
#' @param S_name Character string specifying the infection indicator.
#' @param estimand Character vector indicating which estimands to compute.
#'   Options include `"nat_inf"`, `"doomed"`, and `"pop"`.
#' @param method Character vector of estimation methods to use. Options include
#'   `"gcomp"`, `"ipw"`, `"aipw"`, `"tmle"`, `"bound"`, `"cov_adj_bound"`, and `"sens"`.
#' @param exclusion_restriction Logical vector indicating whether to impose the
#'   exclusion restriction assumption. Results are returned separately for each value.
#' @param cross_world Logical vector indicating whether to impose the cross-world
#'   assumption. Results are returned separately.
#' @param ml Logical; if TRUE, uses SuperLearner-based models. If FALSE (default),
#'   parametric models (e.g., GLMs) are used.
#' @param Y_Z_X_model Optional model specification for outcome regression on treatment
#'   and covariates.
#' @param Y_X_S1_model Optional model for outcome regression among infected individuals.
#' @param Y_X_S0_model Optional model for outcome regression among uninfected individuals.
#' @param S_X_model Optional model for infection regression on covariates.
#' @param S_Z_X_model Optional model for infection regression on treatment and covariates.
#' @param Z_X_model Optional model for treatment regression on covariates.
#' @param Y_Z_X_library SuperLearner libraries for outcome regression on treatment + covariates.
#' @param Y_X_library SuperLearner libraries for outcome regression on covariates.
#' @param S_X_library SuperLearner libraries for infection regression.
#' @param S_Z_X_library SuperLearner libraries for infection regression on treatment + covariates.
#' @param Z_X_library SuperLearner libraries for treatment regression.
#' @param family Outcome family (default `"gaussian"` for continuous growth outcomes).
#' @param v_folds Number of cross-validation folds for SuperLearner (default 3).
#' @param effect_dir Direction of beneficial effect (`"positive"` or `"negative"`).
#'   Used for one-sided inference in bound-based methods.
#' @param epsilon Numeric vector of sensitivity parameters for sensitivity analysis.
#' @param max_resample Maximum number of attempts to resample if the bootstrap sample
#'   violates the monotonicity condition (i.e., infection risk in vaccinated exceeds
#'   unvaccinated).
#' @param return_se Logical; if FALSE, returns influence-function–based standard errors
#'   (if available). If TRUE, standard errors are obtained via bootstrap aggregation.
#' @param two_part_model Logical; whether to use a two-part outcome model.
#'
#' @details
#' A bootstrap sample is drawn with replacement from the observed data.
#' To enforce a monotonicity-type condition, the sample is resampled up to
#' `max_resample` times if the infection rate among vaccinated individuals
#' exceeds that among unvaccinated individuals.
#'
#' Estimators are computed separately for each requested estimand and method.
#' Results are returned as a nested list indexed by estimand and estimator name.
#'
#' @returns A nested list with elements corresponding to each estimand
#' (`"nat_inf"`, `"doomed"`, `"pop"`). Each contains sublists for each estimator,
#' including point estimates, test statistics, p-values, and (optionally)
#' standard errors or sensitivity analysis results.
#'
#' @export
one_boot <- function(
    data,
    Y_name = "Y",
    Z_name = "Z",
    X_name = "X",
    S_name = "S", 
    estimand = c("nat_inf", "doomed", "pop"),
    method = c("gcomp", "ipw", "aipw", "tmle", "bound", "cov_adj_bound", "sens"),
    exclusion_restriction = FALSE,
    cross_world = TRUE,
    ml = FALSE, 
    Y_Z_X_model = NULL,
    Y_X_S1_model = NULL,
    Y_X_S0_model = NULL,
    S_X_model = NULL,
    S_Z_X_model = NULL,
    Z_X_model = NULL,
    Y_Z_X_library = c("SL.glm"),
    Y_X_library = c("SL.glm"),
    S_X_library = c("SL.glm"),
    S_Z_X_library = c("SL.glm"),
    Z_X_library = c("SL.mean"),
    family = "gaussian",
    v_folds = 3,
    effect_dir = "positive",
    epsilon = exp(seq(log(0.5), log(2), length = 50)),
    max_resample = 10,
    return_se = TRUE,
    two_part_model = FALSE
){
  
  n <- dim(data)[1]
  boot_row_idx <- sample(1:n, replace=TRUE)
  boot_data <- data[boot_row_idx,]
  
  # If there are more infections in vaccine arm than placebo arm, resample up to max_resample 
  rhobar_0_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 0])
  rhobar_1_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 1])
  resample <- 0
  
  while(rhobar_0_n <= rhobar_1_n & resample <= max_resample){
    boot_row_idx <- sample(1:n, replace=TRUE)
    boot_data <- data[boot_row_idx,]
    
    rhobar_0_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 0])
    rhobar_1_n <- mean(boot_data[[S_name]][boot_data[[Z_name]] == 1])
    resample <- resample + 1
  }
  
  if(resample > max_resample){
    stop(paste0("Exceeded max_resample of ", max_resample, "for given bootstrap replicate"))
  }
  
  # compute estimators using bootstrap data set
  if(ml){
    
    if(any(method %in% c("aipw", "tmle", "sens")) & return_se == FALSE){
      boot_ml_models <- vaxstrat::fit_ml_models(data = boot_data, 
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
    } 
    
    if(any(method %in% c("gcomp", "ipw"))){
      boot_models <- vaxstrat::fit_models(data = boot_data, 
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
    }
    
  } else{
    # GLMS for all
    if(any(method %in% c("gcomp", "ipw", "aipw", "tmle", "sens"))){
      boot_models <- vaxstrat::fit_models(data = boot_data, 
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
    }
    
    # Otherwise for bounds only no models needed
    
  } 
  
  out <- vector("list", length = length(estimand))
  names(out) <- estimand
  
  # Naturally infected --------------------------------------------------------
  
  if("nat_inf" %in% estimand){
    
    for(er in exclusion_restriction){
      
      # Exclustion restriction -- save as <estimator>_ER 
      er_suffix <- if (er) "_ER" else ""
      
      if("ipw" %in% method){
        estimator <- paste0("ipw", er_suffix)
        out$nat_inf[[estimator]] <- do_ipw_nat_inf(data = boot_data, models = boot_models, S_name = S_name, Y_name = Y_name, Z_name = Z_name, exclusion_restriction = er)
      }
      
      # Cross-world assumption can be toggled for AIPW only at this point
      for(cw in cross_world){
        
        cw_suffix <- if(cw) "_CW" else ""
        
        # cannot have scenario where both are false
        if(er == FALSE & cw == FALSE) next
        
        if("gcomp" %in% method){
          estimator <- paste0("gcomp", er_suffix, cw_suffix)
          out$nat_inf[[estimator]] <- do_gcomp_nat_inf(data = boot_data, models = boot_models, Z_name = Z_name, X_name = X_name, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
        }
        
        # if we want bootstrap SE for AIPW (otherwise return closed form SE when we get point estimate)
        if("aipw" %in% method & return_se == FALSE){
          estimator <- paste0("aipw", er_suffix, cw_suffix)
          if(ml){
            out$nat_inf[[estimator]] <- do_aipw_nat_inf(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
          } else{
            out$nat_inf[[estimator]] <- do_aipw_nat_inf(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, X_name = X_name, return_se = return_se, exclusion_restriction = er, cross_world = cw, two_part_model = two_part_model)
          }
        }
      }
      
      
    }
    
    if("tmle" %in% method & return_se == FALSE){
      if(ml){
        out$nat_inf$tmle <- do_tmle_nat_inf(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$nat_inf$tmle <- do_tmle_nat_inf(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
    }
    
    if("sens" %in% method & return_se == FALSE){
      if(ml){
        out$nat_inf$sens<- do_sens_aipw_nat_inf(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      } else{
        out$nat_inf$sens <- do_sens_aipw_nat_inf(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, epsilon = epsilon, return_se = return_se)
      }
    }
    
    if("bound" %in% method){
      out$nat_inf$bound <- get_bound_nat_inf(data = boot_data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
    }
    
    if("cov_adj_bound" %in% method){
      if(length(X_name) > 1) stop("cov_adj_bound only implemented for single covariate")
      out$nat_inf$cov_adj_bound <- get_cov_adj_bound_nat_inf(data = boot_data, X_name = X_name, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
    }
    
  }
  
  # Doomed --------------------------------------------------------------------
  
  if("doomed" %in% estimand){
    
    if("gcomp" %in% method){
      out$doomed$gcomp <- do_gcomp_doomed(data = boot_data, models = boot_models)
    }
    
    if("ipw" %in% method){
      out$doomed$ipw <- do_ipw_doomed(data = boot_data, models = boot_models, S_name = S_name, Y_name = Y_name, Z_name = Z_name)
    }
    
    if("aipw" %in% method & return_se == FALSE){
      if(ml){
        out$doomed$aipw <- do_aipw_doomed(data = boot_data, models = boot_ml_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      } else{
        out$doomed$aipw <- do_aipw_doomed(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name, S_name = S_name, return_se = return_se)
      }
    }
    
    if("bound" %in% method){
      out$doomed$bound <- get_bound_doomed(data = boot_data, Y_name = Y_name, Z_name = Z_name, S_name = S_name, family = family)
    }
    
  }
  
  # Population ----------------------------------------------------------------
  
  if("pop" %in% estimand){
    
    if("gcomp" %in% method){
      out$pop$gcomp <- do_gcomp_pop(data = boot_data, models = boot_models, Z_name = Z_name, X_name = X_name)
    }
    
    if("ipw" %in% method){
      out$pop$ipw <- do_ipw_pop(data = boot_data, models = boot_models, Y_name = Y_name, Z_name = Z_name)
    }
    
    if("aipw" %in% method & return_se == FALSE){
      if(ml){
        out$pop$aipw <- do_aipw_pop(data = boot_data, models = boot_ml_models, Z_name = Z_name, Y_name = Y_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      } else{
        out$pop$aipw <- do_aipw_pop(data = boot_data, models = boot_models, Z_name = Z_name, Y_name = Y_name, X_name = X_name, return_se = return_se, two_part_model = two_part_model)
      }
    }
    
  }

  return(out)

}

#' Bootstrap standard errors and confidence intervals for vaxstrat estimators
#'
#' Repeats the bootstrap procedure using [one_boot()] to compute bootstrap-based
#' standard errors and confidence intervals for specified estimators.
#'
#' @param data A data.frame containing the observed data.
#' @param Y_name Character string specifying the outcome (growth) variable.
#' @param Z_name Character string specifying the treatment (vaccination) variable.
#' @param X_name Character vector of covariate names.
#' @param S_name Character string specifying the infection indicator.
#' @param n_boot Integer; number of bootstrap replicates (default 1000).
#' @param estimand Character vector indicating which estimands to compute.
#' @param method Character vector of estimation methods to use.
#' @param exclusion_restriction Logical vector indicating whether to impose the
#'   exclusion restriction assumption.
#' @param cross_world Logical vector indicating whether to impose the cross-world assumption.
#' @param ml Logical; if TRUE, uses SuperLearner models.
#' @param Y_Z_X_model Optional model specification for outcome regression.
#' @param Y_X_S1_model Optional model for outcome among infected.
#' @param Y_X_S0_model Optional model for outcome among uninfected.
#' @param S_X_model Optional model for infection regression.
#' @param S_Z_X_model Optional model for infection regression on treatment + covariates.
#' @param Z_X_model Optional model for treatment regression.
#' @param Y_Z_X_library SuperLearner libraries for outcome regression.
#' @param Y_X_library SuperLearner libraries for covariate adjustment.
#' @param S_X_library SuperLearner libraries for infection regression.
#' @param S_Z_X_library SuperLearner libraries for infection regression on treatment + covariates.
#' @param Z_X_library SuperLearner libraries for treatment regression.
#' @param family Outcome family (default `"gaussian"`).
#' @param v_folds Number of cross-validation folds for SuperLearner.
#' @param effect_dir Direction of beneficial effect.
#' @param epsilon Numeric vector of sensitivity parameters.
#' @param return_se Logical; if FALSE, bootstrap is used to compute standard errors.
#' @param two_part_model Logical; whether to use a two-part outcome model.
#'
#' @details
#' This function calls [one_boot()] repeatedly (`n_boot` times) to generate
#' bootstrap replicates of each estimator. Bootstrap standard errors and
#' confidence intervals are computed using the empirical distribution of
#' bootstrap estimates.
#'
#' Results are aggregated separately by estimand and method. For sensitivity
#' analyses and bound estimators, specialized bootstrap aggregation functions
#' are used.
#'
#' @returns A nested list with the same structure as a single call to
#' [one_boot()], but augmented with bootstrap-based standard errors and
#' confidence intervals (stored under `boot_se`).
#'
#' @export
bootstrap_estimates <- function(
    data, 
    Y_name = "Y",
    Z_name = "Z",
    X_name = "X",
    S_name = "S", 
    n_boot = 1000, 
    estimand = c("nat_inf", "doomed", "pop"),
    method = c("gcomp", "ipw", "aipw", "tmle", "bound", "sens"),
    exclusion_restriction = FALSE,
    cross_world = TRUE,
    ml = ml,
    Y_Z_X_model = NULL,
    Y_X_S1_model = NULL, 
    Y_X_S0_model = NULL, 
    S_X_model = NULL,
    S_Z_X_model = NULL,
    Z_X_model = NULL,
    Y_Z_X_library = c("SL.glm"),
    Y_X_library = c("SL.glm"),
    S_X_library = c("SL.glm"),
    S_Z_X_library = c("SL.glm"),
    Z_X_library = c("SL.mean"),
    family = "gaussian",
    v_folds = 3,
    effect_dir = "positive",
    epsilon = exp(seq(log(0.5), log(2), length = 50)),
    return_se = TRUE,
    two_part_model = FALSE
){
  
  # Initial boot_estimates for all viable estimand & method combinations
  boot_estimates <- replicate(n_boot, one_boot(data, 
                                               Y_name = Y_name,
                                               Z_name = Z_name,
                                               S_name = S_name,
                                               X_name = X_name,
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
                                               v_folds = v_folds,
                                               family = family,
                                               epsilon = epsilon,
                                               return_se = return_se,
                                               two_part_model = two_part_model), simplify = FALSE)
    # List to store results
    out <- vector("list", length = length(estimand))
    names(out) <- estimand
  
    # Naturally infected --------------------------------------------------------
    
    if("nat_inf" %in% estimand){
      
      for(er in exclusion_restriction){
        
        # Exclustion restriction -- save as <estimator>_ER 
        er_suffix <- if (er) "_ER" else ""
        
        if("ipw" %in% method){
          estimator <- paste0("ipw", er_suffix)
          out$nat_inf[[estimator]]$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = estimator)
        }
        
        for(cw in cross_world){
          
          # can't have both false
          if(er == FALSE & cw == FALSE) next
          
          cw_suffix <- if (cw) "_CW" else ""
          
          if("gcomp" %in% method){
            estimator <- paste0("gcomp", er_suffix, cw_suffix)
            out$nat_inf[[estimator]]$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = estimator)
          }
          
          # if we want bootstrap SE for AIPW (otherwise return closed form SE when we get point estimate)
          if("aipw" %in% method & return_se == FALSE){
            estimator <- paste0("aipw", er_suffix, cw_suffix)
            out$nat_inf[[estimator]]$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = estimator)
          }
          
        }
        
      }

      if("tmle" %in% method & return_se == FALSE){
        out$nat_inf$tmle$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "nat_inf", method = "tmle")
      }
      
      if("sens" %in% method & return_se == FALSE){
        out$nat_inf$sens$boot_se <- get_boot_se_sens(boot_estimates = boot_estimates, estimand = "nat_inf", method = "sens")
      }
      
      if("bound" %in% method){
        out$nat_inf$bound$boot_se <- get_boot_se_bound(boot_estimates = boot_estimates, estimand = "nat_inf", method = "bound")
      }

      if("cov_adj_bound" %in% method){
        out$nat_inf$cov_adj_bound$boot_se <- get_boot_se_bound(boot_estimates = boot_estimates, estimand = "nat_inf", method = "cov_adj_bound")
      }
      
    }
    
    
    # Doomed --------------------------------------------------------------------
    
    if("doomed" %in% estimand){
      
      if("gcomp" %in% method){
        out$doomed$gcomp$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "doomed", method = "gcomp")
      }
      
      if("ipw" %in% method){
        out$doomed$ipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "doomed", method = "ipw")
      }
      
      if("aipw" %in% method & return_se == FALSE){
        out$doomed$aipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "doomed", method = "aipw")
      }
      
      if("bound" %in% method){
        out$doomed$bound$boot_se <- get_boot_se_bound(boot_estimates = boot_estimates, estimand = "doomed", method = "bound")
      }
      
    }
    
    # Population ----------------------------------------------------------------
    
    if("pop" %in% estimand){
      
      if("gcomp" %in% method){
        out$pop$gcomp$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "pop", method = "gcomp")
      }
      
      if("ipw" %in% method){
        out$pop$ipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "pop", method = "ipw")
      }
      
      if("aipw" %in% method & return_se == FALSE){
        out$pop$aipw$boot_se <- get_boot_se(boot_estimates = boot_estimates, estimand = "pop", method = "aipw")
      }
      
    }

  return(out)
    
}