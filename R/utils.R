#' Unified prediction helper for fitted nuisance models
#'
#' Generates predictions from supported fitted model objects in a consistent
#' way. Supported model classes are \code{"glm"}, \code{"SuperLearner"}, and
#' \code{"cross_fit"}.
#'
#' @param model A fitted model object of class \code{"glm"},
#'   \code{"SuperLearner"}, or \code{"cross_fit"}.
#' @param newdata A \code{data.frame} containing covariate values at which
#'   predictions are computed.
#'
#' @details
#' This function standardizes prediction calls across supported model classes:
#' \itemize{
#'   \item For \code{"glm"} objects, predictions are obtained using
#'   \code{type = "response"}.
#'   \item For \code{"SuperLearner"} objects, predictions are obtained from
#'   \code{predict.SuperLearner()} and extracted from the \code{$pred} slot.
#'   \item For \code{"cross_fit"} objects, predictions are obtained using
#'   \code{\link{predict.cross_fit}}, yielding out-of-fold predictions.
#' }
#'
#' @return A numeric vector of predicted values corresponding to \code{newdata}.
#'
#' @seealso \code{\link{predict.cross_fit}}
#'
#' @export
simple_predict <- function(model, newdata){
  if(!any(c("glm", "SuperLearner", "cross_fit") %in% class(model))){
    stop("method only guaranteed to work correctly with cross_fit, glm, or SuperLearner")
  }
  if( inherits(model, "SuperLearner") ){
    pred <- predict(model, newdata = newdata)$pred
  }else{
    pred <- predict(model, newdata = newdata, type = "response")
  }
  return(pred)
}


#' Obtain out-of-fold predictions from a cross-fitted nuisance model
#'
#' Generates predictions from an object of class \code{"cross_fit"} by applying
#' the fitted model for each fold to the observations assigned to that fold. This
#' ensures that each prediction is obtained from a model that was not trained on
#' the corresponding observation.
#'
#' @param object An object of class \code{"cross_fit"}, typically one nuisance
#'   model component returned by \code{\link{cross_fit_models}} after
#'   reorganization by \code{\link{reorg_cross_fit_models}}.
#' @param newdata A \code{data.frame} containing covariate values at which
#'   predictions are computed.
#' @param fold_assignments Optional \code{data.frame} containing fold
#'   assignments. If \code{NULL}, the fold assignments stored in \code{object}
#'   are used.
#' @param id_name Optional character string giving the identifier column in
#'   \code{newdata}. If \code{NULL}, the identifier stored in \code{object} is
#'   used. If the identifier is not present in \code{newdata}, predictions fall
#'   back to matching on the original row index \code{idx}; this requires
#'   \code{newdata} to have the same row order and number of rows as the data
#'   used to fit the cross-fitted models.
#' @param ... Additional arguments. Currently unused.
#'
#' @details
#' For each cross-fitting fold, this method identifies the observations in
#' \code{newdata} assigned to that fold and predicts using the corresponding
#' fold-specific fitted model. Predictions are returned in the same row order as
#' \code{newdata}.
#'
#' When \code{newdata} contains the original ID column, fold membership is
#' determined using that ID. When \code{newdata} does not contain the original ID
#' column, the function assumes that \code{newdata} is aligned with the original
#' analysis data by row order and uses \code{idx} instead. This is useful for
#' counterfactual prediction datasets, such as copies of the original data with
#' treatment set to \eqn{Z = 1} or \eqn{Z = 0}, where the ID column may have been
#' dropped.
#'
#' @return A numeric vector of predictions with length \code{nrow(newdata)}.
#'
#' @seealso \code{\link{cross_fit_models}},
#'   \code{\link{reorg_cross_fit_models}}, \code{\link{simple_predict}}
#'
#' @method predict cross_fit
#' @export
predict.cross_fit <- function(object, newdata, fold_assignments = NULL, id_name = NULL, ...) {
  
  # Use fold assignments stored in object unless supplied
  if (is.null(fold_assignments)) {
    fold_assignments <- object$fold_assignments
  }
  
  # Default to id_name stored in object
  if (is.null(id_name)) {
    id_name <- object$id_name
  }
  
  # Decide whether to match prediction rows by id or idx
  if (!is.null(id_name) && id_name %in% colnames(newdata)) {
    
    # Standard case: newdata has the original id column
    pred_id_name <- id_name
    fold_id_name <- "id"
    
  } else {
    
    # Fallback case: newdata does not have id_name, so assume rows correspond
    # to the original data order and match using idx
    newdata$idx <- seq_len(nrow(newdata))
    pred_id_name <- "idx"
    fold_id_name <- "idx"
    
    if (nrow(newdata) != nrow(fold_assignments)) {
      stop(
        "newdata does not contain id_name, so prediction falls back to idx. ",
        "This only works when newdata has the same number of rows and order as the data used for cross-fitting."
      )
    }
  }
  
  preds <- rep(NA_real_, nrow(newdata))
  
  for (k in seq_along(object$fold_specific_models)) {
    
    validation_keys <- fold_assignments[[fold_id_name]][
      fold_assignments$validation_fold == k
    ]
    
    validation_rows <- newdata[[pred_id_name]] %in% validation_keys
    
    if (!any(validation_rows)) {
      next
    }
    
    validation_data <- newdata[validation_rows, , drop = FALSE]
    
    preds[validation_rows] <- simple_predict(
      object$fold_specific_models[[k]],
      newdata = validation_data
    )
  }
  
  if (any(is.na(preds))) {
    stop("Some rows did not receive predictions. Check id_name, idx, and fold_assignments.")
  }
  
  return(preds)
}

#' Convert bootstrap results into a tidy data frame
#'
#' Transforms a list of bootstrap outputs into a single data frame for downstream
#' summarization and inference.
#'
#' @param boot_estimates A list of bootstrap results, where each element corresponds
#'   to one bootstrap replicate (e.g., output from \code{one_boot}).
#' @param estimand Character string specifying the estimand of interest
#'   (e.g., \code{"nat_inf"}, \code{"pop"}, \code{"doomed"}).
#' @param method Character string specifying the estimation method
#'   (e.g., \code{"gcomp"}, \code{"ipw"}, \code{"aipw"}, \code{"bound"}, \code{"sens_cw"}).
#'
#' @details
#' Each bootstrap replicate is expected to be a nested list indexed by
#' \code{estimand} and \code{method}. Failed or missing replicates are skipped.
#'
#' @return A data.frame containing stacked bootstrap results, with additional columns:
#' \describe{
#'   \item{boot_id}{Bootstrap replicate index.}
#'   \item{estimand}{Estimand label.}
#'   \item{method}{Estimation method label.}
#' }
#'
make_boot_df <- function(boot_estimates, estimand = "nat_inf", method = "gcomp") {
  res_list <- lapply(seq_along(boot_estimates), function(i) {
    boot_res <- boot_estimates[[i]]
    
    # Check if the desired estimand/method exists in this iteration
    if (!is.null(boot_res[[estimand]]) && !is.null(boot_res[[estimand]][[method]])) {
      row <- boot_res[[estimand]][[method]]
      
      # If it's just a single number, convert to named data frame
      if (is.atomic(row) && length(row) == 1) {
        row <- data.frame(estimate = row)
      } else if(class(row) == "sens_cw"){
        class(row) <- "data.frame"
      } else if (is.null(dim(row))) {
        row <- as.data.frame(t(row))
      }
      
      row$boot_id <- i
      row$estimand <- estimand
      row$method <- method
      return(row)
    } else {
      return(NULL)  # Skip missing or failed bootstraps
    }
  })
  # Combine all into one data frame
  data.frame(do.call(rbind, res_list))
}

#' Bootstrap standard errors and confidence intervals (point estimates)
#'
#' Computes bootstrap-based standard errors and percentile confidence intervals
#' for additive and multiplicative effects.
#'
#' @param boot_estimates A list of bootstrap results.
#' @param estimand Character string specifying the estimand of interest.
#' @param method Character string specifying the estimation method.
#'
#' @return A data.frame with:
#' \describe{
#'   \item{se_additive}{Standard error of the additive effect.}
#'   \item{lower_ci_additive}{2.5th percentile of additive effect.}
#'   \item{upper_ci_additive}{97.5th percentile of additive effect.}
#'   \item{se_log_mult}{Standard error of the log multiplicative effect.}
#'   \item{lower_ci_mult}{Lower bound of multiplicative effect (exponentiated).}
#'   \item{upper_ci_mult}{Upper bound of multiplicative effect (exponentiated).}
#'   \item{se_psi_1}{Standard error of \eqn{\psi_1}.}
#'   \item{se_psi_0}{Standard error of \eqn{\psi_0}.}
#' }
#'
get_boot_se <- function(boot_estimates, estimand = "nat_inf", method = "gcomp"){
  boot_df <- make_boot_df(boot_estimates = boot_estimates,
                          estimand = estimand,
                          method = method)
  
  data.frame(se_additive = stats::sd(boot_df$additive_effect),
             lower_ci_additive = stats::quantile(boot_df$additive_effect, 0.025),
             upper_ci_additive = stats::quantile(boot_df$additive_effect, 0.975),
             se_log_mult = stats::sd(boot_df$log_multiplicative_effect),
             lower_ci_mult = exp(stats::quantile(boot_df$log_multiplicative_effect, 0.025, na.rm = TRUE)), # added NA rm true for some bootstrap replicates
             upper_ci_mult = exp(stats::quantile(boot_df$log_multiplicative_effect, 0.975, na.rm = TRUE)),
             se_psi_1 = stats::sd(boot_df$psi_1),
             lower_ci_psi_1 = stats::quantile(boot_df$psi_1, 0.025),
             upper_ci_psi_1 = stats::quantile(boot_df$psi_1, 0.975),
             se_psi_0 = stats::sd(boot_df$psi_0),
             lower_ci_psi_0 = stats::quantile(boot_df$psi_0, 0.025),
             upper_ci_psi_0 = stats::quantile(boot_df$psi_0, 0.975))
}

#' Bootstrap inference for bounds
#'
#' Computes bootstrap-based standard errors and confidence intervals for lower
#' and upper bounds on additive and multiplicative effects.
#'
#' @param boot_estimates A list of bootstrap results.
#' @param estimand Character string specifying the estimand.
#' @param method Character string specifying the method (typically \code{"bound"}).
#'
#' @return A data.frame containing:
#' \describe{
#'   \item{number_NA_replicates}{Number of bootstrap replicates with missing bounds.}
#'   \item{se_additive_lower}{Standard error of lower bound (additive scale).}
#'   \item{se_additive_upper}{Standard error of upper bound (additive scale).}
#'   \item{lower_ci_additive_lower}{Lower CI for additive lower bound.}
#'   \item{upper_ci_additive_upper}{Upper CI for additive upper bound.}
#'   \item{se_log_mult_lower}{Standard error of log lower multiplicative bound.}
#'   \item{se_log_mult_upper}{Standard error of log upper multiplicative bound.}
#' }
#'
get_boot_se_bound <- function(boot_estimates, estimand = "nat_inf", method = "bound"){
  boot_df <- make_boot_df(boot_estimates = boot_estimates,
                          estimand = estimand,
                          method = method)
  
  number_NA <- length(which(is.na(boot_df$additive_effect_lower)))
  
  data.frame(number_NA_replicates = number_NA,
             se_additive_lower = stats::sd(boot_df$additive_effect_lower, na.rm = TRUE),
             lower_ci_additive_lower = stats::quantile(boot_df$additive_effect_lower, 0.025, na.rm = TRUE),
             upper_ci_additive_lower = stats::quantile(boot_df$additive_effect_lower, 0.975, na.rm = TRUE),
             se_additive_upper = stats::sd(boot_df$additive_effect_upper, na.rm = TRUE),
             lower_ci_additive_upper = stats::quantile(boot_df$additive_effect_upper, 0.025, na.rm = TRUE),
             upper_ci_additive_upper = stats::quantile(boot_df$additive_effect_upper, 0.975, na.rm = TRUE),
             #original
             #se_mult_lower = stats::sd(boot_df$mult_effect_lower),
             #new
             se_log_mult_lower = stats::sd(log(boot_df$mult_effect_lower), na.rm = TRUE),
             lower_ci_mult_lower = stats::quantile(boot_df$mult_effect_lower, 0.025, na.rm = TRUE),
             upper_ci_mult_lower = stats::quantile(boot_df$mult_effect_lower, 0.975, na.rm = TRUE),
             #original
             #se_mult_upper = stats::sd(boot_df$mult_effect_upper),
             #new
             se_log_mult_upper = stats::sd(log(boot_df$mult_effect_upper), na.rm = TRUE),
             lower_ci_mult_upper = stats::quantile(boot_df$mult_effect_upper, 0.025, na.rm = TRUE),
             upper_ci_mult_upper = stats::quantile(boot_df$mult_effect_upper, 0.975, na.rm = TRUE))
}

#' Bootstrap inference for sensitivity analysis
#'
#' Computes bootstrap standard errors and confidence intervals for each value
#' of the sensitivity parameter \eqn{\epsilon}.
#'
#' @param boot_estimates A list of bootstrap results.
#' @param estimand Character string specifying the estimand.
#' @param method Character string specifying the method (typically \code{"sens_cw"}).
#'
#' @return A data.frame with one row per \eqn{\epsilon}, containing:
#' \describe{
#'   \item{epsilon}{Sensitivity parameter value.}
#'   \item{se_additive}{Standard error of additive effect.}
#'   \item{lower_ci_additive}{Lower CI for additive effect.}
#'   \item{upper_ci_additive}{Upper CI for additive effect.}
#'   \item{se_mult}{Standard error of log multiplicative effect.}
#'   \item{lower_ci_mult}{Lower CI (multiplicative scale).}
#'   \item{upper_ci_mult}{Upper CI (multiplicative scale).}
#' }
#'
get_boot_se_sens_cw <- function(boot_estimates, estimand = "nat_inf", method = "sens_cw"){
  boot_df <- make_boot_df(boot_estimates = boot_estimates,
                          estimand = estimand,
                          method = method)
  
  epsilon <- unique(boot_df$epsilon)
  boot_res_list <- vector("list", length = length(epsilon))
  names(boot_res_list) <- paste0("epsilon_", epsilon)
  
  for(e in 1:length(epsilon)){
    # QUESTION handling NAs 
    boot_res_list[[e]] <- data.frame(epsilon = epsilon[e],
                                     se_additive = stats::sd(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], na.rm = TRUE),
                                     lower_ci_additive = stats::quantile(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], 0.025, na.rm = TRUE),
                                     upper_ci_additive = stats::quantile(boot_df$additive_effect[boot_df$epsilon == epsilon[e]], 0.975, na.rm = TRUE),
                                     se_mult = stats::sd(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], na.rm = TRUE),
                                     lower_ci_mult = exp(stats::quantile(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], 0.025, na.rm = TRUE)),
                                     upper_ci_mult = exp(stats::quantile(boot_df$log_multiplicative_effect[boot_df$epsilon == epsilon[e]], 0.975, na.rm = TRUE)))
  }
  
  return(do.call(rbind, boot_res_list))
  
}

#' Example dataset inspired by the PROVIDE study
#'
#' A simulated dataset based on PROVIDE study data for analysis of vaccine on antibiotic use
#'
#' @format A data frame with multiple observations and the following variables:
#' \describe{
#'   \item{wk10_haz}{Height-for-age Z-score at 10 weeks (numeric).}
#'   \item{gender}{Infant gender ("Male", "Female").}
#'   \item{num_hh_sleep}{Number of individuals sleeping in the household (integer).}
#'   \item{rotaarm}{Vaccination indicator (0 = control, 1 = vaccine).}
#'   \item{rotaepi}{Infection indicator (0 = no episode, 1 = episode).}
#'   \item{any_abx_wk52}{Any antibiotic use by 52 weeks (0/1).}
#' }
#'
#' @source Simulated data loosely based on the PROVIDE study.
"provide"

#' Print method for \code{"vaxstrat"} objects
#'
#' Displays formatted estimates and confidence intervals for growth effects.
#'
#' @param x An object of class \code{"vaxstrat"}.
#' @param scale Character string specifying effect scale:
#'   \code{"additive"} (default) or \code{"multiplicative"}.
#' @param ... Additional arguments (not used).
#'
#' @details
#' Results are grouped by estimand (e.g., naturally infected, doomed, population)
#' and estimation method. Confidence intervals are drawn either from bootstrap
#' summaries (if available) or normal approximations.
#'
#' @return Invisibly returns \code{x}.
#'
#' @method print vaxstrat
#' @export
print.vaxstrat <- function(x, scale = "additive", ...) {
  
  # Helper to print one row
  print_row <- function(label, method, est, lower, upper) {
    cat(sprintf("%-25s%-15s%-15.4f%-15.4f%-15.4f\n", label, method, est, lower, upper))
  }
  
  # Helper to get estimates
  extract_estimates <- function(est_obj, method, scale, label_prefix) {
    if (scale == "additive") {
      if (is.null(est_obj$boot_se)) {
        est <- est_obj$pt_est['additive_effect']
        se <- est_obj$pt_est['additive_se']
        print_row(label_prefix, method, est, est - 1.96 * se, est + 1.96 * se)
      } else {
        est <- est_obj$pt_est['additive_effect']
        lower <- est_obj$boot_se$lower_ci_additive
        upper <- est_obj$boot_se$upper_ci_additive
        print_row(label_prefix, method, est, lower, upper)
      }
    } else {
      if (is.null(est_obj$boot_se)) {
        est <- exp(est_obj$pt_est['log_multiplicative_effect'])
        se <- est_obj$pt_est['log_multiplicative_se']
        print_row(label_prefix, method, est,
                  exp(log(est) - 1.96 * se),
                  exp(log(est) + 1.96 * se))
      } else {
        est <- exp(est_obj$pt_est['log_multiplicative_effect'])
        lower <- est_obj$boot_se$lower_ci_mult
        upper <- est_obj$boot_se$upper_ci_mult
        print_row(label_prefix, method, est, lower, upper)
      }
    }
  }
  
  # Helper to print bounds
  print_bounds <- function(est_obj, scale, label_prefix) {
    if (scale == "additive") {
      print_row(label_prefix, "Lower Bound",
                est_obj$pt_est['additive_effect_lower'],
                est_obj$boot_se$lower_ci_additive_lower,
                est_obj$boot_se$upper_ci_additive_lower)
      print_row(label_prefix, "Upper Bound",
                est_obj$pt_est['additive_effect_upper'],
                est_obj$boot_se$lower_ci_additive_upper,
                est_obj$boot_se$upper_ci_additive_upper)
    } else {
      print_row(label_prefix, "Lower Bound",
                est_obj$pt_est['mult_effect_lower'],
                est_obj$boot_se$lower_ci_mult_lower,
                est_obj$boot_se$upper_ci_mult_lower)
      print_row(label_prefix, "Upper Bound",
                est_obj$pt_est['mult_effect_upper'],
                est_obj$boot_se$lower_ci_mult_upper,
                est_obj$boot_se$upper_ci_mult_upper)
    }
  }
  
  # Header
  scale_label <- ifelse(scale == "additive", "Additive", "Multiplicative")
  cat(sprintf("%50s\n", paste("                         Growth Effect Estimation Results:", scale_label)))
  cat(paste(rep("-", 90), collapse = ""), "\n")
  col_names <- c("Estimand", "Method", "Point Est.", "95% CI: Lower", "95% CI: Upper")
  cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", col_names[1], col_names[2], col_names[3], col_names[4], col_names[5]))
  cat(paste(rep("-", 90), collapse = ""), "\n")
  
  # Loop through estimands
  lapply(x, function(i) {
    if (inherits(i, "nat_inf")) {
      
      if(any(c("gcomp_CW", "ipw", "aipw", "bound", "tmle", "aipw_CW") %in% names(i))){
        cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", "Naturally Infected - - - ", "- - - - - - - -", " - - - - - - - ", "- - - - - - - -", " - - - - - - - - - - "))
        
        lapply(i, function(j) {
          if (inherits(j, "gcomp_CW")) extract_estimates(j, "G-Computation (Cross-World only)", scale, "")
          if (inherits(j, "ipw")) extract_estimates(j, "IPW", scale, "")
          if (inherits(j, "aipw_CW")) extract_estimates(j, "AIPW (Cross-World only)", scale, "")
          if (inherits(j, "tmle")) extract_estimates(j, "TMLE", scale, "")
          if (inherits(j, "bound")) print_bounds(j, scale, "")
        })
        
      }
      
      if(any(c("gcomp_ER", "gcomp_ER_CW", "ipw_ER", "aipw_ER", "aipw_ER_CW") %in% names(i))){
        cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", "Naturally Infected (ER) - ", "- - - - - - - -", " - - - - - - - ", "- - - - - - - -", " - - - - - - - - - - "))
        lapply(i, function(j) {
          if (inherits(j, "gcomp_ER")) extract_estimates(j, "G-Computation (ER only)", scale, "")
          if (inherits(j, "gcomp_ER_CW")) extract_estimates(j, "G-Computation (ER + Cross-World)", scale, "")
          if (inherits(j, "ipw_ER")) extract_estimates(j, "IPW", scale, "")
          if (inherits(j, "aipw_ER")) extract_estimates(j, "AIPW (ER only)", scale, "")
          if (inherits(j, "aipw_ER_CW")) extract_estimates(j, "AIPW (ER + Cross-World)", scale, "")
        })
      }
      
    }
    
    if (inherits(i, "doomed")) {
      cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", "Doomed - - - - - - - - -", "- - - - - - - -", " - - - - - - - ", "- - - - - - - -", " - - - - - - - - - - "))
      lapply(i, function(j) {
        if (inherits(j, "gcomp")) extract_estimates(j, "G-Computation", scale, "")
        if (inherits(j, "ipw")) extract_estimates(j, "IPW", scale, "")
        if (inherits(j, "aipw")) extract_estimates(j, "AIPW", scale, "")
        if (inherits(j, "tmle")) extract_estimates(j, "TMLE", scale, "")
        if (inherits(j, "bound")) print_bounds(j, scale, "")
      })
    }
    
    if (inherits(i, "pop")) {
      cat(sprintf("%-25s%-15s%-15s%-15s%-15s\n", "Population - - - - - - - ", "- - - - - - - -", " - - - - - - - ", "- - - - - - - -", " - - - - - - - - - - "))
      lapply(i, function(j) {
        if (inherits(j, "gcomp")) extract_estimates(j, "G-Computation", scale, "")
        if (inherits(j, "ipw")) extract_estimates(j, "IPW", scale, "")
        if (inherits(j, "aipw")) extract_estimates(j, "AIPW", scale, "")
        if (inherits(j, "tmle")) extract_estimates(j, "TMLE", scale, "")
      })
    }
  })
  
  invisible(x)
}


#' Plot sensitivity analysis results
#'
#' Produces a line plot of estimated effects across values of the sensitivity
#' parameter \eqn{\epsilon}.
#'
#' @param object An object of class \code{"sens_cw"}.
#' @param se Logical; if TRUE, includes 95\% confidence bands.
#' @param effect_type Character string specifying effect scale:
#'   \code{"additive"} or \code{"multiplicative"}.
#' @param ... Additional arguments (not used).
#'
#' @details
#' Confidence intervals are constructed using normal approximations based on
#' bootstrap standard errors.
#'
#' @return Produces a \code{ggplot2} plot and returns it invisibly.
#' @export
#' 
#' @method plot sens_cw
plot.sens_cw <- function(
  object, se = TRUE, effect_type = c("additive", "multiplicative"), 
  ...
) {
  effect_type <- match.arg(effect_type)

  if (!inherits(object, "sens_cw")) {
    stop("Object must be of class 'sens_cw'")
  }

  df <- object
  class(df) <- "data.frame"

  if (effect_type == "additive") {
    df$effect <- df$additive_effect
    if (se) {
      df$lower <- df$additive_effect - 1.96 * df$additive_se
      df$upper <- df$additive_effect + 1.96 * df$additive_se
    }
    ylab <- "Additive Effect"
  } else {
    df$effect <- exp(df$log_multiplicative_effect)
    if (se) {
      df$lower <- exp(df$log_multiplicative_effect - 1.96 * df$log_multiplicative_se)
      df$upper <- exp(df$log_multiplicative_effect + 1.96 * df$log_multiplicative_se)
    }
    ylab <- "Multiplicative Effect"
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = epsilon, y = effect)) +
    ggplot2::geom_line() +
    ggplot2::labs(x = expression(epsilon), y = ylab) +
    ggplot2::theme_minimal()

  if (se) {
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2)
  }

  print(p)
}
