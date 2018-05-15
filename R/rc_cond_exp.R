#' Regression Calibration (Conditional Expectation Method)
#'
#' Implements the "conditional expectation" version of regression calibration as
#' described in Rosner et al. (\emph{Stat. Med.} 1989). For the "algebraic"
#' version, see \code{\link{rc_algebraic}}.
#'
#' The true disease model is a GLM:
#'
#' g[E(Y)] = beta_0 + beta_z Z + beta_c^T C + beta_b^T B
#'
#' The measurement error model is:
#'
#' h[E(Z)] = \strong{alpha_d}^T \strong{D} + \strong{alpha_c}^T \strong{C}
#'
#' The procedure is as follows: in the validation study, fit the measurement
#' error model to estimate alpha's; in the main study, calculate E(Z|D,C) and
#' fit the true disease model with those values in place of the unobserved Z's.
#'
#' @param all_data Data frame with data for main study and validation study.
#' @param main Data frame with data for the main study.
#' @param internal Data frame with data for internal validation study.
#' @param external Data frame with data for the external validation study.
#' @param y_var Character string specifying name of Y variable.
#' @param z_var Character string specifying name of Z variable.
#' @param d_var Character string specifying name of D variables.
#' @param c_vars Character vector specifying names of C variables.
#' @param b_vars Character vector specifying names of variables in true disease
#' model but not in measurement error model.
#' @param tdm_covariates Character vector specifying variables in true disease
#' model. The Z variable is automatically included whether you include it in
#' \code{tdm_covariates} or not.
#' @param tdm_family Character string specifying family of true disease model
#' (see \code{\link[stats]{glm}}).
#' @param mem_covariates Character vector specifying variables in measurement
#' error model.
#' @param mem_family Character string specifying family of measurement error
#' model (see \code{\link[stats]{glm}}).
#' @param all_imputed Logical value for whether to use imputed Z's for all
#' subjects, even those with the actual Z's observed (i.e. internal validation
#' subjects).
#' @param boot_var Logical value for whether to calculate a bootstrap
#' variance-covariance matrix.
#' @param boots Numeric value specifying number of bootstrap samples to use.
#'
#'
#' @return
#' If \code{boot_var = TRUE}, list containing parameter estimates and
#' variance-covariance matrix; otherwise just the parameter estimates.
#'
#'
#' @inherit rc_algebraic references
#'
#'
#' @export
rc_cond_exp <- function(all_data = NULL,
                        main = NULL,
                        internal = NULL,
                        external = NULL,
                        y_var,
                        z_var,
                        d_vars = NULL,
                        c_vars = NULL,
                        b_vars = NULL,
                        tdm_covariates = NULL, tdm_family = "gaussian",
                        mem_covariates = NULL, mem_family = "gaussian",
                        all_imputed = FALSE,
                        boot_var = FALSE, boots = 100) {

  # If tdm_covariates is NULL, construct it as (Z, C, B)
  if (is.null(tdm_covariates)) {
    tdm_covariates <- c(z_var, c_vars, b_vars)
  }

  # If mem_covariates is NULL, construct it as (D, C)
  if (is.null(mem_covariates)) {
    mem_covariates <- c(d_vars, c_vars)
  }

  # Make sure z_var is included in tdm_covariates
  if (! z_var %in% tdm_covariates) {
    tdm_covariates <- c(z_var, tdm_covariates)
  }

  # Get full list of covariates
  covariates <- unique(c(tdm_covariates, mem_covariates))

  # If all_data not specified, create it from main, internal, and external
  if (is.null(all_data)) {

    if (! is.null(main)) {
      if (! z_var %in% names(main)) {
        main[, z_var] <- NA
      }
      main <- main[, c(y_var, covariates)]
      all_data <- main
    }

    if (! is.null(internal)) {
      internal <- internal[, c(y_var, covariates)]
      all_data <- rbind(all_data, internal)
    }

    if (! is.null(external)) {
      if (! y_var %in% names(external)) {
        external[, y_var] <- NA
      }
      external <- external[, c(y_var, covariates)]
      all_data <- rbind(all_data, external)
    }

  }

  # Get sample size
  n <- nrow(all_data)

  # Store original all_data data frame
  all_data_original <- all_data

  # Fit MEM using all available data
  mem.formula <- paste(paste(z_var, " ~ ", sep = ""),
                       paste(mem_covariates, collapse = " + "), sep = "")
  if (mem_family == "gaussian") {
    mem.fit <- lm(mem.formula, data = all_data)
    sigsq_delta.hat <- rev(anova(mem.fit)$"Mean Sq")[1]
  } else {
    mem.fit <- glm(mem.formula, data = all_data, family = mem_family)
  }

  # Replace all or just missing Z's with Zhat's, depending on all_imputed input
  if (all_imputed) {

    all_data[, z_var] <- predict(mem.fit, newdata = all_data, type = "response")

  } else {

    locs.missingz <- which(is.na(all_data[, z_var]))
    all_data[locs.missingz, z_var] <-
      predict(mem.fit, newdata = all_data[locs.missingz, ], type = "response")

  }

  # Fit TDM
  tdm.formula <- paste(paste(y_var, " ~ ", sep = ""),
                       paste(tdm_covariates, collapse = " + "),
                       sep = "")
  tdm.fit <- glm(tdm.formula, data = all_data, family = tdm_family)

  # Create labels for parameter estimates
  beta.labels <- c("beta_0", paste(rep("beta_", length(tdm_covariates)),
                                   tdm_covariates,
                                   sep = ""))
  alpha.labels <- c("alpha_0", paste(rep("alpha_", length(mem_covariates)),
                                     mem_covariates,
                                     sep = ""))
  theta.labels <- c(beta.labels, alpha.labels)

  # Extract point estimates for theta
  if (mem_family == "gaussian") {

    # theta = (beta^T, alpha^T, sigsq_delta)^T
    theta.hat <- c(tdm.fit$coef, mem.fit$coef, sigsq_delta.hat)
    theta.labels <- c(theta.labels, "sigsq_delta")
    names(theta.hat) <- theta.labels

  } else {

    # theta = (beta^T, alpha^T)^T
    theta.hat <- c(tdm.fit$coef, mem.fit$coef)
    names(theta.hat) <- theta.labels

  }

  # Add theta.hat to ret.list
  ret.list <- list(theta.hat = theta.hat)

  # Get bootstrap variance estimate if requested
  if (boot_var) {

    # Various data types
    locs.main <-
      which(is.na(all_data_original[, z_var]) &
              complete.cases(all_data_original[, setdiff(covariates, z_var)]))
    locs.internal <-
      which(complete.cases(all_data_original[, c(y_var, covariates)]))
    locs.external <-
      which(is.na(all_data_original[, y_var]) &
              complete.cases(all_data_original[, c(z_var, mem_covariates)]))
    n.main <- length(locs.main)
    n.internal <- length(locs.internal)
    n.external <- length(locs.external)

    # Initialize matrix for theta estimates
    theta.hat.boots <- matrix(NA, ncol = length(theta.hat), nrow = boots)

    # Bootstrap
    for (ii in 1: boots) {

      locs.main.sampled <- sample(locs.main, replace = TRUE)
      locs.internal.sampled <- sample(locs.internal, replace = TRUE)
      locs.external.sampled <- sample(locs.external, replace = TRUE)

      all_data_boot <- all_data_original[c(locs.main.sampled,
                                           locs.internal.sampled,
                                           locs.external.sampled), ]

      theta.hat.boots[ii, ] <- rc_cond_exp(all_data = all_data_boot,
                                           y_var = y_var,
                                           z_var = z_var,
                                           tdm_covariates = tdm_covariates,
                                           tdm_family = tdm_family,
                                           mem_covariates = mem_covariates,
                                           mem_family = mem_family,
                                           all_imputed = all_imputed)

    }

    # Calculate bootstrap variance estimate and add it to ret.list
    boot.variance <- var(theta.hat.boots)
    rownames(boot.variance) <- colnames(boot.variance) <- theta.labels
    ret.list$boot.var <- boot.variance

  }

  # Return ret.list
  if (length(ret.list) == 1) {
    ret.list <- ret.list[[1]]
  }
  return(ret.list)

}
