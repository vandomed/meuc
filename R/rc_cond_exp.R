#' Regression Calibration (Conditional Expectation Method)
#'
#' Implements the "conditional expectation" version of regression calibration as
#' described in Rosner et al. (\emph{Stat. Med.} 1989). For the "algebraic"
#' version, see \code{\link{rc_algebraic}}.
#'
#' The true disease model is a GLM:
#'
#' g[E(Y)] = beta_0 + beta_z Z + \strong{beta_c}^T \strong{C} +
#' \strong{beta_b}^T \strong{B}
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
#' @param d_vars Character string specifying name of D variables.
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
#' @references
#' Kuha, J. (1994) "Corrections for exposure measurement error in logistic
#' regression models with an application to nutritional data." \emph{Stat. Med.}
#' \strong{13}(11): 1135-1148.
#'
#' Lyles, R.H. and Kupper, L.L. (2012) "Approximate and pseudo-likelihood
#' analysis for logistic regression using external validation data to model log
#' exposure." \emph{J. Agric. Biol. Environ. Stat.} \strong{18}(1): 22-38.
#'
#' Rosner, B., Willett, W. and Spiegelman, D. (1989) "Correction of logistic
#' regression relative risk estimates and confidence intervals for systematic
#' within-person measurement error." \emph{Stat. Med.} \strong{8}(9): 1051-69.
#'
#' Spiegelman, D., Carroll, R.J. and Kipnis, V. (2001) "Efficient regression
#' calibration for logistic regression in main study/internal validation study
#' designs with an imperfect reference instrument." \emph{Stat. Med.}
#' \strong{20}(1): 139-160.
#'
#'
#' @export
# # Data for testing
# n.m <- 1000
# n.e <- 1000
# n.i <- 1000
# n <- n.m + n.e + n.i
#
# alphas <- c(0, 0.25, 0.25)
# sigsq_d <- 0.5
#
# betas <- c(0, 0.25, 0.1)
# sigsq_e <- 0.5
#
# d <- rnorm(n)
# c <- rnorm(n)
# z <- alphas[1] + alphas[2] * d + alphas[3] * c + rnorm(n, sd = sqrt(sigsq_d))
# y <- betas[1] + betas[2] * z + betas[3] * c + rnorm(n, sd = sqrt(sigsq_e))
#
# all_data <- data.frame(y = y, z = z, c = c, d = d)
# all_data[1: n.e, 1] <- NA
# all_data[(n.e + 1): (n.e + n.m), 2] <- NA
# main <- internal <- external <- NULL
# main <- internal <- external <- NULL
# y_var <- "y"
# z_var <- "z"
# d_vars <- "d"
# c_vars <- "c"
# b_vars <- NULL
# tdm_covariates <- mem_covariates <- NULL
# tdm_family <- "gaussian"
# mem_family <- "gaussian"
# beta_0_formula <- 1
# all_imputed <- FALSE
# boot_var <- TRUE
# boots <- 100
# fit1 <- rc_algebraic(all_data = all_data,
#                      y_var = "y",
#                      z_var = "z",
#                      d_var = "d",
#                      c_vars = "c")
# fit2 <- rc_cond_exp(all_data = all_data,
#                     y_var = "y",
#                     z_var = "z",
#                     d_var = "d",
#                     c_vars = "c")
# fit1
# fit2

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
    locs.m <- which(! is.na(all_data[, y_var]) & is.na(all_data[, z_var]))
    locs.i <- which(! is.na(all_data[, y_var]) & ! is.na(all_data[, z_var]))
    locs.e <- which(is.na(all_data[, y_var]) & ! is.na(all_data[, z_var]))

    # Initialize matrix for theta estimates
    theta.hat.boots <- matrix(NA, ncol = length(theta.hat), nrow = boots)

    # Bootstrap
    for (ii in 1: boots) {

      theta.hat.boots[ii, ] <- rc_cond_exp(
        all_data = all_data[c(sample(locs.m, replace = TRUE),
                              sample(locs.i, replace = TRUE),
                              sample(locs.e, replace = TRUE)), ],
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
