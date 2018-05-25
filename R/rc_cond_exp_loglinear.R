#' Regression Calibration with Loglinear Measurement Error Model
#'
#' Implements the regression calibration method described by Lyles and Kupper
#' (\emph{J. Agric. Biol. Environ. Stat.} 2012).
#'
#' The true disease model is a GLM:
#'
#' g[E(Y)] = beta_0 + beta_z Z + \strong{beta_c}^T \strong{C} +
#' \strong{beta_b}^T \strong{B}
#'
#' The measurement error model is:
#'
#' log(Z) = alpha_0 + \strong{alpha_d}^T \strong{D} + \strong{alpha_c}^T
#' \strong{C} + d, d ~ N(0, sigsq_d)
#'
#' The procedure is as follows: in the validation study, fit the measurement
#' error model to estimate (\strong{alpha}, sigsq_d); in the main study,
#' calculate E(Z|\strong{D},\strong{C}) = exp(alpha_0 + \strong{alpha_d}^T
#' \strong{D} + \strong{alpha_c}^T \strong{C} + sigq_d / 2 and fit the true
#' disease model with those values in place of the unobserved Z's.
#'
#' @inheritParams rc_cond_exp
#'
#'
#' @inherit rc_cond_exp return
#'
#'
#' @references
#' Lyles, R.H. and Kupper, L.L. (2012) "Approximate and pseudo-likelihood
#' analysis for logistic regression using external validation data to model log
#' exposure." \emph{J. Agric. Biol. Environ. Stat.} \strong{18}(1): 22-38.
#'
#'
#' @export
# # Data for testing
# n.m <- 50000
# n.e <- 50000
# n <- n.m + n.e
#
# alphas <- c(0, 0.25, 0.25)
# sigsq_d <- 0.1
#
# betas <- c(0, 0.25, 0.1)
# sigsq_e <- 0.5
#
# d <- rnorm(n)
# c <- rnorm(n)
# z <- exp(alphas[1] + alphas[2] * d + alphas[3] * c + rnorm(n, sd = sqrt(sigsq_d)))
# y <- betas[1] + betas[2] * z + betas[3] * c + rnorm(n, sd = sqrt(sigsq_e))
#
# all_data <- data.frame(y = y, z = z, c = c, d = d)
# all_data[1: n.e, 1] <- NA
# all_data[(n.e + 1): (n.e + n.m), 2] <- NA
# main <- internal <- external <- NULL
# y_var <- "y"
# z_var <- "z"
# d_vars <- "d"
# c_vars <- "c"
# b_vars <- NULL
# tdm_covariates <- mem_covariates <- NULL
# tdm_family <- "gaussian"
# all_imputed <- FALSE
# boot_var <- TRUE
# boots <- 100
# fit <- rc_cond_exp_loglinear(all_data = all_data,
#                              y_var = "y",
#                              z_var = "z",
#                              d_var = "d",
#                              c_vars = "c")

rc_cond_exp_loglinear <- function(all_data = NULL,
                                  main = NULL,
                                  internal = NULL,
                                  external = NULL,
                                  y_var,
                                  z_var,
                                  d_vars = NULL,
                                  c_vars = NULL,
                                  b_vars = NULL,
                                  tdm_covariates = NULL,
                                  tdm_family = "gaussian",
                                  mem_covariates = NULL,
                                  all_imputed = FALSE,
                                  boot_var = FALSE,
                                  boots = 100) {

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
  all_data_o <- all_data

  # Fit MEM using all available data
  mem.formula <- paste(paste("log(", z_var, ")", " ~ ", sep = ""),
                       paste(mem_covariates, collapse = " + ", sep = ""))
  mem.fit <- lm(mem.formula, data = all_data)
  sigsq_delta.hat <- rev(anova(mem.fit)$"Mean Sq")[1]

  # Replace all or just missing Z's with Zhat's, depending on all_imputed input
  if (all_imputed) {

    mu_logz <- predict(mem.fit, newdata = all_data, type = "response")
    all_data[, z_var] <- exp(mu_logz + sigsq_delta.hat / 2)

  } else {

    locs.missingz <- which(is.na(all_data[, z_var]))
    mu_logz <- predict(mem.fit, newdata = all_data[locs.missingz, ], type = "response")
    all_data[locs.missingz, z_var] <- exp(mu_logz + sigsq_delta.hat / 2)

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
  theta.labels <- c(beta.labels, alpha.labels, "sigsq_delta")

  # Extract point estimates for theta
  theta.hat <- c(tdm.fit$coef, mem.fit$coef, sigsq_delta.hat)
  names(theta.hat) <- theta.labels

  # Add theta.hat to ret.list
  ret.list <- list(theta.hat = theta.hat)

  # Get bootstrap variance estimate if requested
  if (boot_var) {

    # Various data types
    locs.m <- which(complete.cases(all_data_o[, c(y_var, mem_covariates)]) &
                      is.na(all_data_o[, z_var]))
    locs.i <- which(complete.cases(all_data_o[, c(y_var, tdm_covariates, mem_covariates)]))
    locs.e <- which(complete.cases(all_data_o[, c(z_var, mem_covariates)]) &
                      is.na(all_data_o[, y_var]))

    # Initialize matrix for theta estimates
    theta.hat.boots <- matrix(NA, ncol = length(theta.hat), nrow = boots)

    # Bootstrap
    for (ii in 1: boots) {

      theta.hat.boots[ii, ] <- rc_cond_exp_loglinear(
        all_data = all_data_o[c(sample(locs.m, replace = TRUE),
                                sample(locs.i, replace = TRUE),
                                sample(locs.e, replace = TRUE)), ],
        y_var = y_var,
        z_var = z_var,
        tdm_covariates = tdm_covariates,
        tdm_family = tdm_family,
        mem_covariates = mem_covariates,
        all_imputed = all_imputed)

    }

    # Calculate bootstrap variance estimates
    boot.variance <- var(theta.hat.boots)
    rownames(boot.variance) <- colnames(boot.variance) <- theta.labels
    ret.list$boot.var <- boot.variance

    boot.ci <- apply(theta.hat.boots, 2, function(x) quantile(x, probs = c(0.025, 0.975)))
    colnames(boot.ci) <- theta.labels
    ret.list$boot.ci <- boot.ci

  }

  # Return ret.list
  if (length(ret.list) == 1) {
    ret.list <- ret.list[[1]]
  }
  return(ret.list)

}
