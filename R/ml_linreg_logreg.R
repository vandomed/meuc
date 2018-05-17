#' Maximum Likelihood with Linear Regression Disease Model and Logistic
#' Regression Measurement Error Model
#'
#' Calculates maximum likelihood estimates for measurement error/unmeasured
#' confounding scenario where true disease model is linear regression and
#' measurement error model is logistic regression.
#'
#' The true disease model is:
#'
#' Y = beta_0 + beta_z Z + \strong{beta_c}^T \strong{C} + \strong{beta_b}^T
#' \strong{B} + epsilon, epsilon ~ N(0, sigsq_epsilon)
#'
#' The measurement error model is:
#'
#' logit[P(Z = 1)] = alpha_0 + \strong{alpha_d}^T \strong{D} +
#' \strong{alpha_c}^T \strong{C}
#'
#' There should be main study data with (Y, \strong{D}, \strong{C}, \strong{B})
#' as well as internal validation data with
#' (Y, Z, \strong{D}, \strong{C}, \strong{B}) and/or external validation data
#' with (Z, \strong{D}, \strong{C}). Parameters are theoretically identifiable
#' without validation data, but estimation may be impractical in that scenario.
#'
#'
#' @inheritParams ml_linreg_linreg
#'
#'
#' @inherit ml_linreg_linreg return
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
#
# betas <- c(0, 0.25, 0.1, 0.05)
# sigsq_e <- 0.5
#
# d <- rnorm(n)
# c <- rnorm(n)
# b <- rnorm(n)
# z <- rbinom(n, size = 1,
#             prob = (1 + exp(-alphas[1] - alphas[2] * d - alphas[3] * c))^(-1))
# y <- betas[1] + betas[2] * z + betas[3] * c + betas[4] * b + rnorm(n, sd = sqrt(sigsq_e))
#
# all_data <- data.frame(y = y, z = z, c = c, b = b, d = d)
# all_data[1: n.e, 1] <- NA
# all_data[(n.e + 1): (n.e + n.m), 2] <- NA
# main <- internal <- external <- NULL
# y_var <- "y"
# z_var <- "z"
# d_vars <- "d"
# c_vars <- "c"
# b_vars <- "b"
# tdm_covariates <- mem_covariates <- NULL
# estimate_var <- TRUE
#
# fit <- ml_linreg_logreg(all_data = all_data,
#                         y_var = "y",
#                         z_var = "z",
#                         d_vars = "d",
#                         c_vars = "c",
#                         b_vars = "b")

ml_linreg_logreg <- function(all_data = NULL,
                             main = NULL,
                             internal = NULL,
                             external = NULL,
                             y_var,
                             z_var,
                             d_vars = NULL,
                             c_vars = NULL,
                             b_vars = NULL,
                             tdm_covariates = NULL,
                             mem_covariates = NULL,
                             estimate_var = TRUE,
                             ...) {

  # If tdm.covariates and mem.covariates specified, figure out d.vars, c.vars,
  # and b.vars
  if (! is.null(tdm_covariates) & ! is.null(mem_covariates)) {
    tdm_covariates <- setdiff(tdm_covariates, z_var)
    d_vars <- setdiff(mem_covariates, tdm_covariates)
    c_vars <- intersect(tdm_covariates, mem_covariates)
    b_vars <- setdiff(tdm_covariates, mem_covariates)
  }

  # Get dimension of D, C, and B
  kd <- length(d_vars)
  kc <- length(c_vars)
  kb <- length(b_vars)

  # Get covariate lists
  zdcb <- c(z_var, d_vars, c_vars, b_vars)
  dcb <- c(d_vars, c_vars, b_vars)
  dc <- c(d_vars, c_vars)

  # Subset various data types
  if (! is.null(all_data)) {
    main <- all_data[complete.cases(all_data[, c(y_var, dcb)]) & is.na(all_data[, z_var]), ]
    internal <- all_data[complete.cases(all_data[, c(y_var, zdcb)]), ]
    external <- all_data[is.na(all_data[, y_var]) & complete.cases(all_data[, dc]), ]
  }

  n.m <- nrow(main)
  some.m <- n.m > 0
  if (some.m) {
    y.m <- main[, y_var]
    onecb.m <- as.matrix(cbind(rep(1, n.m), main[, c(c_vars, b_vars)]))
    onedc.m <- as.matrix(cbind(rep(1, n.m), main[, c(d_vars, c_vars)]))
  }

  n.i <- nrow(internal)
  some.i <- n.i > 0
  if (some.i) {
    y.i <- internal[, y_var]
    z.i <- internal[, z_var]
    onezcb.i <- as.matrix(cbind(rep(1, n.i), internal[, c(z_var, c_vars, b_vars)]))
    onedc.i <- as.matrix(cbind(rep(1, n.i), internal[, c(d_vars, c_vars)]))
  }

  n.e <- nrow(external)
  some.e <- n.e > 0
  if (some.e) {
    z.e <- external[, z_var]
    onedc.e <- as.matrix(cbind(rep(1, n.i), external[, c(d_vars, c_vars)]))
  }

  # Get number of betas and alphas
  n.betas <- 2 + kc + kb
  n.alphas <- 1 + kd + kc

  # Get indices for parameters being estimated and create labels
  loc.betas <- 1: n.betas
  beta.labels <- paste("beta", c("0", z_var, c_vars, b_vars), sep = "_")

  loc.alphas <- (n.betas + 1): (n.betas + n.alphas)
  alpha.labels <- paste("alpha", c("0", d_vars, c_vars), sep = "_")

  loc.sigsq_e <- n.betas + n.alphas + 1

  theta.labels <- c(beta.labels, alpha.labels, "sigsq_e")

  # Log-likelihood function
  llf <- function(f.theta, estimating.hessian = FALSE) {

    # Extract parameters
    f.betas <- matrix(f.theta[loc.betas], ncol = 1)
    f.beta_z <- f.betas[2]

    f.alphas <- matrix(f.theta[loc.alphas], ncol = 1)

    f.sigsq_e <- f.theta[loc.sigsq_e]

    if (some.m) {

      # Likelihood for main study subjects:
      # L = \sum_z f(Y|Z,C,B) P(Z=z|D,C)
      mu_y.0cb <- onecb.m %*% f.betas[-2]
      mu_y.1cb <- mu_y.0cb + f.beta_z
      p_z.dc <- (1 + exp(-onedc.m %*% f.alphas))^(-1)
      ll.m <- sum(log(
        dnorm(y.m, mean = mu_y.0cb, sd = sqrt(f.sigsq_e)) *
          dbinom(0, size = 1, prob = p_z.dc) +
          dnorm(y.m, mean = mu_y.1cb, sd = sqrt(f.sigsq_e)) *
          dbinom(1, size = 1, prob = p_z.dc)
      ))

    } else {
      ll.m <- 0
    }

    if (some.i) {

      # Likelihood for internal validation subjects:
      # L = f(Y|Z,C,B) f(Z|D,C)
      ll.i <- sum(
        dnorm(y.i, log = TRUE,
              mean = onezcb.i %*% f.betas,
              sd = sqrt(f.sigsq_e)) +
          dbinom(z.i, log = TRUE, size = 1,
                 prob = (1 + exp(-onedc.i %*% f.alphas))^(-1))
      )

    } else {
      ll.i <- 0
    }

    if (some.e) {

      # Likelihood for external validation subjects:
      # L = f(Z|D,C)
      ll.e <- sum(
        dbinom(z.e, log = TRUE, size = 1,
               prob = (1 + exp(-onedc.e %*% f.alphas))^(-1))
      )

    } else {
      ll.e <- 0
    }

    # Return negative log-likelihood
    ll <- ll.m + ll.i + ll.e
    return(-ll)

  }

  # Create list of extra arguments, and assign default starting values and
  # lower values if not specified by user
  extra.args <- list(...)
  if (is.null(extra.args$start)) {
    extra.args$start <- c(rep(0.01, n.betas + n.alphas), 1)
  }
  if (is.null(extra.args$lower)) {
    extra.args$lower <- c(rep(-Inf, n.betas + n.alphas), 1e-4)
  }
  if (is.null(extra.args$control$rel.tol)) {
    extra.args$control$rel.tol <- 1e-6
  }
  if (is.null(extra.args$control$eval.max)) {
    extra.args$control$eval.max <- 1000
  }
  if (is.null(extra.args$control$iter.max)) {
    extra.args$control$iter.max <- 750
  }

  # Obtain ML estimates
  ml.max <- do.call(nlminb, c(list(objective = llf), extra.args))

  # Create list to return
  theta.hat <- ml.max$par
  names(theta.hat) <- theta.labels
  ret.list <- list(theta.hat = theta.hat)

  # If requested, add variance-covariance matrix to ret.list
  if (estimate_var) {
    hessian.mat <- pracma::hessian(f = llf, estimating.hessian = TRUE,
                                   x0 = theta.hat)
    theta.variance <- try(solve(hessian.mat), silent = TRUE)
    if (class(theta.variance) == "try-error") {
      message("Estimated Hessian matrix is singular, so variance-covariance matrix cannot be obtained.")
      ret.list$theta.var <- NULL
    } else {
      colnames(theta.variance) <- rownames(theta.variance) <- theta.labels
      ret.list$theta.var <- theta.variance
    }
  }

  # Add nlminb object and AIC to ret.list
  ret.list$nlminb.object <- ml.max
  ret.list$aic <- 2 * (length(theta.hat) + ml.max$objective)

  # Return ret.list
  return(ret.list)

}
