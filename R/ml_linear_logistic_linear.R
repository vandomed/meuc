#' Maximum Likelihood with Three Models: Linear Regression, Logistic
#' Regression, and Linear Regression
#'
#' Calculates maximum likelihood estimates for measurement error/unmeasured
#' confounding scenario where Y|(Z,X,\strong{C},\strong{B}) is linear regression,
#' X|(Z,\strong{C},\strong{B}) is logistic regression, and
#' Z|(\strong{D},\strong{C}) is linear regression.
#'
#' The true disease model is:
#'
#' Y = beta_0 + beta_x X + beta_z Z + \strong{beta_c}^T \strong{C} +
#' \strong{beta_b}^T \strong{B} + e, e ~ N(0, sigsq_e)
#'
#' The X|(Z,\strong{C},\strong{B}) model is:
#'
#' logit[P(X = 1)] = gamma_0 + gamma_z Z + \strong{gamma_c}^T \strong{C} +
#' \strong{gamma_b}^T \strong{B}
#'
#' The Z|(\strong{D},\strong{C}) model is:
#'
#' Z = alpha_0 + \strong{alpha_d}^T \strong{D} + \strong{alpha_c}^T \strong{C} +
#' d, d ~ N(0, sigsq_d)
#'
#' There should be main study data with
#' (Y, X, \strong{D}, \strong{C}, \strong{B}) as well as internal validation
#' data with (Y, X, Z, \strong{D}, \strong{C}, \strong{B}) and/or external
#' validation data with (Z, X, \strong{D}, \strong{C}). Parameters are
#' theoretically identifiable without validation data, but estimation may be
#' impractical in that scenario.
#'
#'
#' @inheritParams ml_logistic_linear
#' @param x_var Character string specifying name of X variable.
#'
#'
#' @inherit ml_linear_linear return
#'
#'
#' @export

# # Data for testing
# n.m <- 0
# n.e <- 50000
# n.i <- 100
# n <- n.m + n.e + n.i
#
# alphas <- c(0, 0.25, 0.25)
# sigsq_d <- 0.5
#
# gammas <- c(-0.2, 0.2, 0.3, 0.1)
#
# betas <- c(0, 0.25, 0.1, 0.05, 0.15)
# sigsq_e <- 0.5
#
# d <- rnorm(n)
# c <- rnorm(n)
# b <- rnorm(n)
#
# z <- alphas[1] + alphas[2] * d + alphas[3] * c + rnorm(n, sd = sqrt(sigsq_d))
#
# x <- rbinom(n, size = 1,
#             prob = (1 + exp(-gammas[1] - gammas[2] * z - gammas[3] * c - gammas[4] * b))^(-1))
#
# y <- betas[1] + betas[2] * x + betas[3] * z + betas[4] * c + betas[5] * b +
#   rnorm(n, sd = sqrt(sigsq_e))
#
# all_data <- data.frame(y = y, x = x, z = z, c = c, b = b, d = d)
# all_data[1: n.e, 1] <- NA
# all_data[(n.e + 1): (n.e + n.m), 3] <- NA
# main <- internal <- external <- NULL
# y_var <- "y"
# x_var <- "x"
# z_var <- "z"
# d_vars <- "d"
# c_vars <- "c"
# b_vars <- "b"
# estimate_var <- FALSE
#
# fit <- ml_linear_logistic_linear(all_data = all_data,
#                                y_var = "y",
#                                x_var = "x",
#                                z_var = "z",
#                                d_vars = "d",
#                                c_vars = "c",
#                                b_vars = "b",
#                                control = list(trace = 1))

ml_linear_logistic_linear <- function(all_data = NULL,
                                      main = NULL,
                                      internal = NULL,
                                      external = NULL,
                                      y_var,
                                      x_var,
                                      z_var,
                                      d_vars = NULL,
                                      c_vars = NULL,
                                      b_vars = NULL,
                                      integrate_tol = 1e-8,
                                      integrate_tol_hessian = integrate_tol,
                                      estimate_var = FALSE,
                                      ...) {

  # Get dimension of D, C, and B
  kd <- length(d_vars)
  kc <- length(c_vars)
  kb <- length(b_vars)

  # Get covariate lists
  xdcb <- c(x_var, d_vars, c_vars, b_vars)
  xzdcb <- c(x_var, z_var, d_vars, c_vars, b_vars)
  xzdc <- c(x_var, z_var, d_vars, c_vars)

  # Subset various data types
  if (! is.null(all_data)) {
    main <- all_data[complete.cases(all_data[, c(y_var, xdcb)]) & is.na(all_data[, z_var]), ]
    internal <- all_data[complete.cases(all_data[, c(y_var, xzdcb)]), ]
    external <- all_data[is.na(all_data[, y_var]) & complete.cases(all_data[, xzdc]), ]
  }

  n.m <- ifelse(is.null(main), 0, nrow(main))
  some.m <- n.m > 0
  if (some.m) {
    y.m <- main[, y_var]
    x.m <- main[, x_var]
    onexcb.m <- as.matrix(cbind(rep(1, n.m), main[, c(x_var, c_vars, b_vars)]))
    onedcb.m <- as.matrix(cbind(rep(1, n.m), main[, c(d_vars, c_vars, b_vars)]))
    onedc.m <- as.matrix(cbind(rep(1, n.m), main[, c(d_vars, c_vars)]))
  }

  n.i <- ifelse(is.null(internal), 0, nrow(internal))
  some.i <- n.i > 0
  if (some.i) {
    y.i <- internal[, y_var]
    x.i <- internal[, x_var]
    z.i <- internal[, z_var]
    onexzcb.i <- as.matrix(cbind(rep(1, n.i), internal[, c(x_var, z_var, c_vars, b_vars)]))
    onezdcb.i <- as.matrix(cbind(rep(1, n.i), internal[, c(z_var, d_vars, c_vars, b_vars)]))
    onedc.i <- as.matrix(cbind(rep(1, n.i), internal[, c(d_vars, c_vars)]))
  }

  n.e <- ifelse(is.null(external), 0, nrow(external))
  some.e <- n.e > 0
  if (some.e) {
    x.e <- external[, x_var]
    z.e <- external[, z_var]
    onezdcb.e <- as.matrix(cbind(rep(1, n.e), external[, c(z_var, d_vars, c_vars, b_vars)]))
    onedc.e <- as.matrix(cbind(rep(1, n.e), external[, c(d_vars, c_vars)]))
  }

  # Get number of betas, gammas, and alphas
  n.betas <- 3 + kc + kb
  n.gammas <- 2 + kd + kc + kb
  n.alphas <- 1 + kd + kc

  # Get indices for parameters being estimated and create labels
  loc.betas <- 1: n.betas
  beta.labels <- paste("beta", c("0", x_var, z_var, c_vars, b_vars), sep = "_")

  loc.gammas <- (n.betas + 1): (n.betas + n.gammas)
  gamma.labels <- paste("gamma", c("0", z_var, d_vars, c_vars, b_vars), sep = "_")

  loc.alphas <- (n.betas + n.gammas + 1): (n.betas + n.gammas + n.alphas)
  alpha.labels <- paste("alpha", c("0", d_vars, c_vars), sep = "_")

  loc.sigsq_e <- n.betas + n.gammas + n.alphas + 1
  loc.sigsq_d <- n.betas + n.gammas + n.alphas + 2

  theta.labels <- c(beta.labels, gamma.labels, alpha.labels, "sigsq_e", "sigsq_d")

  # Likelihood function for full ML
  if (some.m) {
    lf.full <- function(y,
                        x,
                        z,
                        sigsq_e,
                        mu_z.dc,
                        sigsq_d,
                        beta_z,
                        gamma_z,
                        xcb.term,
                        dcb.term
                        ) {

      z <- matrix(z, nrow = 1)
      dens <- apply(z, 2, function(z) {

        # Transformation
        s <- z / (1 - z^2)

        # E(Y|X,Z,C,B)
        mu_y.xzcb <- xcb.term + beta_z * s

        # P(X=1|Z,D,C,B)
        p_x.zdcb <- (1 + exp(-dcb.term - gamma_z * s))^(-1)

        # f(Y,X,Z|D,C,B) = f(Y|X,Z,C,B) P(X=x|Z,D,C,B) f(Z|D,C)
        dnorm(y, mean = mu_y.xzcb, sd = sqrt(sigsq_e)) *
          dbinom(x, size = 1, prob = p_x.zdcb) *
          dnorm(s, mean = mu_z.dc, sd = sqrt(sigsq_d))

      })

      # Back-transformation
      out <- matrix(dens * (1 + z^2) / (1 - z^2)^2, ncol = ncol(z))
      return(out)

    }
  }

  # Log-likelihood function
  llf <- function(f.theta, estimating.hessian = FALSE) {

    # Extract parameters
    f.betas <- matrix(f.theta[loc.betas], ncol = 1)
    f.beta_z <- f.betas[3]

    f.gammas <- matrix(f.theta[loc.gammas], ncol = 1)
    f.gamma_z <- f.gammas[2]

    f.alphas <- matrix(f.theta[loc.alphas], ncol = 1)

    f.sigsq_e <- f.theta[loc.sigsq_e]
    f.sigsq_d <- f.theta[loc.sigsq_d]

    if (some.m) {

      # Likelihood for main study subjects:
      # L = \int_z f(Y|Z,C,B) P(X=x|Z,D,C,B) f(Z|D,C) dZ

      # Get integration tolerance
      int.tol <- ifelse(estimating.hessian, integrate_tol_hessian, integrate_tol)

      # Terms for integral
      xcb.terms <- onexcb.m %*% f.betas[-3]
      dcb.terms <- onedcb.m %*% f.gammas[-2]
      mu_z.dc <- onedc.m %*% f.alphas

      int.vals <- c()
      for (ii in 1: n.m) {

        # Perform integration
        int.ii <- cubature::hcubature(f = lf.full,
                                      tol = int.tol,
                                      lowerLimit = -1,
                                      upperLimit = 1,
                                      vectorInterface = TRUE,
                                      y = y.m[ii],
                                      x = x.m[ii],
                                      sigsq_e = f.sigsq_e,
                                      mu_z.dc = mu_z.dc[ii],
                                      sigsq_d = f.sigsq_d,
                                      beta_z = f.beta_z,
                                      gamma_z = f.gamma_z,
                                      xcb.term = xcb.terms[ii],
                                      dcb.term = dcb.terms[ii])
        int.vals[ii] <- int.ii$integral
        if (int.ii$integral == 0) {
          print(paste("Integral is 0 for ii = ", ii, sep = ""))
          print(f.theta)
          break
        }

      }

      ll.m <- sum(log(int.vals))

    } else {
      ll.m <- 0
    }

    if (some.i) {

      # Likelihood for internal validation subjects:
      # L = f(Y|X,Z,C,B) P(X=x|Z,D,C,B) f(Z|D,C)
      ll.i <- sum(
        dnorm(y.i, log = TRUE,
              mean = onexzcb.i %*% f.betas,
              sd = sqrt(f.sigsq_e)) +
          dbinom(x.i, log = TRUE, size = 1,
                 prob = (1 + exp(-onezdcb.i %*% f.gammas))^(-1)) +
          dnorm(z.i, log = TRUE,
                mean = onedc.i %*% f.alphas, sd = sqrt(f.sigsq_d))
      )

    } else {
      ll.i <- 0
    }

    if (some.e) {

      # Likelihood for external validation subjects:
      # L = P(X=x|Z,D,C,B) f(Z|D,C)
      ll.e <- sum(
        dbinom(x.e, log = TRUE, size = 1,
               prob = (1 + exp(-onezdcb.e %*% f.gammas))^(-1)) +
          dnorm(z.e, log = TRUE,
                mean = onedc.e %*% f.alphas, sd = sqrt(f.sigsq_d))
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
    extra.args$start <- c(rep(0.01, n.betas + n.gammas + n.alphas), 1, 1)
  }
  if (is.null(extra.args$lower)) {
    extra.args$lower <- c(rep(-Inf, n.betas + n.gammas + n.alphas), 1e-4, 1e-4)
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
