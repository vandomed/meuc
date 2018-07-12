#' Regression Calibration (Algebraic Method)
#'
#' Implements the "algebraic" version of regression calibration as described by
#' Rosner et al. (\emph{Stat. Med.} 1989). For the "conditional expectation"
#' version, see \code{\link{rc_cond_exp}}.
#'
#' The true disease model is a GLM:
#'
#' g[E(Y)] = beta_0 + beta_z Z + \strong{beta_c}^T \strong{C} +
#' \strong{beta_b}^T \strong{B}
#'
#' The measurement error model is:
#'
#' E(Z) = alpha_d D + \strong{alpha_c}^T \strong{C}
#'
#' And the naive disease model is:
#'
#' g[E(Y)] = beta*_0 + beta*_Z D + \strong{beta*_C}^T \strong{C} +
#' \strong{beta*_B}^T \strong{B}
#'
#' The procedure involves fitting the naive disease model using main study
#' data, fitting the measurement error model using validation data, and
#' solving a system of equations to get the regression calibration estimates.
#'
#'
#' @inheritParams rc_cond_exp
#'
#' @param d_var Character string specifying name of D variable.
#' @param beta_0_formula If \code{1}, formula for true disease model intercept
#' is:
#'
#' beta_0.hat = betastar_0.hat - alpha_0.hat beta_Z.hat
#'
#' If \code{2}, formula is:
#'
#' beta_0.hat = betastar_0.hat - alpha_0.hat beta_Z.hat - 1/2 beta_Z.hat^2
#' sigma_delta^2
#'
#' Formula 1 yields the same beta_0.hat as the "conditional expectation" view of
#' regression calibration (see \code{\link{rc_cond_exp}}). When the disease
#' model is logistic regression and the measurement error model is linear
#' regression, formula 1 is appropriate if Y is rare and Z|(D,\strong{C}) is
#' normal, and formula 2 is appropriate if beta_Z^2 sigma_delta^2 is small
#' (Kuha, Stat. Med. 1994). If neither criteria is met, regression calibration
#' may be unreliable.
#'
#' @param delta_var Logical value for whether to calculate a Delta method
#' variance-covariance matrix.
#'
#'
#' @return
#' If no variance estimates are requested, a named numeric vector of parameter
#' estimates. If one or more variance estimates are requested, a list that also
#' contains a variance-covariance matrix for each variance estimator.
#'
#'
#' @inherit rc_cond_exp references
#'
#'
#' @export
#'
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
# y_var <- "y"
# z_var <- "z"
# d_var <- "d"
# c_vars <- "c"
# b_vars <- NULL
# tdm_covariates <- mem_covariates <- NULL
# tdm_family <- "gaussian"
# beta_0_formula <- 1
# delta_var <- TRUE
# boot_var <- TRUE
# boots <- 100
# fit <- rc_algebraic(all_data = all_data,
#                     y_var = "y",
#                     z_var = "z",
#                     d_var = "d",
#                     c_vars = "c")

rc_algebraic <- function(all_data = NULL,
                         main = NULL,
                         internal = NULL,
                         external = NULL,
                         y_var,
                         z_var,
                         d_var = NULL,
                         c_vars = NULL,
                         b_vars = NULL,
                         tdm_covariates = NULL,
                         tdm_family = "gaussian",
                         mem_covariates = NULL,
                         beta_0_formula = 1,
                         delta_var = TRUE,
                         boot_var = FALSE, boots = 100,
                         alpha = 0.05) {

  # If beta_0_formula = 2, check that it is reasonable
  if (beta_0_formula == 2) {
    if (tdm_family != "binomial") {
      warning("The beta_0_formula = 2 intercept is appropriate when the true
              disease model is logistic regression and the outcome is rare. The
              tdm_family you selected is not logistic regression, so, depending
              on your specific scenario, you may want to re-run with
              beta_0_formula = 1.")
    }
  }

  # If tdm_covariates and mem_covariates specified, figure out d_var, c_vars,
  # and b_vars
  if (! is.null(tdm_covariates) & ! is.null(mem_covariates)) {
    tdm_covariates <- setdiff(tdm_covariates, z_var)
    d_var <- setdiff(mem_covariates, tdm_covariates)
    c_vars <- intersect(tdm_covariates, mem_covariates)
    b_vars <- setdiff(tdm_covariates, mem_covariates)
  }

  # Get full list of covariates
  covariates <- c(z_var, d_var, c_vars, b_vars)

  # Get dimension of C and B
  kc <- length(c_vars)
  kb <- length(b_vars)

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

  # Fit naive TDM to get betastar.hat vector
  naive.formula <- paste(paste(y_var, " ~ ", sep = ""),
                         paste(c(d_var, c_vars, b_vars), collapse = " + "),
                         sep = "")
  naive.fit <- glm(naive.formula, data = all_data, family = tdm_family)
  betastar.hat <- naive.fit$coef

  # Fit MEM using all available data to get alpha.hat vector
  mem.formula <- paste(paste(z_var, " ~ ", sep = ""),
                       paste(c(d_var, c_vars), collapse = " + "),
                       sep = "")
  mem.fit <- lm(mem.formula, data = all_data)
  sigsq_d.hat <- rev(anova(mem.fit)$"Mean Sq")[1]
  alpha.hat <- mem.fit$coef

  # Create labels for parameter estimates
  beta.labels <- c("beta_0", paste(rep("beta_", 1 + kc + kb),
                                   c(z_var, c_vars, b_vars),
                                   sep = ""))
  alpha.labels <- c("alpha_0", paste(rep("alpha_", 1 + kc),
                                     c(d_var, c_vars),
                                     sep = ""))
  theta.labels <- c(beta.labels, alpha.labels, "sigsq_d")

  # theta = (beta^T, alpha^T, sigsq_d)^T = f(betastar, alpha, sigsq_d)
  f <- function(x) {

    # Extract betastar, alpha, and sigsq_d
    f.betastar <- x[1: length(betastar.hat)]
    f.alpha <- x[(length(betastar.hat) + 1): (length(x) - 1)]
    f.sigsq_d <- x[length(x)]

    # Calculate f.beta depending on value of beta_0_formula input
    if (beta_0_formula == 1) {

      # beta_0 = betastar_0 - alpha_0 beta_Z

      # Construct A matrix
      f.A <- cbind(c(1, rep(0, 1 + kc + kb)),
                   c(f.alpha, rep(0, kb)),
                   rbind(matrix(0, nrow = 2, ncol = kc + kb), diag(kc + kb)))

      # Calculate beta = A^(-1) betastar
      f.beta <- solve(f.A) %*% f.betastar

    } else if (beta_0_formula == 2) {

      # beta_0 = betastar_0 - alpha_0 beta_Z - 1/2 beta_Z^2 sigsq_d

      # Construct A matrix
      f.A <- cbind(c(f.alpha[-1], rep(0, kb)),
                   rbind(matrix(0, nrow = 1, ncol = kc + kb), diag(kc + kb)))

      # Calculate (beta_Z, beta_C^T, beta_B^T) = A^(-1) betastar[-1]
      f.beta.nointercept <- solve(f.A) %*% f.betastar[-1]

      # Calculate beta_0
      f.beta_Z <- f.beta.nointercept[1]
      f.beta_0 <- f.betastar[1] - f.alpha[1] * f.beta_Z -
        1/2 * f.beta_Z^2 * f.sigsq_d

      # Construct beta = (beta_0, beta_Z, beta_C^T, beta_B^T)
      f.beta <- c(f.beta_0, f.beta.nointercept)

    }

    # Return theta = (beta^T, alpha^T, sigsq_d)^T
    f.theta <- c(f.beta, f.alpha, f.sigsq_d)
    return(f.theta)

  }

  # Obtain point estimate for theta, and add to ret.list
  theta.hat <- f(c(betastar.hat, alpha.hat, sigsq_d.hat))
  names(theta.hat) <- theta.labels
  ret.list <- list(theta.hat = theta.hat)

  # Calculate Delta-method variance estimate if requested
  if (delta_var) {

    # Estimate f'(betastar.hat, alpha.hat, sigsq_d)
    fprime <- jacobian(func = f,
                       x = c(betastar.hat, alpha.hat, sigsq_d.hat))

    # Construct V.hat(betastar.hat, alpha.hat, sigsq_d.hat)
    Sigma <- bdiag(vcov(naive.fit), vcov(mem.fit),
                   2 * sigsq_d.hat^2 / mem.fit$df.residual)

    # Calculate f'(.) V(.) f'(.)^T
    delta.variance <- as.matrix(fprime %*% Sigma %*% t(fprime))

    # Attach labels and add to ret.list
    colnames(delta.variance) <- rownames(delta.variance) <- theta.labels
    ret.list$delta.var <- delta.variance

  }

  # Calculate bootstrap variance estimate if requested
  if (boot_var) {

    # Various data types
    locs.m <- which(complete.cases(all_data[, c(y_var, mem_covariates)]) &
                      is.na(all_data[, z_var]))
    locs.i <- which(complete.cases(all_data[, c(y_var, covariates)]))
    locs.e <- which(complete.cases(all_data[, c(z_var, mem_covariates)]) &
                      is.na(all_data[, y_var]))

    # Initialize matrix for theta estimates
    theta.hat.boots <- matrix(NA, ncol = length(theta.hat), nrow = boots)

    # Bootstrap
    for (ii in 1: boots) {

      theta.hat.boots[ii, ] <- rc_algebraic(
        all_data = all_data[c(sample(locs.m, replace = TRUE),
                              sample(locs.i, replace = TRUE),
                              sample(locs.e, replace = TRUE)), ],
        y_var = y_var,
        z_var = z_var,
        d_var = d_var,
        c_vars = c_vars,
        b_vars = b_vars,
        tdm_family = tdm_family,
        beta_0_formula = beta_0_formula)

    }

    # Calculate bootstrap variance estimates
    boot.variance <- var(theta.hat.boots)
    rownames(boot.variance) <- colnames(boot.variance) <- theta.labels
    ret.list$boot.var <- boot.variance

    boot.ci <- apply(theta.hat.boots, 2, function(x)
      quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
    colnames(boot.ci) <- theta.labels
    ret.list$boot.ci <- boot.ci

  }

  # Convert ret.list to vector if length is 1
  if (length(ret.list) == 1) {
    ret.list <- ret.list[[1]]
  }

  # Return object
  return(ret.list)

}
