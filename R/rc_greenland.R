#' Regression Calibration with Internal Validation Data and One Surrogate
#'
#' Implements an efficient version of regression calibration for
#' main study/internal validation designs with one surrogate, as described by
#' Spiegelman et al. (\emph{Stat. Med.} 2001). Uses ideas from Greenland
#' (\emph{Stat. Med.} 1988).
#'
#' The true disease model is a GLM:
#'
#' g[E(Y)] = beta_0 + beta_z Z + beta_c^T C + beta_b^T B
#'
#' The measurement error model is:
#'
#' E(Z) = alpha_d D + alpha_c^T C
#'
#' And the naive disease model is:
#'
#' g[E(Y)] = beta*_0 + beta*_Z D + \strong{beta*_C}^T \strong{C} +
#' \strong{beta*_B}^T \strong{B}
#'
#' The procedure involves obtaining two sets of \strong{beta} estimates: one
#' by fitting the true disease model with internal validation data, and the
#' other by treating the validation study as external and doing a main
#' study/external validation study correction. The two estimates are weighted by
#' the inveres of their estimated variances.
#'
#'
#' @inheritParams rc_algebraic
#'
#'
#' @return
#' List containing parameter estimates and variance estimates.
#'
#'
#' @references
#' Greenland, S. (1988) "Variance estimation for epidemiologic effect estimates
#' under misclassification." \emph{Stat. Med.} \strong{7}: 745-757.
#'
#' Spiegelman, D., Carroll, R.J., and Kipnis, V. (2001) "Efficient regression
#' calibration for logistic regression in main study/internal validation study
#' designs with an imperfect reference instrument." \emph{Stat. Med.}
#' \strong{20}: 139-160.
#'
#'
#' @export
# # Data for testing
# n.m <- 1000
# n.i <- 1000
# n <- n.m + n.i
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
# all_data[1: n.m, 2] <- NA
# main <- subset(all_data, is.na(z))
# internal <- subset(all_data, ! is.na(z))
# y_var <- "y"
# z_var <- "z"
# d_var <- "d"
# c_vars <- "c"
# b_vars <- NULL
# tdm_family <- "gaussian"
# beta_0_formula <- 1
#
# fit <- rc_algebraic(all_data = all_data,
#                     y_var = "y",
#                     z_var = "z",
#                     d_var = "d",
#                     c_vars = "c")
rc_greenland <- function(all_data = NULL,
                         main = NULL,
                         internal = NULL,
                         y_var,
                         z_var,
                         d_var,
                         c_vars = NULL,
                         b_vars = NULL,
                         tdm_family = "gaussian") {

  # Make sure dat is a matrix
  dat <- as.matrix(dat)

  # Create main and internal datasets if necessary
  if (! is.null(all_data)) {

    main <- subset(all_data, is.na(all_data[, z_var]))
    n.main <- nrow(main)

    internal <- subset(all_data, ! is.na(all_data[, z_var]))
    n.internal <- nrow(internal)

  }

  # Estimate betas using internal validation data
  tdm.formula <- paste(paste(y_var, " ~ ", sep = ""),
                       paste(c(z_var, c_vars, b_vars), collapse = " + "), sep = "")
  tdm.fit.i <- glm(tdm.formula, data = internal, family = tdm_family)
  theta.i <- tdm.fit.i$coef
  var.i <- vcov(tdm.fit.i)

  # Estimate betas using main study/pretend external validation data
  external <- internal
  external[, y_var] <- NA
  tdm.fit.e <- rc_algebraic(all_data = rbind(main, external),
                            y_var = y_var,
                            z_var = z_var,
                            d_var = d_var,
                            c_vars = c_vars,
                            tdm_family = tdm_family,
                            delta_var = TRUE)
  theta.e <- tdm.fit.e$theta.hat[1: length(theta.i)]
  var.e <- tdm.fit.e$delta.var

  # Determine weights, estimate parameters, and estimate variance
  c1 <- diag(var.e) / (diag(var.i) + diag(var.e))
  theta <- c1 * theta.i + (1 - c1) * theta.e
  theta.var <- as.matrix(diag(c1) %*% var.i %*% diag(c1) + diag(1 - c1) %*% var.e %*% diag(1 - c1))

  # Create and return ret.list
  ret.list <- list(theta = theta, theta.var = theta.var)
  return(ret.list)

}
