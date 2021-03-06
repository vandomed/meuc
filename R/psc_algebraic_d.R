#' Propensity Score Calibration with Extra D Variable to Relax Surrogacy
#' Assumption
#'
#' Implements the "algebraic" version of propensity score calibration as
#' described by Sturmer et al. (\emph{Am. J. Epidemiol.} 2005), but with an
#' extra D variable in the MEM, which allows the error-prone propensity score
#' to be included in the TDM rather than assumed uninformative of the outcome.
#'
#' The true disease model is a GLM:
#'
#' g[E(Y)] = beta_0 + beta_x X + beta_g G + beta_gstar Gstar
#'
#' where G = P(X|\strong{C},\strong{Z}), with \strong{C} but not \strong{Z}
#' available in the main study, and Gstar = P(X|\strong{C}). Logistic regression
#' models are used to model P(X = 1|\strong{Z},\strong{C}) and
#' P(X = 1|\strong{C}):
#'
#' logit[P(X = 1)] = gamma_0 + \strong{gamma_z}^T \strong{Z} +
#' \strong{gamma_c}^T \strong{C}
#' logit[P(X = 1)] = gamma_0^* + \strong{gamma_c}^*T \strong{C}
#'
#' A linear model is used to map fitted error-prone propensity scores Gstar to
#' the expected gold standard propensity score G logistic regression):
#'
#' E(G) = lambda_0 + lambda_d D + lambda_x X + lambda_g Gstar
#'
#' There should be main study data with
#' (Y, X, D, \strong{C}) as well as external validation data with
#' (X, D, \strong{C}, \strong{Z}).
#'
#'
#' @inheritParams psc_algebraic
#' @param d_var Character string specifying name of D variable.
#'
#'
#' @references
#' Sturmer, T., Schneeweiss, S., Avorn, J. and Glynn, R.J. (2005) "Adjusting
#' effect estimates for unmeasured confounding with validation data using
#' propensity score calibration." \emph{Am. J. Epidemiol.} \strong{162}(3):
#' 279-289.
#'
#'
#' @export
# # Data for testing
# n.m <- 1000
# n.e <- 1000
# n <- n.m + n.e
#
# gammas <- c(0, 0.25, 0.25)
# betas <- c(0, 0.25, 0.1, 0.2)
# sigsq_e <- 0.5
#
# c <- rnorm(n)
# d <- rnorm(n)
# z <- 0.25 + 0.2 * d + rnorm(n, sd = sqrt(1))
# x <- rbinom(n, size = 1, prob = (1 + exp(-gammas[1] - gammas[2] * z - gammas[3] * c))^(-1))
# y <- betas[1] + betas[2] * x + betas[3] * z + betas[4] * c + rnorm(n, sd = sqrt(sigsq_e))
#
# all_data <- data.frame(y = y, x = x, z = z, d = d, c = c)
# all_data[1: n.e, 1] <- NA
# all_data[(n.e + 1): n, 3] <- NA
# main <- internal <- external <- NULL
# y_var <- "y"
# x_var <- "x"
# d_var <- "d"
# gs_vars <- c("z", "c")
# ep_vars <- "c"
# tdm_family <- "binomial"
# ep_data <- "validation"
# delta_var <- TRUE
#
# fit <- psc_algebraic_d(all_data = all_data,
#                        y_var = "y",
#                        x_var = "x",
#                        d_var = "d",
#                        gs_vars = c("z", "c"),
#                        ep_vars = "c",
#                        boot_var = TRUE)

psc_algebraic_d <- function(all_data = NULL,
                            main = NULL,
                            internal = NULL,
                            external = NULL,
                            y_var,
                            x_var,
                            d_var,
                            gs_vars,
                            ep_vars,
                            tdm_family = "gaussian",
                            ep_data = "validation",
                            delta_var = TRUE,
                            boot_var = FALSE, boots = 100,
                            alpha = 0.05) {

  # Get full list of covariates
  covariates <- unique(c(x_var, d_var, gs_vars, ep_vars))

  # Get list of covariates subject to missingness
  covariates.missing <- setdiff(gs_vars, ep_vars)

  # If all_data not specified, create it from main, internal, and external
  if (is.null(all_data)) {

    if (! all(covariates.missing %in% names(main))) {
      main[, covariates.missing] <- NA
    }
    main <- main[, c(y_var, covariates)]
    all_data <- main

    if (! y_var %in% names(external)) {
      external[, y_var] <- NA
    }
    external <- external[, c(y_var, covariates)]
    all_data <- rbind(all_data, external)

  }

  # Find validation data with (X, Z, D, C)
  locs.val <- which(complete.cases(all_data[, covariates]))
  val_data <- all_data[locs.val, ]

  # Fit error-prone propensity score model
  ep.formula <- paste(paste(x_var, "~"), paste(ep_vars, collapse = " + "))

  if (ep_data == "validation") {

    fit.ep <- glm(ep.formula, data = val_data, family = "binomial")
    all_data$gstar <- predict(fit.ep, newdata = all_data, type = "response")

  } else if (ep_data == "separate") {

    fit.ep.val <- glm(ep.formula, data = val_data, family = "binomial")
    fit.ep.main <- glm(ep.formula, data = all_data[-locs.val, ], family = "binomial")

    all_data$gstar <- NA
    all_data$gstar[locs.val] <- predict(fit.ep.val, newdata = val_data, type = "response")
    all_data$gstar[-locs.val] <- predict(fit.ep.main, newdata = all_data[-locs.val, ], type = "response")

  } else {

    fit.ep <- glm(ep.formula, data = all_data, family = "binomial")
    all_data$gstar <- predict(fit.ep, newdata = all_data, type = "response")

  }

  # Fit gold standard propensity score model
  gs.formula <- paste(paste(x_var, "~"), paste(gs_vars, collapse = " + "))
  fit.gs <- glm(gs.formula, data = val_data, family = "binomial")
  all_data$g <- predict(fit.gs, newdata = all_data, type = "response")

  # Call rc_algebraic
  fit.rc <- rc_algebraic(
    all_data = all_data,
    y_var = y_var,
    z_var = "g",
    d_var = d_var,
    c_vars = c(x_var, "gstar"),
    tdm_family = tdm_family,
    delta_var = delta_var
  )
  if (delta_var) {
    ret.list <- list(theta.hat = fit.rc$theta.hat,
                     delta.var = fit.rc$delta.var)
  } else {
    ret.list <- list(theta.hat = fit.rc)
  }

  # Get bootstrap variance estimate if requested
  if (boot_var) {

    # Various data types
    locs.m <- which(! is.na(all_data[, y_var]) & ! complete.cases(all_data[, covariates]))
    locs.i <- which(! is.na(all_data[, y_var]) & complete.cases(all_data[, covariates]))
    locs.e <- which(is.na(all_data[, y_var]) & complete.cases(all_data[, covariates]))

    # Initialize matrix for theta estimates
    theta.hat.boots <- matrix(NA, ncol = length(ret.list$theta.hat), nrow = boots)

    # Bootstrap
    for (ii in 1: boots) {

      theta.hat.boots[ii, ] <- psc_algebraic_d(
        all_data = all_data[c(sample(locs.m, replace = TRUE),
                              sample(locs.i, replace = TRUE),
                              sample(locs.e, replace = TRUE)), ],
        y_var = y_var,
        x_var = x_var,
        d_var = d_var,
        gs_vars = gs_vars,
        ep_vars = ep_vars,
        tdm_family = tdm_family,
        ep_data = ep_data,
        delta_var = FALSE
      )

    }

    # Calculate bootstrap variance estimates
    boot.variance <- var(theta.hat.boots)
    rownames(boot.variance) <- colnames(boot.variance) <- names(ret.list$theta.hat)
    ret.list$boot.var <- boot.variance

    boot.ci <- apply(theta.hat.boots, 2, function(x)
      quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
    colnames(boot.ci) <- names(ret.list$theta.hat)
    ret.list$boot.ci <- boot.ci

  }

  # Return ret.list
  if (length(ret.list) == 1) {
    ret.list <- ret.list[[1]]
  }
  return(ret.list)

}



# psc_algebraic_d <- function(all_data = NULL,
#                             main = NULL,
#                             internal = NULL,
#                             external = NULL,
#                             y_var,
#                             x_var,
#                             d_var,
#                             gs_vars,
#                             ep_vars,
#                             tdm_family = "gaussian",
#                             ep_data = "validation",
#                             delta_var = TRUE,
#                             boot_var = FALSE, boots = 100) {
#
#   # Get full list of covariates
#   covariates <- unique(c(x_var, d_var, gs_vars, ep_vars))
#
#   # Get list of covariates subject to missingness
#   covariates.missing <- setdiff(gs_vars, ep_vars)
#
#   # If all_data not specified, create it from main, internal, and external
#   if (is.null(all_data)) {
#
#     if (! all(covariates.missing %in% names(main))) {
#       main[, covariates.missing] <- NA
#     }
#     main <- main[, c(y_var, covariates)]
#     all_data <- main
#
#     if (! y_var %in% names(external)) {
#       external[, y_var] <- NA
#     }
#     external <- external[, c(y_var, covariates)]
#     all_data <- rbind(all_data, external)
#
#   }
#
#   # Find validation data with (X, Z, D, C)
#   locs.val <- which(complete.cases(all_data[, covariates]))
#   val_data <- all_data[locs.val, ]
#
#   # Fit error-prone propensity score model
#   ep.formula <- paste(paste(x_var, "~"), paste(ep_vars, collapse = " + "))
#
#   if (ep_data == "validation") {
#
#     fit.ep <- glm(ep.formula, data = val_data, family = "binomial")
#     all_data$gstar <- predict(fit.ep, newdata = all_data, type = "response")
#
#   } else if (ep_data == "separate") {
#
#     fit.ep.val <- glm(ep.formula, data = val_data, family = "binomial")
#     fit.ep.main <- glm(ep.formula, data = all_data[-locs.val, ], family = "binomial")
#
#     all_data$gstar <- NA
#     all_data$gstar[locs.val] <- predict(fit.ep.val, newdata = val_data, type = "response")
#     all_data$gstar[-locs.val] <- predict(fit.ep.main, newdata = all_data[-locs.val, ], type = "response")
#
#   } else {
#
#     fit.ep <- glm(ep.formula, data = all_data, family = "binomial")
#     all_data$gstar <- predict(fit.ep, newdata = all_data, type = "response")
#
#   }
#
#   # Fit gold standard propensity score model
#   gs.formula <- paste(paste(x_var, "~"), paste(gs_vars, collapse = " + "))
#   fit.gs <- glm(gs.formula, data = val_data, family = "binomial")
#   all_data$g <- predict(fit.gs, newdata = all_data, type = "response")
#
#   # Call rc_algebraic
#   fit.rc <- rc_algebraic(
#     all_data = all_data,
#     y_var = y_var,
#     z_var = "g",
#     d_var = d_var,
#     c_vars = c(x_var, "gstar"),
#     tdm_family = tdm_family,
#     delta_var = delta_var
#   )
#
#   return(fit.rc)
#
# }
