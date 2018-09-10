#' Propensity Score Calibration (Conditional Expectation Method)
#'
#' Implements the "conditional expectation" version of propensity score
#' calibration as described by Sturmer et al. (\emph{Am. J. Epidemiol.} 2005).
#' For the "algebraic" version, see \code{\link{psc_algebraic}}.
#'
#' The true disease model is a GLM:
#'
#' g[E(Y)] = beta_0 + beta_x X + beta_g G
#'
#' where G = P(X|\strong{C},\strong{Z}), with \strong{C} but not \strong{Z}
#' available in the main study.
#'
#' In a validation study with (X, \strong{C}, \strong{Z}), logistic regression
#' is used to obtain fitted probabilities for G as well as an error-prone
#' version G* = P(X|\strong{C}). A linear model is fitted to map from G* to
#' E(G|G*). Finally, in the main study, G*'s are calculated, then E(G|G*),
#' and the disease model is fit for Y vs. (X, E(G|G*)).
#'
#' @inheritParams rc_cond_exp
#'
#' @param x_var Character string specifying name of X variable.
#' @param gs_vars Character vector specifying names of variables for gold
#' standard propensity score.
#' @param ep_vars Character vector specifying names of variables for error-prone
#' propensity score.
#' @param surrogacy Logical value for whether to assume surrogacy, which means
#' that the error-prone propensity score is not informative of Y given X and the
#' gold standard propensity score. Have to assume surrogacy if validation data
#' is external.
#' @param ep_data Character string controlling what data is used to fit the
#' error-prone propensity score model. Choices are \code{"validation"} for
#' validation study data, \code{"all"} for main study and validation study data,
#' and \code{"separate"} for validation data for first step and main study data
#' for second step.
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
# n.m <- 5000
# n.e <- 5000
# n <- n.m + n.e
#
# gammas <- c(0, 0.25, 0.25)
# betas <- c(0, 0.25, 0.1, 0.2)
# sigsq_e <- 0.5
#
# c <- rnorm(n)
# z <- rnorm(n)
# x <- rbinom(n, size = 1, prob = (1 + exp(-gammas[1] - gammas[2] * c - gammas[3] * z))^(-1))
# y <- betas[1] + betas[2] * z + betas[3] * x + betas[4] * c + rnorm(n, sd = sqrt(sigsq_e))
#
# all_data <- data.frame(y = y, z = z, x = x, c = c)
# all_data[1: n.e, 1] <- NA
# all_data[(n.e + 1): n, 2] <- NA
# main <- internal <- external <- NULL
# y_var <- "y"
# x_var <- "x"
# gs_vars <- c("z", "c")
# ep_vars <- "c"
# tdm_family <- "gaussian"
# surrogacy <- TRUE
# ep_data <- "validation"
# beta_0_formula <- 1
# delta_var <- TRUE
# boot_var <- TRUE
# boots <- 100
#
# fit <- psc_cond_exp(all_data = all_data,
#                     y_var = "y",
#                     x_var = "x",
#                     gs_vars = c("z", "c"),
#                     ep_vars = "c")

psc_cond_exp <- function(all_data = NULL,
                         main = NULL,
                         internal = NULL,
                         external = NULL,
                         y_var,
                         x_var,
                         gs_vars,
                         ep_vars,
                         tdm_family = "gaussian",
                         surrogacy = TRUE,
                         ep_data = "separate",
                         boot_var = FALSE, boots = 100,
                         alpha = 0.05) {

  # Get full list of covariates
  covariates <- unique(c(x_var, gs_vars, ep_vars))

  # Get list of covariates subject to missingness
  covariates.missing <- setdiff(gs_vars, ep_vars)

  # If all_data not specified, create it from main, internal, and external
  if (is.null(all_data)) {

    if (! is.null(main)) {
      if (! all(covariates.missing %in% names(main))) {
        main[, covariates.missing] <- NA
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

  # Subset validation data with (X, Z, C) and data with outcome variable
  locs.val <- which(complete.cases(all_data[, covariates]))
  val_data <- all_data[locs.val, ]

  locs.y <- which(! is.na(all_data[, y_var]))
  y_data <- all_data[locs.y, ]

  locs.main <- which(! is.na(all_data[, y_var]) & ! complete.cases(all_data[, covariates]))
  main <- all_data[locs.main, ]

  # Figure out type of validation data
  val.type <- ifelse(length(locs.main) == length(locs.y), "external", "internal")

  # Fit error-prone propensity score model
  ep.formula <- paste(paste(x_var, "~"), paste(ep_vars, collapse = " + "))
  if (ep_data %in% c("validation", "separate")) {
    fit.ep <- glm(ep.formula, data = val_data, family = "binomial")
  } else if (ep_data == "all") {
    fit.ep <- glm(ep.formula, data = all_data, family = "binomial")
  }
  fitted.gstar <- predict(fit.ep, newdata = val_data, type = "response")

  # Fit gold standard propensity score model
  gs.formula <- paste(paste(x_var, "~"), paste(gs_vars, collapse = " + "))
  fit.gs <- glm(gs.formula, data = val_data, family = "binomial")
  fitted.g <- fit.gs$fitted

  # Fit MEM for G vs. (X, G*)
  fit.mem <- lm(fitted.g ~ val_data[, x_var] + fitted.gstar)

  # Add G*'s and G's to y_data for internal validation subjects
  y_data$gstar <- NA
  y_data$g <- NA
  if (val.type == "internal") {
    y_data$gstar[locs.val] <- fitted.gstar
    y_data$g[locs.val] <- fitted.g
  }

  # Re-fit error-prone propensity score model if requested
  if (ep_data == "separate") {
    fit.ep <- glm(ep.formula, data = main, family = "binomial")
  }

  # Calculate G*'s and E(G)'s for main study subjects
  y_data$gstar[locs.main] <- predict(fit.ep, newdata = main, type = "response")
  y_data$g[locs.main] <- as.matrix(cbind(1, y_data[locs.main, c(x_var, "gstar")])) %*% fit.mem$coef

  # Fit TDM and return beta estimate
  if (surrogacy) {
    tdm.formula <- paste(paste(y_var, "~"), paste(c(x_var, "g"), collapse = " + "))
    beta.labels <- paste("beta", c("0", x_var, "g"), sep = "_")
  } else {
    tdm.formula <- paste(paste(y_var, "~"), paste(c(x_var, "g", "gstar"), collapse = " + "))
    beta.labels <- paste("beta", c("0", x_var, "g", "gstar"), sep = "_")
  }
  tdm.fit <- glm(tdm.formula, data = y_data, family = tdm_family)
  beta.hat <- tdm.fit$coef
  names(beta.hat) <- beta.labels

  # Add beta.hat to ret.list
  ret.list <- list(beta.hat = beta.hat)

  # Get bootstrap variance estimate if requested
  if (boot_var) {

    # Various data types
    locs.m <- which(! is.na(all_data[, y_var]) & ! complete.cases(all_data[, covariates]))
    locs.i <- which(! is.na(all_data[, y_var]) & complete.cases(all_data[, covariates]))
    locs.e <- which(is.na(all_data[, y_var]) & complete.cases(all_data[, covariates]))

    # Initialize matrix for theta estimates
    beta.hat.boots <- matrix(NA, ncol = length(beta.hat), nrow = boots)

    # Bootstrap
    for (ii in 1: boots) {

      beta.hat.boots[ii, ] <- psc_cond_exp(
        all_data = all_data[c(sample(locs.m, replace = TRUE),
                              sample(locs.i, replace = TRUE),
                              sample(locs.e, replace = TRUE)), ],
        y_var = y_var,
        x_var = x_var,
        gs_vars = gs_vars,
        ep_vars = ep_vars,
        tdm_family = tdm_family,
        surrogacy = surrogacy,
        ep_data = ep_data
      )

    }

    # Calculate bootstrap variance estimates
    boot.variance <- var(beta.hat.boots)
    rownames(boot.variance) <- colnames(boot.variance) <- beta.labels
    ret.list$boot.var <- boot.variance

    boot.ci <- apply(beta.hat.boots, 2, function(x)
      quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
    colnames(boot.ci) <- beta.labels
    ret.list$boot.ci <- boot.ci

  }

  # Return ret.list
  if (length(ret.list) == 1) {
    ret.list <- ret.list[[1]]
  }
  return(ret.list)

}


# Archived 6/5/18 - was in middle of modifying

# psc <- function(all_data = NULL,
#                 main = NULL,
#                 internal = NULL,
#                 external = NULL,
#                 y_var,
#                 x_var,
#                 gs_vars,
#                 ep_vars,
#                 tdm_family = "gaussian",
#                 surrogacy = TRUE,
#                 ep_data = "validation",
#                 all_imputed = FALSE,
#                 boot_var = FALSE, boots = 100,
#                 alpha = 0.05) {
#
#   # Get full list of covariates
#   covariates <- unique(c(x_var, gs_vars, ep_vars))
#
#   # Get list of covariates subject to missingness
#   covariates.missing <- setdiff(gs_vars, ep_vars)
#
#   # If all_data not specified, create it from main, internal, and external
#   if (is.null(all_data)) {
#
#     if (! is.null(main)) {
#       if (! all(covariates.missing %in% names(main))) {
#         main[, covariates.missing] <- NA
#       }
#       main <- main[, c(y_var, covariates)]
#       all_data <- main
#     }
#
#     if (! is.null(internal)) {
#       internal <- internal[, c(y_var, covariates)]
#       all_data <- rbind(all_data, internal)
#     }
#
#     if (! is.null(external)) {
#       if (! y_var %in% names(external)) {
#         external[, y_var] <- NA
#       }
#       external <- external[, c(y_var, covariates)]
#       all_data <- rbind(all_data, external)
#     }
#
#   }
#
#   # Get validation data with (X, Z, C)
#   val_data <- all_data[complete.cases(all_data[, covariates]), ]
#
#   # Fit error-prone propensity score model
#   ep.formula <- paste(paste(x_var, "~"), paste(ep_vars, collapse = " + "))
#   if (ep_data %in% c("validation", "separate")) {
#     fit.ep <- glm(ep.formula, data = val_data, family = "binomial")
#   } else if (ep_data == "all") {
#     fit.ep <- glm(ep.formula, data = all_data, family = "binomial")
#   }
#
#   # Fit gold standard propensity score model
#   gs.formula <- paste(paste(x_var, "~"), paste(gs_vars, collapse = " + "))
#   fit.gs <- glm(gs.formula, data = val_data, family = "binomial")
#
#   # Fit MEM for G vs. (X, G*)
#   fitted.gstar <- predict(fit.ep, newdata = val_data, type = "response")
#   fitted.g <- fit.gs$fitted
#   fit.mem <- lm(fitted.g ~ val_data[, x_var] + fitted.gstar)
#
#   # Calculate G's for data with Y
#   y_data <- all_data[! is.na(all_data[, y_var]), ]
#
#   # Calculate G's for subjects with (C, Z), unless all_imputed is TRUE
#   y_data$g <- NA
#   if (! all_imputed) {
#     y_data$g <- predict(fit.gs, newdata = y_data, type = "response")
#   }
#
#   # Re-fit error-prone propensity score model if requested
#   if (ep_data == "separate") {
#     fit.ep <- glm(ep.formula, data = y_data, family = "binomial")
#   }
#
#   # Calculate E(G|X,G*)
#   locs <- which(is.na(y_data$g))
#   y_data$gstar <- NA
#   y_data$gstar[locs] <- predict(fit.ep, newdata = y_data[locs, ], type = "response")
#   y_data$g[locs] <- as.matrix(cbind(1, y_data[locs, c(x_var, "gstar")])) %*% fit.mem$coef
#
#   # Fit TDM and return beta estimate
#   if (surrogacy) {
#     tdm.formula <- paste(paste(y_var, "~"), paste(c(x_var, "g"), collapse = " + "))
#     beta.labels <- paste("beta", c("0", x_var, "g"), sep = "_")
#   } else {
#     tdm.formula <- paste(paste(y_var, "~"), paste(c(x_var, "g", "gstar"), collapse = " + "))
#     beta.labels <- paste("beta", c("0", x_var, "g", "gstar"), sep = "_")
#   }
#   tdm.fit <- glm(tdm.formula, data = y_data, family = tdm_family)
#   beta.hat <- tdm.fit$coef
#   names(beta.hat) <- beta.labels
#
#   # Add beta.hat to ret.list
#   ret.list <- list(beta.hat = beta.hat)
#
#   # Get bootstrap variance estimate if requested
#   if (boot_var) {
#
#     # Various data types
#     locs.m <- which(! is.na(all_data[, y_var]) & ! complete.cases(all_data[, covariates]))
#     locs.i <- which(! is.na(all_data[, y_var]) & complete.cases(all_data[, covariates]))
#     locs.e <- which(is.na(all_data[, y_var]) & complete.cases(all_data[, covariates]))
#
#     # Initialize matrix for theta estimates
#     beta.hat.boots <- matrix(NA, ncol = length(beta.hat), nrow = boots)
#
#     # Bootstrap
#     for (ii in 1: boots) {
#
#       beta.hat.boots[ii, ] <- psc(
#         all_data = all_data[c(sample(locs.m, replace = TRUE),
#                               sample(locs.i, replace = TRUE),
#                               sample(locs.e, replace = TRUE)), ],
#         y_var = y_var,
#         x_var = x_var,
#         gs_vars = gs_vars,
#         ep_vars = ep_vars,
#         tdm_family = tdm_family,
#         surrogacy = surrogacy,
#         ep_data = ep_data,
#         all_imputed = all_imputed
#       )
#
#     }
#
#     # Calculate bootstrap variance estimates
#     boot.variance <- var(beta.hat.boots)
#     rownames(boot.variance) <- colnames(boot.variance) <- beta.labels
#     ret.list$boot.var <- boot.variance
#
#     boot.ci <- apply(beta.hat.boots, 2, function(x)
#       quantile(x, probs = c(alpha / 2, 1 - alpha / 2)))
#     colnames(boot.ci) <- beta.labels
#     ret.list$boot.ci <- boot.ci
#
#   }
#
#   # Return ret.list
#   if (length(ret.list) == 1) {
#     ret.list <- ret.list[[1]]
#   }
#   return(ret.list)
#
# }
