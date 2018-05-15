#' Propensity Score Calibration
#'
#' Implements propensity score calibration as described by Sturmer et al.
#' (\emph{Am. J. Epidemiol.} 2005). Requires validation data.
#'
#' The disease model is a GLM:
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
#' gold standard propensity score. If validation data is external, have to
#' assume surrogacy.
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
# n.e <- 100
# n <- n.m + n.e
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
# all_data[(n.e + 1): n, 2] <- NA
# main <- internal <- external <- NULL
# y_var <- "y"
# z_var <- "z"
# d_var <- "d"
# c_vars <- "c"
# b_vars <- NULL
# tdm_covariates <- mem_covariates <- NULL
# tdm_family <- "gaussian"
# mem_family <- "gaussian"
# beta_0_formula <- 1
# delta_var <- TRUE
# boot_var <- TRUE
# boots <- 100
psc <- function(all_data = NULL,
                main = NULL,
                internal = NULL,
                external = NULL,
                y_var,
                x_var,
                gs_vars,
                ep_vars,
                tdm_family = "gaussian",
                surrogacy = TRUE,
                weighted = FALSE,
                ep_refit = FALSE) {

  # Using conditional expectation view of RC rather than algebraic view, because
  # Delta-method variance esimator wouldn't be valid anyway, and conditional is
  # more flexible.

  # Get full list of covariates
  covariates <- unique(c(x_var, gs_vars, ep_vars))

  # Get list of covariates subject to missingness
  covariates.missing <- setdiff(gs_vars, ep_vars)

  # covariates <- unique(c(z.var, surrogate.var, tdm.covariates))

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



  # Fit MEM using all available data to get alpha.hat vector
  mem.formula <- paste(paste(z.var, " ~ ", sep = ""), paste(c(surrogate.var, tdm.covariates), collapse = " + "), sep = "")
  mem.fit <- glm(mem.formula, data = all_data, family = mem.family)
  alpha.hat <- mem.fit$coef

  # Fit naive TDM to get betastar.hat vector
  tdm.naive.formula <- paste(paste(y_var, " ~ ", sep = ""), paste(c(surrogate.var, tdm.covariates), collapse = " + "), sep = "")
  tdm.naive.fit <- glm(tdm.naive.formula, data = all_data, family = tdm_family)
  betastar.hat <- tdm.naive.fit$coef


  # Create matrix and data frame versions of all_data
  dat.mat <- as.matrix(dat)
  dat.df <- as.data.frame(dat)

  # Fit propensity score models using validation data, whether it is internal, external, or both
  locs.val <- which(complete.cases(dat.mat[, x.gs.columns]))
  locs.main <- which(complete.cases(dat.mat[, x.ep.columns]) & !complete.cases(dat.mat[, x.gs.columns]) & !is.na(dat.mat[, y.column]))
  val.mat <- dat.mat[locs.val, ]
  ps.ep.fit <- glm(val.mat[, a.column] ~ val.mat[, x.ep.columns], family = "binomial")
  ps.gs.fit <- glm(val.mat[, a.column] ~ val.mat[, x.gs.columns], family = "binomial")

  # Fit MEM relating gold-standard propensity score to exposure and error-prone propensity score
  mem.fit <- lm(ps.gs.fit$fitted ~ val.mat[, a.column] + ps.ep.fit$fitted)

  # Create ps.ep variable for all observations
  dat.df$ps.ep <- (1 + exp(- cbind(rep(1, nrow(dat.mat)), dat.mat[, x.ep.columns]) %*% ps.ep.fit$coef))^-1

  # Create ps.gs variable, directly from ps.gs.fit where possible
  dat.df$ps.gs <- mem.fit$coef[1] + mem.fit$coef[2] * dat.df[, a.column] + mem.fit$coef[3] * dat.df$ps.ep
  dat.df$ps.gs[locs.val] <- ps.gs.fit$fitted

  # Fit TDM and return beta estimate
  if (surrogacy) {
    if (weighted) {
      dat.new <- data.frame(y = c(dat.df[locs.main, y.column], dat.df[locs.val, y.column]),
                            z = c(rep(NA, length(locs.main)), dat.df$ps.gs[locs.val]),
                            c = c(dat.df[locs.main, a.column], dat.df[locs.val, a.column]),
                            d = c(dat.df$ps.ep[locs.main], dat.df$ps.ep[locs.val]))
      beta <- rc.greenland(dat = dat.new, y.column = 1, z.column = 2, c.columns = 3, d.column = 4)$theta[c(1, 3, 2)]
      beta.labels <- c("b.0", "b.a", "b.gs")
      names(beta) <- beta.labels
    } else {
      tdm.fit <- glm(dat.df[, y.column] ~ dat.df[, a.column] + dat.df$ps.gs, family = tdm_family)
      beta <- tdm.fit$coef
      beta.labels <- c("b.0", "b.a", "b.gs")
      names(beta) <- beta.labels
    }
  } else {
    tdm.fit <- glm(dat.df[, y.column] ~ dat.df[, a.column] + dat.df$ps.gs + dat.df$ps.ep, family = tdm_family)
    beta <- tdm.fit$coef
    beta.labels <- c("b.0", "b.a", "b.gs", "b.ep")
    names(beta) <- beta.labels
  }
  return(beta)

}
