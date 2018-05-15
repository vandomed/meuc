#' Propensity Score Calibration for Linear Regression
#'
#' Implements propensity score calibration as described by Sturmer et al.
#' (\strong{Am. J. Epidemiol.} 2005). Requires validation data.
#'
#' The disease model is:
#'
#' E(Y) = beta_0 + beta_x X + beta_g G
#'
#' where G = P(X|\strong{C},\strong{Z}), with \strong{C} but not \strong{Z}
#' available in the main study.
#'
#' In a validation dataset with (X, \strong{C}, \strong{Z}), logistic regression
#' is used to obtain fitted probabilities for G as well as an error-prone
#' propensity score G* = P(X|\strong{C}). A linear model is fitted to map from
#' G* to E(G|G*). Finally, in the main study, G*'s are calculated, then E(G|G*),
#' and the disease model is fit for Y vs. (X, E(G|G*)).
#'
#'
#'
#'
#' @references
#' Sturmer, T., Schneeweiss, S., Avorn, J. and Glynn, R.J. (2005) "Adjusting
#' effect estimates for unmeasured confounding with validation data using
#' propensity score calibration." \emph{Am. J. Epidemiol.} \strong{162}(3):
#' 279-289.
#'
#'
# Basic propensity score calibration function
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

  # Try using data frames for everything

  # Regression models:

  # TDM: g[E(Y)] = beta_0 + beta_X X + beta_P e(X_GS)
  # MEM: E[e(X_GS)] = lambda_0 + lambda_X X + lambda_p e(X_EP)

  # Propensity score models:

  # e(GS) = P(X = 1 | X_GS) = (1 + exp(-gamma_0 - gamma_X^T X_GS))^(-1)
  # e(EP) = P(X = 1 | X_EP) = (1 + exp(-gamma_0* - gamma_X*^T X_EP))^(-1)

  # Using conditional expectation view of RC rather than algebraic view, because
  # Delta-method variance esimator wouldn't be valid anyway, and conditional is
  # more flexible.

  # Ensure that input datasets are data frames, not matrices
  if (! is.null(all_data) & class(all_data) == "matrix") {
    all_data <- as.data.frame(all_data)
  }
  if (! is.null(main) & class(main) == "matrix") {
    main <- as.data.frame(main)
  }
  if (! is.null(internal) & class(internal) == "matrix") {
    internal <- as.data.frame(internal)
  }
  if (! is.null(external) & class(external) == "matrix") {
    external <- as.data.frame(external)
  }

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
