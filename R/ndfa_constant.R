#' Normal Discriminant Function Approach for Estimating Odds Ratio with Exposure
#' Potentially Subject to Additive Normal Measurement Error (Constant Odds Ratio
#' Version)
#'
#' See \code{\link{ndfa}}.
#'
#'
#' @param y Numeric vector of Y values.
#' @param xtilde List of numeric vectors with Xtilde values.
#' @param c Numeric matrix with \strong{C} values (if any), with one row for
#' each subject. Can be a vector if there is only 1 covariate.
#' @param merror Logical value for whether there is measurement error.
#' @param ... Additional arguments to pass to \code{\link[stats]{nlminb}}.
#'
#'
#' @return
#' List containing:
#' \enumerate{
#'  \item Numeric vector of parameter estimates.
#'  \item Object returned by \code{\link[stats]{nlminb}} containing information
#'  about likelihood maximization.
#'  \item Akaike information criterion (AIC).
#' }
#'
#'
#' @references
#' Lyles, R.H., Guo, Y., and Hill, A.N. (2009) "A fresh look at the discriminant
#' function approach for estimating crude or adjusted odds ratios."
#' \emph{Am. Stat} \strong{63}(4): 320--327.
#'
#'
#' @export
# data(dat_ndfa)
# dat <- dat_ndfa$dat
# y <- dat$y
# xtilde <- dat_ndfa$reps
# c <- dat$c
# merror <- TRUE
ndfa_constant <- function(y,
                          xtilde,
                          c = NULL,
                          merror = FALSE,
                          ...) {

  # Check that inputs are valid
  if (! is.logical(merror)) {
    stop("The input 'merror' should be TRUE or FALSE.")
  }

  # Get name of y input
  y.varname <- deparse(substitute(y))
  if (length(grep("$", y.varname, fixed = TRUE)) > 0) {
    y.varname <- substr(y.varname,
                        start = which(unlist(strsplit(y.varname, "")) == "$") + 1,
                        stop = nchar(y.varname))
  }

  # Get information about covariates C
  if (is.null(c)) {
    c.varnames <- NULL
    n.cvars <- 0
  } else {
    c.varname <- deparse(substitute(c))
    if (class(c) != "matrix") {
      c <- as.matrix(c)
    }
    n.cvars <- ncol(c)
    c.varnames <- colnames(c)
    if (is.null(c.varnames)) {
      if (n.cvars == 1) {
        if (length(grep("$", c.varname, fixed = TRUE)) > 0) {
          c.varname <- substr(c.varname,
                              start = which(unlist(strsplit(c.varname, "")) == "$") + 1,
                              stop = nchar(c.varname))
        }
        c.varnames <- c.varname
      } else {
        c.varnames <- paste("c", 1: n.cvars, sep = "")
      }
    }
  }

  # Sample size
  n <- length(y)

  # Get number of gammas
  n.gammas <- 2 + n.cvars

  # Construct (1, Y, C) matrix
  oneyc <- cbind(rep(1, n), y, c)

  # If no measurement error and xtilde is a list, just use first measurements
  if (! merror & class(xtilde) == "list") {
    xtilde <- sapply(xtilde, function(x) x[1])
  }

  # Separate out subjects with replicates
  class.xtilde <- class(xtilde)
  if (class.xtilde == "list") {
    k <- sapply(xtilde, length)
    which.r <- which(k > 1)
    n.r <- length(which.r)
    some.r <- n.r > 0
    if (some.r) {

      # Replicates
      k.r <- k[which.r]
      oneyc.r <- oneyc[which.r, , drop = FALSE]
      xtilde.r <- xtilde[which.r]

    }
    n <- n - n.r
    some.s <- n > 0
    if (some.s) {

      # Singles
      oneyc <- oneyc[-which.r, , drop = FALSE]
      xtilde <- unlist(xtilde[-which.r])

    }
  } else {
    some.r <- FALSE
    some.s <- TRUE
  }

  # Get indices for parameters being estimated and create labels
  loc.gammas <- 1: n.gammas
  gamma.labels <- paste("gamma", c("0", y.varname, c.varnames), sep = "_")

  loc.sigsq <- n.gammas + 1
  loc.sigsq_m <- n.gammas + 2

  theta.labels <- c(gamma.labels, "sigsq")
  if (merror) {
    theta.labels <- c(theta.labels, "sigsq_m")
  }

  # Log-likelihood function
  llf <- function(f.theta) {

    # Extract parameters
    f.gammas <- matrix(f.theta[loc.gammas], ncol = 1)
    f.sigsq <- f.theta[loc.sigsq]

    if (merror) {
      f.sigsq_m <- f.theta[loc.sigsq_m]
    } else {
      f.sigsq_m <- 0
    }

    # Likelihood:
    # L = f(Xtilde|Y,C)

    if (some.r) {

      mu_xtilde.yc <- oneyc.r %*% f.gammas
      ll.r <- sum(
        mapply(
          FUN = function(k, mu_xtilde.yc, xtilde) {
            dmvnorm(x = xtilde, log = TRUE,
                    mean = rep(mu_xtilde.yc, k),
                    sigma = f.sigsq + diag(f.sigsq_m, k))
          },
          k = k.r,
          mu_xtilde.yc = mu_xtilde.yc,
          xtilde = xtilde.r
        )
      )

    } else {
      ll.r <- 0
    }

    if (some.s) {

      # E(Xtilde|Y,C)
      mu_xtilde.yc <- oneyc %*% f.gammas

      # Log-likelihood
      ll.s <- sum(dnorm(x = xtilde, log = TRUE,
                        mean = mu_xtilde.yc,
                        sd = sqrt(f.sigsq + f.sigsq_m)))

    } else {
      ll.s <- 0
    }

    # Return negative log-likelihood
    ll <- ll.r + ll.s
    return(-ll)

  }

  # Create list of extra arguments, and assign default starting values and
  # lower values if not specified by user
  extra.args <- list(...)
  if (is.null(extra.args$start)) {
    if (! merror) {
      extra.args$start <- c(rep(0.01, n.gammas), 1)
    } else {
      extra.args$start <- c(rep(0.01, n.gammas), 1, 1)
    }
  }
  if (is.null(extra.args$lower)) {
    if (! merror) {
      extra.args$lower <- c(rep(-Inf, n.gammas), 1e-3)
    } else {
      extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 2))
    }
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
  ml.estimates <- ml.max$par

  # Obtain point estimate for log-odds ratio
  gamma_y.hat <- ml.estimates[2]
  sigsq.hat <- ml.estimates[loc.sigsq]
  logOR.hat <- gamma_y.hat / sigsq.hat

  # Obtain variance estimates and perform bias adjustment
  hessian.mat <- pracma::hessian(f = llf, x0 = ml.estimates)
  theta.variance <- try(solve(hessian.mat), silent = TRUE)
  if (class(theta.variance) == "try-error") {
    message("Estimated Hessian matrix is singular, so variance-covariance matrix cannot be obtained and bias adjustment cannot be applied.")
    theta.variance <- NULL
    logOR.var <- logOR.var <- logOR_adj.hat <- logOR_adj.var <- NA
  } else {
    fprime <- matrix(c(1 / sigsq.hat, -gamma_y.hat / sigsq.hat^2), nrow = 1)
    colnames(theta.variance) <- rownames(theta.variance) <- theta.labels
    logOR.var <- fprime %*%
      theta.variance[c(2, loc.sigsq), c(2, loc.sigsq)] %*% t(fprime)
    sigsq.var <- theta.variance[loc.sigsq, loc.sigsq]
    logOR_adj.hat <- logOR.hat - gamma_y.hat * sigsq.var / sigsq.hat^3
    logOR_adj.var <- logOR.var * (logOR_adj.hat / logOR.hat)^2
    if (sign(logOR.hat) != sign(logOR_adj.hat)) {
      message("Bias adjustment flipped the sign of the log-OR estimate, so you may want to use the non-bias adjusted version.")
    }
  }

  # Create vector of estimates to return
  estimates <- c(ml.estimates,
                 logOR.hat, logOR.var,
                 logOR_adj.hat, logOR_adj.var)
  names(estimates) <- c(theta.labels,
                        "logOR.hat", "logOR.var",
                        "logOR_adj.hat", "logOR_adj.var")

  # Create list to return
  ret.list <- list(estimates = estimates,
                   theta.var = theta.variance,
                   nlminb.object = ml.max,
                   aic = 2 * (length(ml.estimates) + ml.max$objective))
  return(ret.list)

}
