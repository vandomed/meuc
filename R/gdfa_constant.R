#' Gamma Discriminant Function Approach for Estimating Odds Ratio with Exposure
#' Potentially Subject to Multiplicative Lognormal Measurement Error (Constant
#' Odds Ratio Version)
#'
#' See \code{\link{gdfa}}.
#'
#'
#' @param y Numeric vector of Y values.
#' @param xtilde Numeric vector (or list of numeric vectors, if there are
#' replicates) of Xtilde values.
#' @param c Numeric matrix with \strong{C} values (if any), with one row for
#' each subject. Can be a vector if there is only 1 covariate.
#' @param merror Logical value for whether there is measurement error.
#' @param integrate_tol Numeric value specifying the \code{tol} input to
#' \code{\link[cubature]{hcubature}}.
#' @param integrate_tol_hessian Same as \code{integrate_tol}, but for use when
#' estimating the Hessian matrix only. Sometimes more precise integration
#' (i.e. smaller tolerance) helps prevent cases where the inverse Hessian is not
#' positive definite.
#' @param estimate_var Logical value for whether to return variance-covariance
#' matrix for parameter estimates.
#' @param fix_posdef Logical value for whether to repeatedly reduce
#' \code{integrate_tol_hessian} by factor of 5 and re-estimate Hessian to try
#' to avoid non-positive definite variance-covariance matrix.
#' @param ... Additional arguments to pass to \code{\link[stats]{nlminb}}.
#'
#'
#' @return List containing:
#' \enumerate{
#' \item Numeric vector of parameter estimates.
#' \item Variance-covariance matrix.
#' \item Returned \code{\link[stats]{nlminb}} object from maximizing the
#' log-likelihood function.
#' \item Akaike information criterion (AIC).
#' }
#'
#'
#' @references
#' Whitcomb, B.W., Perkins, N.J., Zhang, Z., Ye, A., and Lyles, R. H. (2012)
#' "Assessment of skewed exposure in case-control studies with pooling."
#' \emph{Stat. Med.} \strong{31}: 2461--2472.
#'
#'
#' @export
gdfa_constant <- function(y,
                          xtilde,
                          c = NULL,
                          merror = TRUE,
                          integrate_tol = 1e-8,
                          integrate_tol_hessian = integrate_tol,
                          estimate_var = TRUE,
                          fix_posdef = TRUE,
                          ...) {

  # Check that inputs are valid
  if (! is.null(constant_or) && ! is.logical(constant_or)) {
    stop("The input 'contant_or' should be set to TRUE, FALSE, or NULL.")
  }
  if (! is.logical(merror)) {
    stop("The input 'merror' should be TRUE or FALSE.")
  }
  if (! (is.numeric(integrate_tol) & inside(integrate_tol, c(1e-32, Inf)))) {
    stop("The input 'integrate_tol' must be a numeric value greater than 1e-32.")
  }
  if (! (is.numeric(integrate_tol_hessian) & inside(integrate_tol_hessian, c(1e-32, Inf)))) {
    stop("The input 'integrate_tol_hessian' must be a numeric value greater than 1e-32.")
  }
  if (! is.logical(estimate_var)) {
    stop("The input 'estimate_var' should be TRUE or FALSE.")
  }
  if (! is.logical(fix_posdef)) {
    stop("The input 'fix_posdef' should be TRUE or FALSE.")
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
  n.gammas <- 1 + n.cvars

  # Construct (1, C) matrix
  onec <- cbind(rep(1, n), c)

  # If no measurement error and xtilde is a list, just use first measurements
  if (! merror & class(xtilde) == "list") {
    xtilde <- sapply(xtilde, function(x) x[1])
  }

  # Separate into subjects with precisely measured X, replicate Xtilde's, and
  # single imprecise Xtilde
  if (merror) {

    which.p <- NULL

    class.xtilde <- class(xtilde)
    if (class.xtilde == "list") {
      k <- sapply(xtilde, length)
    } else {
      k <- rep(1, n)
    }

    which.r <- which(k > 1)
    n.r <- length(which.r)
    some.r <- n.r > 0
    if (some.r) {
      k.r <- k[which.r]
      y.r <- y[which.r]
      onec.r <- onec[which.r, , drop = FALSE]
      xtilde.r <- xtilde[which.r]
    }

    which.i <- which(k == 1)
    n.i <- length(which.i)
    some.i <- n.i > 0
    if (some.i) {
      y.i <- y[which.i]
      onec.i <- onec[which.i, , drop = FALSE]
      xtilde.i <- unlist(xtilde[which.i])
    }

  } else {

    which.p <- 1: n
    some.p <- TRUE
    y.p <- y
    onec.p <- onec
    x.p <- xtilde

    some.r <- FALSE
    some.i <- FALSE

  }

  # Get indices for parameters being estimated and create labels
  loc.gammas <- 1: n.gammas
  gamma.labels <- paste("gamma", c("0", c.varnames), sep = "_")

  loc.bs <- (n.gammas + 1): (n.gammas + 2)
  loc.sigsq_m <- n.gammas + 3

  theta.labels <- c(gamma.labels, "b1", "b0")
  if (merror) {
    theta.labels <- c(theta.labels, "sigsq_m")
  }

  # Likelihood function for singles and replicates
  if (some.i | some.r) {

    lf <- function(xtilde,
                   x,
                   shape,
                   scale,
                   sigsq_m) {

      # f(XtildeX|Y,C_1,...,C_g)
      x <- matrix(x, nrow = 1)
      f_xtildex.yc <- apply(x, 2, function(z) {

        # Transformation
        s <- z / (1 - z)

        # Density
        1 / prod(xtilde) *
          prod(dnorm(x = log(xtilde),
                     mean = log(s) - 1/2 * sigsq_m,
                     sd = sqrt(sigsq_m))) *
          dgamma(x = s, shape = shape, scale = scale)

      })

      # Back-transformation
      out <- matrix(f_xtildex.yc / (1 - x)^2, ncol = ncol(x))

    }

  }

  # Log-likelihood function
  llf <- function(f.theta, estimating.hessian = FALSE) {

    # Extract parameters
    f.gammas <- matrix(f.theta[loc.gammas], ncol = 1)
    f.b1 <- f.theta[loc.bs[1]]
    f.b0 <- f.theta[loc.bs[2]]

    if (merror) {
      f.sigsq_m <- f.theta[loc.sigsq_m]
    } else {
      f.sigsq_m <- 0
    }

    if (some.p) {

      # Likelihood for subjects with precisely measured X:
      # L = f(X|Y,C)
      ll.p <- sum(dgamma(x = x.p, log = TRUE,
                         shape = exp(onec.p %*% f.gammas),
                         scale = ifelse(y.p == 1, f.b1, f.b0)))

    } else {
      ll.p <- 0
    }

    # Set skip.rest flag to FALSE
    skip.rest <- FALSE

    if (some.r) {

      # Likelihood for subjects with replicates:
      # L = int_X f(Xtilde|X) f(X|Y,C) dX

      # Shape and scale parameters to feed integral
      alphas <- exp(onec.r %*% f.gammas)
      betas <- ifelse(y.r == 1, f.b1, f.b0)

      # Get integration tolerance
      if (estimating.hessian) {
        int_tol <- integrate_tol_hessian
      } else {
        int_tol <- integrate_tol
      }

      int.vals <- c()
      for (ii in 1: length(xtilde.r)) {

        int.ii <- hcubature(
          f = lf,
          tol = integrate_tol,
          vectorInterface = TRUE,
          lowerLimit = 0, upperLimit = 1,
          xtilde = xtilde.r[[ii]],
          shape = alphas[ii],
          scale = betas[ii],
          sigsq_m = f.sigsq_m
        )
        int.vals[ii] <- int.ii$integral

        # If integral 0, set skip.rest to TRUE to skip further LL calculations
        if (int.ii$integral == 0) {
          print(paste("Integral is 0 for ii = ", ii, sep = ""))
          print(f.theta)
          skip.rest <- TRUE
          break
        }

      }
      ll.r <- sum(log(int.vals))

    } else {
      ll.r <- 0
    }

    if (some.i & ! skip.rest) {

      # Likelihood for subjects with single Xtilde:
      # L = int_X f(Xtilde|X) f(X|Y,C) dX

      # Shape and scale parameters to feed integral
      alphas <- exp(onec.i %*% f.gammas)
      betas <- ifelse(y.i == 1, f.b1, f.b0)

      # Get integration tolerance
      if (estimating.hessian) {
        int_tol <- integrate_tol_hessian
      } else {
        int_tol <- integrate_tol
      }

      int.vals <- c()
      for (ii in 1: length(xtilde.i)) {
        int.ii <- hcubature(
          f = lf,
          tol = integrate_tol,
          vectorInterface = TRUE,
          lowerLimit = 0, upperLimit = 1,
          xtilde = xtilde.i[ii],
          shape = alphas[ii],
          scale = betas[ii],
          sigsq_m = f.sigsq_m
        )
        int.vals[ii] <- int.ii$integral

        # If integral 0, set skip.rest to TRUE to skip further LL calculations
        if (int.ii$integral == 0) {
          print(paste("Integral is 0 for ii = ", ii, sep = ""))
          print(f.theta)
          skip.rest <- TRUE
          break
        }

      }
      ll.i <- sum(log(int.vals))

    } else {
      ll.i <- 0
    }

    # Return negative log-likelihood
    ll <- ll.p + ll.r + ll.i
    return(-ll)

  }

  # Create list of extra arguments, and assign default starting values and
  # lower values if not specified by user
  extra.args <- list(...)
  if (is.null(extra.args$start)) {
    if (! merror) {
      extra.args$start <- c(rep(0.01, n.gammas), rep(1, 2))
    } else {
      extra.args$start <- c(rep(0.01, n.gammas), rep(1, 3))
    }
  }
  if (is.null(extra.args$lower)) {
    if (! merror) {
      extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 2))
    } else {
      extra.args$lower <- c(rep(-Inf, n.gammas), rep(1e-3, 3))
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
  b1.hat <- ml.estimates[loc.bs[1]]
  b0.hat <- ml.estimates[loc.bs[2]]
  logOR.hat <- 1 / b0.hat - 1 / b1.hat

  # Obtain variance estimates
  if (estimate_var) {

    # Estimate Hessian
    hessian.mat <- pracma::hessian(f = llf, estimating.hessian = TRUE,
                                   x0 = ml.estimates)
    theta.variance <- try(solve(hessian.mat), silent = TRUE)
    if (class(theta.variance) == "try-error" ||
        ! all(eigen(x = theta.variance, only.values = TRUE)$values > 0)) {

      # Repeatedly divide integrate_tol_hessian by 5 and re-try
      while (integrate_tol_hessian > 1e-15 & fix_posdef) {
        integrate_tol_hessian <- integrate_tol_hessian / 5
        hessian.mat <- pracma::hessian(f = llf, estimating.hessian = TRUE,
                                       x0 = ml.estimates)
        theta.variance <- try(solve(hessian.mat), silent = TRUE)
        if (class(theta.variance) != "try-error" &&
            all(eigen(x = theta.variance, only.values = TRUE)$values > 0)) {
          break
        }

      }
    }

    if (class(theta.variance) == "try-error" ||
        ! all(eigen(x = theta.variance, only.values = TRUE)$values > 0)) {

      message("Estimated Hessian matrix is singular, so variance-covariance matrix cannot be obtained.")
      theta.variance <- NULL
      logOR.var <- NA

    } else {

      fprime <- matrix(c(1 / b1.hat^2, -1 / b0.hat^2), nrow = 1)
      colnames(theta.variance) <- rownames(theta.variance) <- theta.labels
      logOR.var <- fprime %*% theta.variance[loc.bs, loc.bs] %*% t(fprime)

    }

  } else {

    theta.variance <- NULL
    logOR.var <- NA

  }

  # Create vector of estimates to return
  estimates <- c(ml.estimates, logOR.hat, logOR.var)
  names(estimates) <- c(theta.labels, "logOR.hat", "logOR.var")

  # Create list to return
  ret.list <- list(estimates = estimates,
                   theta.var = theta.variance,
                   nlminb.object = ml.max,
                   aic = 2 * (length(ml.estimates) + ml.max$objective))
  return(ret.list)

}
