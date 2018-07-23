#' Gamma Discriminant Function Approach for Estimating Odds Ratio with Exposure
#' Potentially Subject to Multiplicative Lognormal Measurement Error
#'
#' Assumes exposure given covariates and outcome is a constant-scale Gamma
#' regression. Exposure measurements can be assumed precise or subject to
#' multiplicative lognormal measurement error. Parameters are estimated using
#' maximum likelihood.
#'
#'
#' @param y Numeric vector of Y values.
#' @param xtilde Numeric vector (or list of numeric vectors, if there are
#' replicates) of Xtilde values.
#' @param c Numeric matrix with \strong{C} values (if any), with one row for
#' each subject. Can be a vector if there is only 1 covariate.
#' @param constant_or Logical value for whether to assume a constant odds ratio
#' for X, which means that gamma_y = 0. If \code{NULL}, model is fit with and
#' without this assumption, and likelihood ratio test is performed to test it.
#' @param merror Logical value for whether there is measurement error.
#' @param integrate_tol Numeric value specifying the \code{tol} input to
#' \code{\link{hcubature}}.
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
#' If \code{constant_or = NULL}, two such lists are returned (one under a
#' constant odds ratio assumption and one not), along with a likelihood ratio
#' test for \code{H0: gamma_y = 0}, which is equivalent to
#' \code{H0: odds ratio is constant}.
#'
#'
#' @references
#' Whitcomb, B.W., Perkins, N.J., Zhang, Z., Ye, A., and Lyles, R. H. (2012)
#' "Assessment of skewed exposure in case-control studies with pooling."
#' \emph{Stat. Med.} \strong{31}: 2461--2472.
#'
#'
#' @examples
#' # Load data frame with (Y, X, Xtilde, C) values for 250 subjects and list
#' # of Xtilde values where 25 subjects have replicates. Xtilde values are
#' # affected by measurement error. True log-OR = 0.5 and sigsq_m = 0.5.
#' data(dat_gdfa)
#' dat <- dat_gdfa$dat
#' reps <- dat_gdfa$reps
#'
#' # Estimate log-OR for X and Y adjusted for C using true X values
#' # (unobservable truth).
#' fit.unobservable <- gdfa(
#'   y = dat$y,
#'   xtilde = dat$x,
#'   c = dat$c,
#'   merror = FALSE
#' )
#' fit.unobservable$estimates
#'
#' # Estimate log-OR for X and Y adjusted for C using observed Xtilde values,
#' # ignoring measurement error.
#' fit.naive <- gdfa(
#'   y = dat$y,
#'   xtilde = dat$xtilde,
#'   c = dat$c,
#'   merror = FALSE
#' )
#' fit.naive$estimates
#'
#' # Repeat, but accounting for measurement error. Takes a few minutes to run
#' # due to numerical integration.
#' \dontrun{
#' fit.corrected <- gdfa(
#'   y = dat$y,
#'   xtilde = reps,
#'   c = dat$c,
#'   merror = TRUE,
#'   integrate_tol = 1e-4,
#'   control = list(trace = 1)
#' )
#' fit.corrected$estimates
#' }
#'
#' # Same as previous, but allowing for non-constant odds ratio.
#' \dontrun{
#' fit.nonconstant <- gdfa(
#'   y = dat$y,
#'   xtilde = reps,
#'   c = dat$c,
#'   constant_or = FALSE,
#'   merror = TRUE,
#'   integrate_tol = 1e-4,
#'   control = list(trace = 1)
#' )
#' fit.nonconstant$estimates
#' }
#'
#' # Perform likelihood ratio test for H0: odds ratio is constant.
#' \dontrun{
#' lrt <- gdfa(
#'   y = dat$y,
#'   xtilde = reps,
#'   c = dat$c,
#'   constant_or = NULL,
#'   merror = TRUE,
#'   integrate_tol = 1e-4,
#'   control = list(trace = 1)
#' )
#' lrt$fit.constant$estimates
#' }
#'
#'
#' @export
gdfa <- function(y,
                 xtilde,
                 c = NULL,
                 constant_or = TRUE,
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

  if (is.null(constant_or)) {

    # Fit model with constant odds ratio
    fit.constant <- gdfa_constant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror,
      integrate_tol = integrate_tol,
      integrate_tol_hessian = integrate_tol_hessian,
      estimate_var = estimate_var,
      fix_posdef = fix_posdef,
      ...
    )

    # Fit model with non-constant odds ratio
    fit.nonconstant <- gdfa_nonconstant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror,
      integrate_tol = integrate_tol,
      integrate_tol_hessian = integrate_tol_hessian,
      estimate_var = estimate_var,
      fix_posdef = fix_posdef,
      ...
    )

    # Likelihood ratio test for H0: gamma_y = 0, which is equivalent to
    # H0: constant odds ratio
    d <- 2 * (-fit.nonconstant$nlminb.object$objective + fit.constant$nlminb.object$objective)
    p <- pchisq(q = d, df = 1, lower.tail = FALSE)
    if (p < 0.05) {
      message <- "H0: Odds ratio is constant rejected at alpha = 0.05."
    } else {
      message <- "H0: Odds ratio is constant not rejected at alpha = 0.05."
    }
    message(message)
    lrt <- list(d = d, p = p, message = message)

    ret.list <- list(fit.constant = fit.constant,
                     fit.nonconstant = fit.nonconstant,
                     lrt = lrt)

  } else if (constant_or) {

    ret.list <- gdfa_constant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror,
      integrate_tol = integrate_tol,
      integrate_tol_hessian = integrate_tol_hessian,
      estimate_var = estimate_var,
      fix_posdef = fix_posdef,
      ...
    )

  } else {

    ret.list <- gdfa_nonconstant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror,
      integrate_tol = integrate_tol,
      integrate_tol_hessian = integrate_tol_hessian,
      estimate_var = estimate_var,
      fix_posdef = fix_posdef,
      ...
    )

  }

  return(ret.list)

}
