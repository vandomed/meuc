#' Normal Discriminant Function Approach for Estimating Odds Ratio with Exposure
#' Potentially Subject to Additive Normal Measurement Error
#'
#' Assumes exposure given covariates and outcome is a normal-errors linear
#' regression. Exposure measurements can be assumed precise or subject to
#' additive normal measurement error. Some replicates are required for
#' identifiability. Parameters are estimated using maximum likelihood.
#'
#'
#' @param y Numeric vector of Y values.
#' @param xtilde List of numeric vectors with Xtilde values.
#' @param c Numeric matrix with \strong{C} values (if any), with one row for
#' each subject. Can be a vector if there is only 1 covariate.
#' @param constant_or Logical value for whether to assume a constant odds ratio
#' for X, which means that sigsq_1 = sigsq_0. If \code{NULL}, model is fit with
#' and without this assumption, and a likelihood ratio test is performed to test
#' it.
#' @param merror Logical value for whether there is measurement error.
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
#' test for \code{H0: sigsq_1 = sigsq_0}, which is equivalent to
#' \code{H0: odds ratio is constant}.
#'
#'
#' @references
#' Lyles, R.H., Guo, Y., and Hill, A.N. (2009) "A fresh look at the discriminant
#' function approach for estimating crude or adjusted odds ratios."
#' \emph{Am. Stat} \strong{63}(4): 320--327.
#'
#'
#' @examples
#' # Load data frame with (Y, X, Xtilde, C) values for 1,000 subjects and list
#' # of Xtilde values where 25 subjects have replicates. Xtilde values are
#' # affected by measurement error. True log-OR = 0.5, sigsq = 1, sigsq_m = 0.5.
#' data(dat_ndfa)
#' dat <- dat_ndfa$dat
#' reps <- dat_ndfa$reps
#'
#' # Estimate log-OR for X and Y adjusted for C using true X values
#' # (unobservable truth).
#' fit.unobservable <- ndfa(
#'   y = dat$y,
#'   xtilde = dat$x,
#'   c = dat$c,
#'   merror = FALSE
#' )
#' fit.unobservable$estimates
#'
#' # Estimate log-OR for X and Y adjusted for C using observed Xtilde values,
#' # ignoring measurement error.
#' fit.naive <- ndfa(
#'   y = dat$y,
#'   xtilde = dat$xtilde,
#'   c = dat$c,
#'   merror = FALSE
#' )
#' fit.naive$estimates
#'
#' # Repeat, but accounting for measurement error.
#' fit.corrected <- ndfa(
#'   y = dat$y,
#'   xtilde = reps,
#'   c = dat$c,
#'   merror = TRUE
#' )
#' fit.corrected$estimates
#'
#' # Same as previous, but allowing for non-constant odds ratio.
#' fit.nonconstant <- ndfa(
#'   y = dat$y,
#'   xtilde = reps,
#'   c = dat$c,
#'   constant_or = FALSE,
#'   merror = TRUE
#' )
#' fit.nonconstant$estimates
#'
#' # Perform likelihood ratio test for H0: odds ratio is constant.
#' lrt <- ndfa(
#'   y = dat$y,
#'   xtilde = reps,
#'   c = dat$c,
#'   constant_or = NULL,
#'   merror = TRUE
#' )
#'
#'
#' @export
ndfa <- function(y,
                 xtilde,
                 c = NULL,
                 constant_or = TRUE,
                 merror = FALSE,
                 ...) {

  # Check that inputs are valid
  if (! is.null(constant_or) && ! is.logical(constant_or)) {
    stop("The input 'constant_or' should be set to TRUE, FALSE, or NULL.")
  }
  if (! is.logical(merror)) {
    stop("The input 'merror' should be TRUE or FALSE.")
  }

  if (is.null(constant_or)) {

    # Fit model with constant odds ratio
    fit.constant <- ndfa_constant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror
    )

    # Fit model with non-constant odds ratio
    fit.nonconstant <- ndfa_nonconstant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror
    )

    # Likelihood ratio test for H0: sigsq_1 = sigsq_0, which is equivalent to
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

    ret.list <- ndfa_constant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror
    )

  } else {

    ret.list <- ndfa_nonconstant(
      y = y,
      xtilde = xtilde,
      c = c,
      merror = merror
    )

  }

  return(ret.list)

}
