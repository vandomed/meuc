% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_logistic_logistic_linear.R
\name{ml_logistic_logistic_linear}
\alias{ml_logistic_logistic_linear}
\title{Maximum Likelihood with Three Models: Logistic Regression, Logistic
Regression, and Linear Regression}
\usage{
ml_logistic_logistic_linear(all_data = NULL, main = NULL,
  internal = NULL, external = NULL, y_var, x_var, z_var,
  d_vars = NULL, c_vars = NULL, b_vars = NULL,
  integrate_tol = 1e-08, integrate_tol_hessian = integrate_tol,
  estimate_var = TRUE, fix_posdef = FALSE, ...)
}
\arguments{
\item{all_data}{Data frame with data for main study and validation study.}

\item{main}{Data frame with data for the main study.}

\item{internal}{Data frame with data for internal validation study.}

\item{external}{Data frame with data for the external validation study.}

\item{y_var}{Character string specifying name of Y variable.}

\item{x_var}{Character string specifying name of X variable.}

\item{z_var}{Character string specifying name of Z variable.}

\item{d_vars}{Character string specifying name of \strong{D} variables.}

\item{c_vars}{Character vector specifying names of \strong{C} variables.}

\item{b_vars}{Character vector specifying names of variables in true disease
model but not in measurement error model.}

\item{integrate_tol}{Numeric value specifying \code{tol} input to
\code{\link[cubature]{hcubature}} for numerical integration.}

\item{integrate_tol_hessian}{Same as \code{integrate_tol}, but for use when
estimating the Hessian matrix only. Sometimes using a smaller value than for
likelihood maximization helps prevent cases where the inverse Hessian is not
positive definite.}

\item{estimate_var}{Logical value for whether to return variance-covariance
matrix for parameter estimates.}

\item{fix_posdef}{Logical value for whether to repeatedly reduce
\code{integrate_tol_hessian} by factor of 5 and re-estimate Hessian to try
to avoid non-positive definite variance-covariance matrix.}

\item{...}{Additional arguments to pass to \code{\link[stats]{nlminb}}.}
}
\value{
List containing:
\enumerate{
\item Numeric vector of parameter estimates.
\item Variance-covariance matrix (if \code{estimate_var = TRUE}).
\item Returned \code{\link[stats]{nlminb}} object from maximizing the
log-likelihood function.
\item Akaike information criterion (AIC).
}
}
\description{
Calculates maximum likelihood estimates for measurement error/unmeasured
confounding scenario where Y|(Z,X,\strong{C},\strong{B}) is logistic
regression, X|(Z,\strong{C},\strong{B}) is logistic regression, and
Z|(\strong{D},\strong{C}) is linear regression.
}
\details{
The true disease model is:

logit[P(Y = 1)] = beta_0 + beta_x X + beta_z Z + \strong{beta_c}^T \strong{C}
+ \strong{beta_b}^T \strong{B}

The X|(Z,\strong{C},\strong{B}) model is:

logit[P(X = 1)] = gamma_0 + gamma_z Z + \strong{gamma_c}^T \strong{C} +
\strong{gamma_b}^T \strong{B}

The Z|(\strong{D},\strong{C}) model is:

Z = alpha_0 + \strong{alpha_d}^T \strong{D} + \strong{alpha_c}^T \strong{C} +
d, d ~ N(0, sigsq_d)

There should be main study data with
(Y, X, \strong{D}, \strong{C}, \strong{B}) as well as internal validation
data with (Y, X, Z, \strong{D}, \strong{C}, \strong{B}) and/or external
validation data with (Z, X, \strong{D}, \strong{C}).
}
