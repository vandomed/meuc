% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rc_cond_exp.R
\name{rc_cond_exp}
\alias{rc_cond_exp}
\title{Regression Calibration (Conditional Expectation Method)}
\usage{
rc_cond_exp(all_data = NULL, main = NULL, internal = NULL,
  external = NULL, y_var, z_var, d_vars = NULL, c_vars = NULL,
  b_vars = NULL, tdm_covariates = NULL, tdm_family = "gaussian",
  mem_covariates = NULL, mem_family = "gaussian",
  all_imputed = FALSE, boot_var = FALSE, boots = 100, alpha = 0.05)
}
\arguments{
\item{all_data}{Data frame with data for main study and validation study.}

\item{main}{Data frame with data for the main study.}

\item{internal}{Data frame with data for internal validation study.}

\item{external}{Data frame with data for the external validation study.}

\item{y_var}{Character string specifying name of Y variable.}

\item{z_var}{Character string specifying name of Z variable.}

\item{d_vars}{Character string specifying name of \strong{D} variables.}

\item{c_vars}{Character vector specifying names of \strong{C} variables.}

\item{b_vars}{Character vector specifying names of variables in true disease
model but not in measurement error model.}

\item{tdm_covariates}{Character vector specifying variables in true disease
model. The Z variable is automatically included whether you include it in
\code{tdm_covariates} or not.}

\item{tdm_family}{Character string specifying family of true disease model
(see \code{\link[stats]{glm}}).}

\item{mem_covariates}{Character vector specifying variables in measurement
error model.}

\item{mem_family}{Character string specifying family of measurement error
model (see \code{\link[stats]{glm}}).}

\item{all_imputed}{Logical value for whether to use imputed Z's for all
subjects, even those with the actual Z's observed (i.e. internal validation
subjects).}

\item{boot_var}{Logical value for whether to calculate a bootstrap
variance-covariance matrix.}

\item{boots}{Numeric value specifying number of bootstrap samples to use.}

\item{alpha}{Significance level for percentile bootstrap confidence interval.}
}
\value{
If \code{boot_var = TRUE}, list containing parameter estimates,
variance-covariance matrix, and percentile bootstrap confidence intervals;
otherwise just the parameter estimates.
}
\description{
Implements the "conditional expectation" version of regression calibration as
described by Rosner et al. (\emph{Stat. Med.} 1989). For the "algebraic"
version, see \code{\link{rc_algebraic}}.
}
\details{
The true disease model is a GLM:

g[E(Y)] = beta_0 + beta_z Z + \strong{beta_c}^T \strong{C} +
\strong{beta_b}^T \strong{B}

The measurement error model is:

h[E(Z)] = alpha_0 + \strong{alpha_d}^T \strong{D} + \strong{alpha_c}^T
\strong{C}

The procedure is as follows: in the validation study, fit the measurement
error model to estimate alpha's; in the main study, calculate
E(Z|\strong{D},\strong{C}) and fit the true disease model with those values
in place of the unobserved Z's.
}
\references{
Kuha, J. (1994) "Corrections for exposure measurement error in logistic
regression models with an application to nutritional data." \emph{Stat. Med.}
\strong{13}(11): 1135-1148.

Lyles, R.H. and Kupper, L.L. (2012) "Approximate and pseudo-likelihood
analysis for logistic regression using external validation data to model log
exposure." \emph{J. Agric. Biol. Environ. Stat.} \strong{18}(1): 22-38.

Rosner, B., Willett, W. and Spiegelman, D. (1989) "Correction of logistic
regression relative risk estimates and confidence intervals for systematic
within-person measurement error." \emph{Stat. Med.} \strong{8}(9): 1051-69.

Spiegelman, D., Carroll, R.J. and Kipnis, V. (2001) "Efficient regression
calibration for logistic regression in main study/internal validation study
designs with an imperfect reference instrument." \emph{Stat. Med.}
\strong{20}(1): 139-160.
}
