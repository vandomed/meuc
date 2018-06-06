#' Fit Poolwise Regression Models
#'
#' Functions for implementing corrective methods for covariate measurement
#' error/missing data/unmeasured confounding. Methods include regression
#' calibration, propensity score calibration, and maximum likelihood based on
#' either two or three specified regression models.
#'
#'
#' \tabular{ll}{
#' Package: \tab meuc \cr
#' Type: \tab Package \cr
#' Version: \tab 1.1.1 \cr
#' Date: \tab 2018-06-06 \cr
#' License: \tab GPL-3 \cr
#' }
#'
#' @author Dane R. Van Domelen \cr \email{vandomed@@gmail.com}
#'
#'
#' @references
#' Acknowledgment: This material is based upon work supported by the National
#' Science Foundation Graduate Research Fellowship under Grant No. DGE-0940903.
#'
#'
#' @docType package
#' @importFrom cubature adaptIntegrate
#' @importFrom dplyr %>%
#' @import ggplot2
#' @importFrom Matrix bdiag
#' @importFrom mvtnorm dmvnorm
#' @importFrom numDeriv jacobian
#' @importFrom pracma hessian
#' @import stats
#' @name meuc
NULL
