% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShiVa.R
\name{fit_OU_mean_var}
\alias{fit_OU_mean_var}
\title{Fit OU Model with Shifts in Mean and Variance}
\usage{
fit_OU_mean_var(
  tree,
  Y,
  alpha,
  shifts_mean,
  shifts_var,
  max.steps = 1000,
  t = 0.01,
  thres = 0.01,
  measurement_error = FALSE,
  max.num.shifts = Inf
)
}
\arguments{
\item{tree}{A phylogenetic tree of class \code{phylo}.}

\item{Y}{A numeric vector of continuous trait values at the tips of the tree.}

\item{alpha}{A non-negative numeric value specifying the strength of selection in the OU process.}

\item{shifts_mean}{An integer vector indicating the indices of branches where shifts in the optimal trait value occur.}

\item{shifts_var}{An integer vector indicating the indices of branches where shifts in evolutionary variance occur.}

\item{max.steps}{Maximum number of optimization steps. Default is 1000.}

\item{t}{Step size for the optimizer. Default is 0.01.}

\item{thres}{Convergence threshold for change in log-likelihood. Default is 0.01.}

\item{measurement_error}{Logical. If \code{TRUE}, a separate measurement error variance is estimated and added to the diagonal of the covariance matrix.}

\item{max.num.shifts}{Maximum allowed number of shifts (combined for mean and variance). Default is \code{Inf}.}
}
\value{
A list containing:
\item{tree}{The phylogenetic tree.}
\item{Y}{The trait values.}
\item{shifts_mean}{Indices of branches with non-zero shifts in optimal trait value.}
\item{shifts_var}{Indices of branches with non-zero shifts in evolutionary variance.}
\item{beta}{Estimated shift magnitudes for optima values.}
\item{gamma}{Estimated shift magnitudes for variance.}
\item{sigma2}{Estimated base evolutionary variance.}
\item{b0}{Estimated intercept (ancestral trait value).}
\item{sigma2_error}{Estimated measurement error variance (only returned if \code{measurement_error = TRUE}).}
\item{loglik}{Log-likelihood of the fitted model.}
\item{BIC}{BBIC for model selection.}
\item{mBIC}{mBIC for accounting shift sparsity and shared support.}
\item{pBIC}{pBIC incorporating determinant of projected design matrix.}
\item{fitted.values}{Fitted trait values based on the estimated model.}
\item{Sigma}{Estimated trait covariance matrix under the fitted model.}
}
\description{
Fits an Ornstein-Uhlenbeck (OU) model with user-specified shifts in both optimal trait values (mean) and evolutionary variance along a phylogeny. The method uses numerical optimization to estimate shift magnitudes, base variance, and intercept, and can optionally incorporate measurement error in trait values.
}
