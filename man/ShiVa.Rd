% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ShiVa.R
\name{ShiVa}
\alias{ShiVa}
\title{ShiVa: Automatic Shift Detection in Mean and Variance}
\usage{
ShiVa(
  Y,
  tree,
  alpha = NULL,
  t = 0.01,
  lambda1_list = NULL,
  lambda2_list = exp(1:10 * 0.4 - 6),
  criterion = "BIC",
  max.steps = 300,
  nfolds = 8,
  top_k = 10,
  measurement_error = FALSE,
  lambda.type = "lambda.1se",
  max.num.shifts = Inf
)
}
\arguments{
\item{Y}{A numeric vector of trait values at the tips of the phylogenetic tree.}

\item{tree}{A phylogenetic tree of class \code{phylo}.}

\item{alpha}{(Optional) A non-negative numeric value specifying the OU selection strength. If \code{NULL}, it is estimated via maximum likelihood using \code{phylolm()}.}

\item{t}{Step size for optimization. Default is 0.01.}

\item{lambda1_list}{A numeric vector of candidate \eqn{\lambda_1} values for penalizing mean shifts.}

\item{lambda2_list}{A numeric vector of candidate \eqn{\lambda_2} values for penalizing variance shifts. Default is \code{exp(1:10 * 0.4 - 6)}.}

\item{criterion}{Model selection criterion to use. Options are \code{"BIC"}, \code{"mBIC"}, or \code{"pBIC"}. Default is \code{"BIC"}.}

\item{max.steps}{Maximum number of optimization steps. Default is 300.}

\item{nfolds}{Number of folds for cross-validation in tuning \code{lambda1}. Default is 8.}

\item{top_k}{Number of top candidate models (based on criterion) to refine using backward correction. Default is 10.}

\item{measurement_error}{Logical. If \code{TRUE}, estimates a measurement error variance term. Default is \code{FALSE}.}

\item{lambda.type}{Cross-validation rule for selecting \code{lambda1}. Options are \code{"lambda.min"} or \code{"lambda.1se"}. Default is \code{"lambda.1se"}.}

\item{max.num.shifts}{Maximum number of allowed shifts (in both mean and variance). Default is \code{Inf}.}
}
\value{
A list with the same structure as \code{get_mean_var_shifts_model_selection}:
\item{best_model}{The final selected OU model object.}
\item{score_summary}{A data frame summarizing candidate models and backward-corrected scores.}
}
\description{
Performs automatic detection of evolutionary shifts in both optimal trait values (mean) and diffusion variance under an Ornstein-Uhlenbeck (OU) process. This function serves as a wrapper for \code{get_mean_var_shifts_model_selection}, with the added ability to automatically estimate the selection strength parameter \code{alpha} if not provided.
}
