% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_func_K2.R
\name{init_est_b}
\alias{init_est_b}
\title{Initial estimation of the distributional parameters of \eqn{\beta_i}}
\usage{
init_est_b(
  x,
  y,
  z = NULL,
  s_max = 4,
  remove_intercept = TRUE,
  iter_gmm = TRUE,
  gmm_res = NULL,
  seed = 2023
)
}
\arguments{
\item{x}{regressor (N-by-1)}

\item{y}{dependent variable (N-by-p_x)}

\item{z}{control variables (N-by-p_z)}

\item{s_max}{Maximum order of moments used in estimation}

\item{remove_intercept}{whether subtract estimated intercept term in calculation of y_tilde}

\item{iter_gmm}{Whether to use iterative GMM (If FALSE, use OW-GMM with initial estimates)}

\item{gmm_res}{GMM estimation results from moment_est_gmm()}

\item{seed}{seed for random number generator to control the randomly generated result}
}
\value{
A list contains
\item{theta_hat}{estimated (pi, b_L, b_H) based on least squares optimization }
\item{sol_status}{whether estimated var(\beta_i) > 0 and we can solve for (pi, b_L, b_H)}
\item{moment_gmm_res}{A list contains the results of GMM moment estimation}
}
\description{
Initial estimation of the distributional parameters of \eqn{\beta_i}
}
