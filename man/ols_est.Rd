% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ols_est.R
\name{ols_est}
\alias{ols_est}
\title{OLS Estimation}
\usage{
ols_est(x, y, z, remove_intercept = TRUE)
}
\arguments{
\item{x}{regressor (N-by-1)}

\item{y}{dependent variable (N-by-p_x)}

\item{z}{control variables (N-by-p_z)}

\item{remove_intercept}{whether subtract estimated intercept term in calculation of y_tilde}
}
\value{
A list contains estimated coefficients and inferential statistics
\item{phi}{Estimated phi}
\item{se}{Estimated  se for phi_hat}
\item{y_tilde}{\eqn{y - z \gamma}}
\item{xi_hat}{OLS residuals}
\item{remove_intecept}{Save the parameter remove_intercept for future reference}
}
\description{
OLS Estimation of \eqn{\phi} which consists of the first moment of \eqn{\beta_i}
and homogeneous slope coefficients \eqn{\gamma}
}
