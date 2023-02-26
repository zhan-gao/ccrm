#' OLS Estimation
#'
#' OLS Estimation of \eqn{\phi} which consists of the first moment of \eqn{\beta_i}
#' and homogeneous slope coefficients \eqn{\gamma}
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-p_x)
#' @param z control variables (N-by-p_z)
#' @param remove_intercept whether subtract estimated intercept term in calculation of y_tilde
#'
#' @return A list contains estimated coefficients and inferential statistics
#' \item{phi}{Estimated phi}
#' \item{se}{Estimated  se for phi_hat}
#' \item{y_tilde}{\eqn{y - z \gamma}}
#' \item{xi_hat}{OLS residuals}
#' \item{remove_intecept}{Save the parameter remove_intercept for future reference}
#'
#' @import Matrix
#'
#' @export
#'

ols_est <- function(x, y, z, remove_intercept = TRUE) {

    n <- length(y)

    # OLS
    w <- cbind(1, x, z)
    phi_hat <- solve(t(w) %*% w, t(w) %*% y)
    # Residuals
    xi_hat <- y - w %*% phi_hat
    # Variance and s.e.
    Q_ww <- (t(w) %*% w) / n
    D_hat <- Matrix::Diagonal(n)
    diag(D_hat) <- xi_hat^2
    V_wxi <- (t(w) %*% D_hat %*% w) / n
    V_hat <- (solve(Q_ww) %*% V_wxi %*% solve(Q_ww)) / n
    se_hat <- sqrt(Matrix::diag(V_hat))

    if (remove_intercept) {
        gamma_hat <- phi_hat[-(2:(1 + ncol(as.matrix(x))))]
        y_tilde <- y - c(cbind(1, z) %*% gamma_hat)
    } else {
        gamma_hat <- phi_hat[-(1:(1 + ncol(as.matrix(x))))]
        y_tilde <- y - c(z %*% gamma_hat)
    }

    return(list(
        phi = phi_hat,
        se = se_hat,
        y_tilde = y_tilde,
        xi = xi_hat,
        remove_intercept = remove_intercept
    ))
}
