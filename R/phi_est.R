phi_est <- function(x, y, z) {

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

    gamma_hat <- phi_hat[-(1:2)]
    y_tilde <- y - c(z %*% gamma_hat)

    return(list(
        phi = phi_hat,
        se = se_hat,
        y_tilde = y_tilde,
        xi = xi_hat
    ))
}
