#' Estimation: two categories
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-1)
#' @param theta_init initial value for distribution parameters
#' @param s_max
#' @param second_step
#'
#' @return A list contains estimated coefficients and inferential statistics
#' \item{theta_b}{Estimated distributional parameter p, b_L, b_H}
#' \item{theta_m}{Estimated moments of beta_i}
#' \item{theta_u}{Estimated moments of u_i}
#' \item{gamma}{estimated coefficients of control varibles, including intercept}
#' \item{V_theta}{variance}
#'
#' @export
#'
ccrm_est_K2 <- function(x, y, z, theta_init = NULL, s_max = 4, weight_mat = NULL, second_step = TRUE) {

    n <- length(x)

    s1 <- s_max - 1
    s2 <- s_max - 2
    s3 <- s_max - 3

    if (is.null(theta_init)) {
        theta_temp <- init_est_b(x, y, z, s_max)
        theta_init <- c(theta_temp$moment_u_hat, theta_temp$theta_hat)
        weight_mat <- theta_temp$weight_mat
    }

    # OLS estimate and replace gamma by gamma_hat
    if (!is.null(z)) {

        # OLS
        w <- cbind(1, x, z)
        phi_hat <- solve(t(w) %*% w, t(w) %*% y)
        # Residuals
        xi_hat <- c(y - w %*% phi_hat)
        # Variance and s.e.
        Q_ww <- (t(w) %*% w) / n
        D_hat <- Matrix::Diagonal(n)
        diag(D_hat) <- xi_hat^2
        V_wxi <- (t(w) %*% D_hat %*% w) / n
        V_hat <- (solve(Q_ww) %*% V_wxi %*% solve(Q_ww)) / n
        se_hat <- sqrt(Matrix::diag(V_hat))

        gamma_hat <- phi_hat[-(1:2)]
        y <- y - c(z %*% gamma_hat)

        gamma_hat <- phi_hat[-2]
        gamma_se_hat <- se_hat[-2]
    }

    # The following estimation is on y_tilde = a + x_i b_i + u_i
    # we don't replace a by a_hat because including the intercept
    #   stabilizes the finite sample performance

    # construct x matrix and xy matrix with different order
    x_mat <- sapply(0:s_max, function(i) {
        x^i
    })

    mxy1_mat <- y * x_mat
    mxy2_mat <- y^2 * x_mat
    mxy3_mat <- y^3 * x_mat

    mxy1 <- colMeans(y * x_mat)
    mxy2 <- colMeans(y^2 * x_mat)
    mxy3 <- colMeans(y^3 * x_mat)
    m <- colMeans(x_mat)

    # return moment matrix
    h_moment_mat_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3) p, b_L, b_H)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        p <- theta[4]
        b_L <- theta[5]
        b_H <- theta[6]

        Eb_1 <- p * b_L + (1 - p) * b_H
        Eb_2 <- p * b_L^2 + (1 - p) * b_H^2
        Eb_3 <- p * b_L^3 + (1 - p) * b_H^3

        h1 <-
            a * x_mat[, 1:(s1 + 1)] + Eb_1 * x_mat[, 2:(s1 + 2)] - mxy1_mat[, 1:(s1 + 1)]
        h2 <-
            (a^2 + sigma_2) * x_mat[, 1:(s2 + 1)] + 2 * a * Eb_1 * x_mat[, 2:(s2 + 2)] + Eb_2 * x_mat[, 3:(s2 + 3)] -
            mxy2_mat[, 1:(s2 + 1)]
        h3 <-
            (a^3 + 3 * a * sigma_2 + sigma_3) * x_mat[, 1:(s3 + 1)] + Eb_3 * x_mat[, 4:(s3 + 4)] +
            3 * x_mat[, 2:(s3 + 2)] * Eb_1 * (a^2 + sigma_2) + 3 * a * x_mat[, 3:(s3 + 3)] * Eb_2 - mxy3_mat[, 1:(s3 + 1)]

        cbind(h1, h2, h3)
    }

    h_moment_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3) p, b_L, b_H)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        p <- theta[4]
        b_L <- theta[5]
        b_H <- theta[6]

        b1 <- p * b_L + (1 - p) * b_H
        b2 <- p * b_L^2 + (1 - p) * b_H^2
        b3 <- p * b_L^3 + (1 - p) * b_H^3

        h1 <- m[2:(s1 + 2)] * b1 + m[1:(s1 + 1)] * a - mxy1[1:(s1 + 1)]
        h2 <- m[3:(s2 + 3)] * b2 + 2 * m[2:(s2 + 2)] * (a * b1) + m[1:(s2 + 1)] * (a^2 + sigma_2) - mxy2[1:(s2 + 1)]
        h3 <- m[4:(s3 + 4)] * b3 + 3 * m[3:(s3 + 3)] * b2 * a + 3 * m[2:(s3 + 2)] * b1 * (a^2 + sigma_2) +
                m[1:(s3 + 1)] * (a^3 + 3 * a * sigma_2 + sigma_3) - mxy3[1:(s3 + 1)]

        c(h1, h2, h3)

    }


    jac_h_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3) p, b_L, b_H)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        p <- theta[4]
        b_L <- theta[5]
        b_H <- theta[6]

        b1 <- p * b_L + (1 - p) * b_H
        b2 <- p * b_L^2 + (1 - p) * b_H^2
        b3 <- p * b_L^3 + (1 - p) * b_H^3

        jac_m <- rbind(
            cbind(m[1:(s1 + 1)], 0, 0, m[2:(s1 + 2)], 0, 0),
            cbind(
                2 * a * m[1:(s2 + 1)] + 2 * m[2:(s2 + 2)] * b1,
                m[1:(s2 + 1)],
                0,
                2 * a * m[2:(s2 + 2)],
                m[3:(s2 + 3)],
                0
            ),
            cbind(
                3 * a^2 * m[1:(s3 + 1)] + 6 * a * m[2:(s3 + 2)] * b1 +
                    3 * m[3:(s3 + 3)] * b2 + 3 * sigma_2 * m[1:(s3 + 1)],
                3 * m[2:(s3 + 2)] * b1 + 3 * a * m[1:(s3 + 1)],
                m[1:(s3 + 1)],
                3 * m[2:(s3 + 2)] * (a^2 + sigma_2),
                3 * a * m[3:(s3 + 3)],
                m[4:(s3 + 4)]
            )
        )

        jac_b <- matrix(c(1, 0, 0, 0, 0, 0,
                          0, 1, 0, 0, 0, 0,
                          0, 0, 1, 0, 0, 0,
                          0, 0, 0, b_L - b_H, p , 1-p,
                          0, 0, 0, b_L^2 - b_H^2, 2 * p * b_L , 2 * (1 - p) * b_H,
                          0, 0, 0, b_L^3 - b_H^3, 3 * p * b_L^2, 3 * (1 - p) * b_H^2),
                          6, 6, byrow = TRUE)

        jac_m %*% jac_b
    }

    # Estimate the weighting
    if (is.null(weight_mat)) {
        h_mat <- h_moment_mat_fn(theta_init)
        h_mean <- colMeans(h_mat)
        W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
    } else {
        W <- weight_mat
    }

    q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(t(h_fn_value) %*% W %*% h_fn_value)
    }

    grad_q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(2 * t(jac_h_fn(theta)) %*% W %*% h_fn_value)
    }

    # OPTIMIZATION
    eval_g_ineq <- function(theta) {
        list(
            "constraints" = c(theta[5] - theta[6]),
            "jacobian"= c(0, 0, 0, 0, 1, -1)
        )
    }

    # opts <- list(
    #     "algorithm" = "NLOPT_LD_LBFGS",
    #     "xtol_rel" = 1.0e-8
    # )
    local_opts <- list("algorithm" = "NLOPT_LD_MMA",
                       "xtol_rel"  = 1.0e-10)
    opts <- list(
        "algorithm" = "NLOPT_LD_AUGLAG",
        "xtol_rel" = 1.0e-15,
        "maxeval" = 1e6,
        "local_opts" = local_opts
    )

    nlopt_sol <- nloptr::nloptr(theta_init,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, 0, -Inf, 0, -Inf, -Inf),
                                ub = c(Inf, Inf, Inf, 1, Inf, Inf),
                                eval_g_ineq = eval_g_ineq,
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution

    # Two-step GMM
    if(second_step) {
        h_mat <- h_moment_mat_fn(theta_hat)
        h_mean <- colMeans(h_mat)
        W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
        nlopt_sol <- nloptr::nloptr(theta_hat,
                                    eval_f = q_fn,
                                    eval_grad_f = grad_q_fn,
                                    lb = c(-Inf, 0, -Inf, 0, -Inf, -Inf),
                                    ub = c(Inf, Inf, Inf, 1, Inf, Inf),
                                    eval_g_ineq = eval_g_ineq,
                                    opts = opts
        )
        theta_hat <- nlopt_sol$solution
    }

    # test statistics
    h_mat <- h_moment_mat_fn(theta_hat)
    h_mean <- colMeans(h_mat)
    W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
    D <- jac_h_fn(theta_hat)

    # Compute the moments
    p_hat <- theta_hat[4]
    b_L_hat <- theta_hat[5]
    b_H_hat <- theta_hat[6]
    Eb_1_hat <- p_hat * b_L_hat + (1 - p_hat) * b_H_hat
    Eb_2_hat <- p_hat * b_L_hat^2 + (1 - p_hat) * b_H_hat^2
    Eb_3_hat <- p_hat * b_L_hat^3 + (1 - p_hat) * b_H_hat^3

    # Compute the variance-covariance matrix of the estimated parameters
    # ERROR CATCHING
    # If pi ~ 0.5 and b_L ~ B_H, then D can be ill-posed (nearly rank deficient)
    if (!check_inverse(t(D) %*% W %*% D)) {
        V_theta <- NULL
        theta_se <- rep(NaN, length(theta_hat))
        kappa2_se <- NaN
        warning("Jacobian of moment function ill-posed.")
    } else {

        if(is.null(z)) {
            V_meat_mat <- t(h_mat) %*% (h_mat) / n
            sandwich_left <- solve(t(D) %*% W %*% D) %*% t(D) %*% W
            V_theta <- sandwich_left %*% V_meat_mat %*% t(sandwich_left) / n
        } else {
            L_mat <- diag(ncol(w))[-c(1, 2), ]
            G_gamma <- - rbind(t(x_mat[, 1:(s1 + 1)]) %*% z,
                             2 * t(mxy1_mat[, 1:(s2 + 1)]) %*% z,
                             3 * t(mxy2_mat[, 1:(s3 + 1)]) %*% z) / n
            meat_mat <- h_mat + t(G_gamma %*% L_mat %*% solve(Q_ww) %*% t(w * xi_hat))
            V_meat_mat <- t(meat_mat) %*% (meat_mat) / n
            sandwich_left <- solve(t(D) %*% W %*% D) %*% t(D) %*% W
            V_theta <- sandwich_left %*% V_meat_mat %*% t(sandwich_left) / n
        }

        H_grad_vec <- t(c(
            0, 0, 0,
            ((b_H_hat - b_L_hat)^2*(- b_H_hat^2*p_hat^2 + 2*b_H_hat^2*p_hat - b_H_hat^2 + b_L_hat^2*p_hat^2))/(b_L_hat^2*p_hat - b_H_hat^2*p_hat + b_H_hat^2)^2,
            -(2*b_H_hat*p_hat*(b_H_hat - b_L_hat)*(p_hat - 1)*(b_H_hat - b_H_hat*p_hat + b_L_hat*p_hat))/(b_L_hat^2*p_hat - b_H_hat^2*p_hat + b_H_hat^2)^2,
            (2*b_L_hat*p_hat*(b_H_hat - b_L_hat)*(p_hat - 1)*(b_H_hat - b_H_hat*p_hat + b_L_hat*p_hat))/(b_L_hat^2*p_hat - b_H_hat^2*p_hat + b_H_hat^2)^2
        ))
        kappa2_se <- c(sqrt(H_grad_vec %*% V_theta %*% t(H_grad_vec)))
    }


    if(is.null(z)) {
        gamma = NULL
        gamma_se = NULL
    }

    list(
        theta_b = theta_hat[4:6],
        theta_m = c(Eb_1_hat, p_hat * (1 - p_hat) * (b_L_hat - b_H_hat)^2, Eb_2_hat, Eb_3_hat),
        theta_m_gmm = theta_temp$moment_hat,
        theta_u = theta_hat[2:3],
        theta_se = sqrt(diag(V_theta)),
        gamma = gamma_hat,
        gamma_se = gamma_se_hat,
        V_theta = V_theta,
        kappa2 = Eb_1_hat^2 / Eb_2_hat,
        kappa2_se = kappa2_se
    )
}

# ----------First step estimation by solving equations----------
moment_est <- function(x, y) {
    n <- length(x)

    m1 <- mean(x)
    m2 <- mean(x^2)
    m3 <- mean(x^3)
    m4 <- mean(x^4)
    my1 <- mean(y)
    my2 <- mean(y^2)
    my3 <- mean(y^3)
    my1x1 <- mean(x * y)
    my2x2 <- mean(x^2 * y^2)
    my3x1 <- mean(x^1 * y^3)

    # First Moment
    a_mat <- matrix(c(
        1, m1,
        m1, m2
    ), 2, 2, byrow = TRUE)
    b_vec <- c(my1, my1x1)
    res_temp <- solve(a_mat, b_vec)
    a <- res_temp[1]
    b1 <- res_temp[2]

    # Second Moment
    a_mat <- matrix(c(
        1, m2,
        m2, m4
    ), 2, 2, byrow = TRUE)
    b_vec <- c(
        my2 - a^2 - 2 * a * m1 * b1,
        my2x2 - (a^2) * m2 - 2 * a * m3 * b1
    )
    res_temp <- solve(a_mat, b_vec)
    sigma_2 <- res_temp[1]
    b2 <- res_temp[2]


    # Third Moment
    a_mat <- matrix(c(
        1, m3,
        m1, m4
    ), 2, 2, byrow = TRUE)
    b_vec <- c(
        my3   - 3 * m2 * a * b2 - 3 * m1 * (a^2 + sigma_2) * b1 - (a^3 + 3 * a * sigma_2),
        my3x1 - 3 * m3 * a * b2 - 3 * m2 * (a^2 + sigma_2) * b1 - (a^3 + 3 * a * sigma_2) * m1
    )
    res_temp <- solve(a_mat, b_vec)
    sigma_3 <- res_temp[1]
    b3  <- res_temp[2]

    if (sigma_2 < 0) sigma_2 <- 0
    if (b2 < 0) b2 <- 0

    c(a, sigma_2, sigma_3, b1, b2, b3)
}



# ----------First step estimation GMM----------

#' Estimation of moments of beta_i
#'
#' @param x
#' @param y
#' @param z
#' @param s_max
#' @param second_step
#'
#' @export
#'
moment_est_gmm <- function(x, y, z, s_max, second_step = TRUE) {

    # with intercept, over-identification, GMM framework

    n <- length(x)

    s1 <- s_max - 1
    s2 <- s_max - 2
    s3 <- s_max - 3


    # OLS estimate and replace gamma by gamma_hat
    if (!is.null(z)) {

        # OLS
        w <- cbind(1, x, z)
        phi_hat <- solve(t(w) %*% w, t(w) %*% y)
        # Residuals
        xi_hat <- c(y - w %*% phi_hat)
        # Variance and s.e.
        Q_ww <- (t(w) %*% w) / n
        D_hat <- Matrix::Diagonal(n)
        diag(D_hat) <- xi_hat^2
        V_wxi <- (t(w) %*% D_hat %*% w) / n
        V_hat <- (solve(Q_ww) %*% V_wxi %*% solve(Q_ww)) / n
        se_hat <- sqrt(Matrix::diag(V_hat))

        gamma_hat <- phi_hat[-(1:2)]
        y <- y - c(z %*% gamma_hat)

        gamma_hat <- phi_hat[-2]
        gamma_se_hat <- se_hat[-2]
    }

    # construct x matrix and xy matrix with different order
    x_mat <- sapply(0:s_max, function(i) {
        x^i
    })

    mxy1_mat <- y * x_mat
    mxy2_mat <- y^2 * x_mat
    mxy3_mat <- y^3 * x_mat

    mxy1 <- colMeans(y * x_mat)
    mxy2 <- colMeans(y^2 * x_mat)
    mxy3 <- colMeans(y^3 * x_mat)
    m <- colMeans(x_mat)

    # return moment matrix
    h_moment_mat_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3) Eb_1, Eb_2, Eb_3)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        Eb_1 <- theta[4]
        Eb_2 <- theta[5]
        Eb_3 <- theta[6]

        h1 <-
            a * x_mat[, 1:(s1 + 1)] + Eb_1 * x_mat[, 2:(s1 + 2)] - mxy1_mat[, 1:(s1 + 1)]
        h2 <-
            (a^2 + sigma_2) * x_mat[, 1:(s2 + 1)] + 2 * a * Eb_1 * x_mat[, 2:(s2 + 2)] + Eb_2 * x_mat[, 3:(s2 + 3)] -
            mxy2_mat[, 1:(s2 + 1)]
        h3 <-
            (a^3 + 3 * a * sigma_2 + sigma_3) * x_mat[, 1:(s3 + 1)] + Eb_3 * x_mat[, 4:(s3 + 4)] +
            3 * x_mat[, 2:(s3 + 2)] * Eb_1 * (a^2 + sigma_2) + 3 * a * x_mat[, 3:(s3 + 3)] * Eb_2 - mxy3_mat[, 1:(s3 + 1)]

        cbind(h1, h2, h3)
    }

    h_moment_fn <- function(theta) {
        # theta = (alpha, sigma^2, E(u_i^3) Eb_1, Eb_2, Eb_3)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        b1 <- theta[4]
        b2 <- theta[5]
        b3 <- theta[6]

        h1 <- m[2:(s1 + 2)] * b1 + m[1:(s1 + 1)] * a - mxy1[1:(s1 + 1)]
        h2 <- m[3:(s2 + 3)] * b2 + 2 * m[2:(s2 + 2)] * (a * b1) + m[1:(s2 + 1)] * (a^2 + sigma_2) - mxy2[1:(s2 + 1)]
        h3 <- m[4:(s3 + 4)] * b3 + 3 * m[3:(s3 + 3)] * b2 * a + 3 * m[2:(s3 + 2)] * b1 * (a^2 + sigma_2) +
                m[1:(s3 + 1)] * (a^3 + 3 * a * sigma_2 + sigma_3) - mxy3[1:(s3 + 1)]

        c(h1, h2, h3)

    }


    jac_h_fn <- function(theta) {

         # theta = (alpha, sigma^2, E(u_i^3) b1, b2, b3)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        b1 <- theta[4]
        b2 <- theta[5]
        b3 <- theta[6]

        rbind(
            cbind(m[1:(s1 + 1)], 0, 0, m[2:(s1 + 2)], 0, 0),
            cbind(
                2 * a * m[1:(s2 + 1)] + 2 * m[2:(s2 + 2)] * b1,
                m[1:(s2 + 1)],
                0,
                2 * a * m[2:(s2 + 2)],
                m[3:(s2 + 3)],
                0
            ),
            cbind(
                3 * a^2 * m[1:(s3 + 1)] + 6 * a * m[2:(s3 + 2)] * b1 +
                    3 * m[3:(s3 + 3)] * b2 + 3 * sigma_2 * m[1:(s3 + 1)],
                3 * m[2:(s3 + 2)] * b1 + 3 * a * m[1:(s3 + 1)],
                m[1:(s3 + 1)],
                3 * m[2:(s3 + 2)] * (a^2 + sigma_2),
                3 * a * m[3:(s3 + 3)],
                m[4:(s3 + 4)]
            )
        )
    }



    # Estimate the weighting matrix
    theta_init <- moment_est(x, y)
    h_mat <- h_moment_mat_fn(theta_init)
    h_mean <- colMeans(h_mat)
    W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))

    q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(t(h_fn_value) %*% W %*% h_fn_value)
    }

    grad_q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(2 * t(jac_h_fn(theta)) %*% W %*% h_fn_value)
    }


    # OPTIMIZATION
    opts <- list(
        "algorithm" = "NLOPT_LD_LBFGS",
        "xtol_rel" = 1.0e-8
    )
    nlopt_sol <- nloptr::nloptr(theta_init,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, 0, -Inf, -Inf, 0, -Inf),
                                ub = c(Inf, Inf, Inf, Inf, Inf, Inf),
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution

    # Two-step GMM
    if(second_step) {
        h_mat <- h_moment_mat_fn(theta_hat)
        h_mean <- colMeans(h_mat)
        W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
        nlopt_sol <- nloptr::nloptr(theta_hat,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, 0, -Inf, -Inf, 0, -Inf),
                                ub = c(Inf, Inf, Inf, Inf, Inf, Inf),
                                opts = opts
        )
        theta_hat <- nlopt_sol$solution
    }

    # test statistics
    h_mat <- h_moment_mat_fn(theta_hat)
    h_mean <- colMeans(h_mat)
    W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
    D <- jac_h_fn(theta_hat)

    # Compute the variance-covariance matrix of the estimated parameters
    # ERROR CATCHING
    if (!check_inverse(t(D) %*% W %*% D)) {
        V_theta <- NULL
        warning("Jacobian of moment function ill-posed.")
    } else {

        if(is.null(z)) {
            V_meat_mat <- t(h_mat) %*% (h_mat) / n
            sandwich_left <- solve(t(D) %*% W %*% D) %*% t(D) %*% W
            V_theta <- sandwich_left %*% V_meat_mat %*% t(sandwich_left) / n
        } else {
            L_mat <- diag(ncol(w))[-c(1, 2), ]
            G_gamma <- rbind(t(x_mat[, 1:(s1 + 1)]) %*% z,
                             2 * t(mxy1_mat[, 1:(s2 + 1)]) %*% z,
                             3 * t(mxy2_mat[, 1:(s3 + 1)]) %*% z) / n
            meat_mat <- h_mat + t(G_gamma %*% L_mat %*% solve(Q_ww) %*% t(w * xi_hat))
            V_meat_mat <- t(meat_mat) %*% (meat_mat) / n
            sandwich_left <- solve(t(D) %*% W %*% D) %*% t(D) %*% W
            V_theta <- sandwich_left %*% V_meat_mat %*% t(sandwich_left) / n
        }

        kappa2 <- theta_hat[5] - theta_hat[4]^2
        H_grad_vec <- t(c(0, 0, 0, -2 * theta_hat[4], 1, 0))
        kappa2_se <- c(sqrt(H_grad_vec %*% V_theta %*% t(H_grad_vec)))
    }



    list(
        theta = theta_hat,
        theta_se = sqrt(diag(V_theta)),
        V_theta = V_theta,
        weight_mat = W,
        kappa2 = kappa2,
        kappa2_se = kappa2_se
    )
}

# ----------Second step estimation----------
init_est_b <- function(x,
                       y,
                       z,
                       s_max) {

    moment_est_result <- moment_est_gmm(x, y, z, s_max)
    moment_est_result_parameter <- moment_est_result$theta

    # first_step_est: estimated parameters only
    # alpha, sigma_2, Eb_1, Eb_2, b3
    a <- moment_est_result_parameter[1]
    sigma_2 <- moment_est_result_parameter[2]
    sigma_3 <- moment_est_result_parameter[3]
    b1 <- moment_est_result_parameter[4]
    b2 <- moment_est_result_parameter[5]
    b3 <- moment_est_result_parameter[6]


    f_fn <- function(theta) {
        f_value <- c(
            theta[1] * theta[2] + (1 - theta[1]) * theta[3],
            theta[1] * theta[2]^2 + (1 - theta[1]) * theta[3]^2,
            theta[1] * theta[2]^3 + (1 - theta[1]) * theta[3]^3
        ) - c(b1, b2, b3)
        f_value
    }
    jac_f_fn <- function(theta) {
        jac_mat <- matrix(c(
            theta[2] - theta[3], theta[1], 1 - theta[1],
            theta[2]^2 - theta[3]^2, 2 * theta[1] * theta[2], 2 * (1 - theta[1]) * theta[3],
            theta[2]^3 - theta[3]^3, 3 * theta[1] * theta[2]^2, 3 * (1 - theta[1]) * theta[3]^2
        ),
        3, 3,
        byrow = TRUE
        )
        jac_mat
    }
    g_fn <- function(theta) {
        # t(f) * f
        f_value <- f_fn(theta)
        sum(f_value^2)
    }
    grad_g_fn <- function(theta) {
        grad_g <- 2 * t(jac_f_fn(theta)) %*% f_fn(theta)
        as.numeric(grad_g)
    }

    gap <- initial_value_gap(b1, b2)
    theta_start <- c(0.5, b1 - gap, b1 + gap)

    eval_g_ineq <- function(theta) {
        list(
            "constraints" = c(theta[2] - theta[3]),
            "jacobian"= c(0, 1, -1)
        )

    }

    # opts <- list(
    #     "algorithm" = "NLOPT_LD_LBFGS",
    #     "xtol_rel" = 1.0e-8
    # )
    local_opts <- list("algorithm" = "NLOPT_LD_MMA",
                       "xtol_rel"  = 1.0e-10)
    opts <- list(
        "algorithm" = "NLOPT_LD_AUGLAG",
        "xtol_rel" = 1.0e-15,
        "maxeval" = 1e6,
        "local_opts" = local_opts
    )
    nlopt_sol <- nloptr::nloptr(theta_start,
                                eval_f = g_fn,
                                eval_grad_f = grad_g_fn,
                                lb = c(0, -Inf, -Inf),
                                ub = c(1, Inf, Inf),
                                eval_g_ineq = eval_g_ineq,
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution

    c(a, sigma_2, sigma_3, theta_hat, b1, b2, b3)
    list(
        theta_hat = theta_hat,
        moment_hat = c(b1, b2, b3), # From GMM
        moment_u_hat = c(a, sigma_2, sigma_3),
        weight_mat = moment_est_result$weight_mat
    )
}
