# --------------------------------------------------------------------
# Bete version of the K = 3 functions. Not fully posished and tested.
# --------------------------------------------------------------------

#' Estimation: Three categories
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-1)
#' @param theta_init initial value for distribution parameters
#' @param s_max
#' @param weight_mat
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
ccrm_est_K3 <- function(x, y, z, theta_init = NULL, s_max = 6, weight_mat = NULL, second_step = TRUE) {

    n <- length(x)

    s1 <- s_max - 1
    s2 <- s_max - 2
    s3 <- s_max - 3
    s4 <- s_max - 4
    s5 <- s_max - 5

    if (is.null(theta_init)) {
        res_temp <- init_est_b_3(x, y, z, s_max)
        theta_init <- c(res_temp$moment_u_hat, res_temp$theta_hat)
        weight_mat <- res_temp$weight_mat
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
    mxy4_mat <- y^4 * x_mat
    mxy5_mat <- y^5 * x_mat


    mxy1 <- colMeans(y * x_mat)
    mxy2 <- colMeans(y^2 * x_mat)
    mxy3 <- colMeans(y^3 * x_mat)
    mxy4 <- colMeans(y^4 * x_mat)
    mxy5 <- colMeans(y^5 * x_mat)
    m <- colMeans(x_mat)

    # return moment matrix
    h_moment_mat_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3), E(u_i^4), E(u_i^5), p_L, p_M, p_H, b_L, b_M, b_H)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        sigma_4 <- theta[4]
        sigma_5 <- theta[5]
        p_L <- theta[6]
        p_M <- theta[7]
        p_H <- 1 - theta[6] - theta[7]
        b_L <- theta[8]
        b_M <- theta[9]
        b_H <- theta[10]

        Eb_1 <- p_L * b_L + p_M * b_M + p_H * b_H
        Eb_2 <- p_L * b_L^2 + p_M * b_M^2 + p_H * b_H^2
        Eb_3 <- p_L * b_L^3 + p_M * b_M^3 + p_H * b_H^3
        Eb_4 <- p_L * b_L^4 + p_M * b_M^4 + p_H * b_H^4
        Eb_5 <- p_L * b_L^5 + p_M * b_M^5 + p_H * b_H^5

        h1 <-
            a * x_mat[, 1:(s1 + 1)] + Eb_1 * x_mat[, 2:(s1 + 2)] - mxy1_mat[, 1:(s1 + 1)]
        h2 <-
            (a^2 + sigma_2) * x_mat[, 1:(s2 + 1)] + 2 * a * Eb_1 * x_mat[, 2:(s2 + 2)] + Eb_2 * x_mat[, 3:(s2 + 3)] -
            mxy2_mat[, 1:(s2 + 1)]
        h3 <-
            (a^3 + 3 * a * sigma_2 + sigma_3) * x_mat[, 1:(s3 + 1)] + Eb_3 * x_mat[, 4:(s3 + 4)] +
            3 * x_mat[, 2:(s3 + 2)] * Eb_1 * (a^2 + sigma_2) + 3 * a * x_mat[, 3:(s3 + 3)] * Eb_2 - mxy3_mat[, 1:(s3 + 1)]

        h4 <-
            (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) * x_mat[, 1:(s4 + 1)] +
            4 * x_mat[, 2:(s4 + 2)] * Eb_1 * (a^3 + 3 * a * sigma_2 + sigma_3)  +
            6 * x_mat[, 3:(s4 + 3)] * Eb_2 * (a^2 + sigma_2) +
            4 * a * x_mat[, 4:(s4 + 4)] * Eb_3 +
            Eb_4 * x_mat[, 5:(s4 + 5)] - mxy4_mat[, 1:(s4 + 1)]

        h5 <-
            (a^5 + 10 * a^3 * sigma_2 + 10 * a^2 * sigma_3 + 5 * a * sigma_4 + sigma_5) * x_mat[, 1:(s5 + 1)] +
            5 * x_mat[, 2:(s5 + 2)] * Eb_1 * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) +
            10 * x_mat[, 3:(s5 + 3)] * Eb_2 * (a^3 + 3 * a * sigma_2 + sigma_3) +
            10 * x_mat[, 4:(s5 + 4)] * Eb_3 * (a^2 + sigma_2) +
            5 * a * x_mat[, 5:(s5 + 5)] * Eb_4 +
            Eb_5 * x_mat[, 6:(s5 + 6)] - mxy5_mat[, 1:(s5 + 1)]

        cbind(h1, h2, h3, h4, h5)
    }

    h_moment_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3), E(u_i^4), E(u_i^5), p_L, p_M, p_H, b_L, b_M, b_H)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        sigma_4 <- theta[4]
        sigma_5 <- theta[5]
        p_L <- theta[6]
        p_M <- theta[7]
        p_H <- 1 - theta[6] - theta[7]
        b_L <- theta[8]
        b_M <- theta[9]
        b_H <- theta[10]

        b1 <- p_L * b_L + p_M * b_M + p_H * b_H
        b2 <- p_L * b_L^2 + p_M * b_M^2 + p_H * b_H^2
        b3 <- p_L * b_L^3 + p_M * b_M^3 + p_H * b_H^3
        b4 <- p_L * b_L^4 + p_M * b_M^4 + p_H * b_H^4
        b5 <- p_L * b_L^5 + p_M * b_M^5 + p_H * b_H^5

        h1 <- m[2:(s1 + 2)] * b1 + m[1:(s1 + 1)] * a - mxy1[1:(s1 + 1)]
        h2 <- m[3:(s2 + 3)] * b2 + 2 * m[2:(s2 + 2)] * (a * b1) + m[1:(s2 + 1)] * (a^2 + sigma_2) - mxy2[1:(s2 + 1)]
        h3 <- m[4:(s3 + 4)] * b3 + 3 * m[3:(s3 + 3)] * b2 * a + 3 * m[2:(s3 + 2)] * b1 * (a^2 + sigma_2) +
            m[1:(s3 + 1)] * (a^3 + 3 * a * sigma_2 + sigma_3) - mxy3[1:(s3 + 1)]
        h4 <- m[5:(s4 + 5)] * b4 +
            4 * m[4:(s4 + 4)] * b3 * a +
            6 * m[3:(s4 + 3)] * b2 * (a^2 + sigma_2) +
            4 * m[2:(s4 + 2)] * b1 * (a^3 + 3 * a * sigma_2 + sigma_3) +
            m[1:(s4 + 1)] * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) - mxy4[1:(s4 + 1)]
        h5 <- m[6:(s5 + 6)] * b5 +
            5 * m[5:(s5 + 5)] * b4 * a +
            10 * m[4:(s5 + 4)] * b3 * (a^2 + sigma_2) +
            10 * m[3:(s5 + 3)] * b2 * (a^3 + 3 * a * sigma_2 + sigma_3) +
            5 * m[2:(s5 + 2)] * b1 * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) +
            m[1:(s5 + 1)] * (a^5 + 10 * a^3 * sigma_2 + 10 * a^2 * sigma_3 + 5 * a * sigma_4 + sigma_5) - mxy5[1:(s5 + 1)]

        c(h1, h2, h3, h4, h5)

    }


    jac_h_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3), E(u_i^4), E(u_i^5), p_L, p_M, p_H, b_L, b_M, b_H)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        sigma_4 <- theta[4]
        sigma_5 <- theta[5]
        p_L <- theta[6]
        p_M <- theta[7]
        p_H <- 1 - theta[6] - theta[7]
        b_L <- theta[8]
        b_M <- theta[9]
        b_H <- theta[10]

        b1 <- p_L * b_L + p_M * b_M + p_H * b_H
        b2 <- p_L * b_L^2 + p_M * b_M^2 + p_H * b_H^2
        b3 <- p_L * b_L^3 + p_M * b_M^3 + p_H * b_H^3
        b4 <- p_L * b_L^4 + p_M * b_M^4 + p_H * b_H^4
        b5 <- p_L * b_L^5 + p_M * b_M^5 + p_H * b_H^5

        jac_m <- rbind(
            cbind(m[1:(s1 + 1)], 0, 0, 0, 0, m[2:(s1 + 2)], 0, 0, 0, 0),
            cbind(
                2 * a * m[1:(s2 + 1)] + 2 * m[2:(s2 + 2)] * b1,
                m[1:(s2 + 1)],
                0,
                0,
                0,
                2 * a * m[2:(s2 + 2)],
                m[3:(s2 + 3)],
                0,
                0,
                0
            ),
            cbind(
                3 * a^2 * m[1:(s3 + 1)] + 6 * a * m[2:(s3 + 2)] * b1 +
                    3 * m[3:(s3 + 3)] * b2 + 3 * sigma_2 * m[1:(s3 + 1)],
                3 * m[2:(s3 + 2)] * b1 + 3 * a * m[1:(s3 + 1)],
                m[1:(s3 + 1)],
                0,
                0,
                3 * m[2:(s3 + 2)] * (a^2 + sigma_2),
                3 * a * m[3:(s3 + 3)],
                m[4:(s3 + 4)],
                0,
                0
            ),
            cbind(
                (4 * a^3 + 12 * sigma_2 * a + 4 * sigma_3) * m[1:(s4 + 1)] +
                    12 * (a^2 + sigma_2)* m[2:(s4 + 2)] * b1 +
                    12 * a * m[3:(s4 + 3)] * b2 +
                    4 * m[4:(s4 + 4)] * b3,
                6 * a^2 * m[1:(s4 + 1)] +
                    12 * a * m[2:(s4 + 2)] * b1 +
                    6 * m[3:(s4 + 3)] * b2,
                4 * a * m[1:(s4 + 1)] + 4 * m[2:(s4 + 2)] * b1,
                m[1:(s4 + 1)],
                0,
                4 * m[2:(s4 + 2)] * (a^3 + 3 * a * sigma_2 + sigma_3),
                6 * m[3:(s4 + 3)] * (a^2 + sigma_2),
                4 * a * m[4:(s4 + 4)],
                m[5:(s4 + 5)],
                0
            ),
            cbind(
                (5 * a^4 + 30 * a^2 * sigma_2 + 20 * a * sigma_3 + 5 * sigma_4) * m[1:(s5 + 1)] +
                    5 * (4 * a^3 + 12 * a * sigma_2 + 4 * sigma_3) * m[2:(s5 + 2)] * b1 +
                    30 * (a^2 + sigma_2) * m[3:(s5 + 3)] * b2 +
                    20 * a * m[4:(s5 + 4)] * b3 +
                    5 * m[5:(s5 + 5)] * b4,
                10 * a^3 * m[1:(s5 + 1)] +
                    30 * a^2 * m[2:(s5 + 2)] * b1 +
                    30 * a * m[3:(s5 + 3)] * b2 +
                    10 * m[4:(s5 + 4)] * b3,
                10 * a^2 * m[1:(s5 + 1)] +
                    20 * a * m[2:(s5 + 2)] * b1 +
                    10 * m[3:(s5 + 3)] * b2,
                5 * a * m[1:(s5 + 1)] +
                    5 * m[2:(s5 + 2)] * b1,
                m[1:(s5 + 1)],
                5 * m[2:(s5 + 2)] * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4),
                10 * m[3:(s5 + 3)] * (a^ 3 + 3 * a * sigma_2 + sigma_3),
                10 * m[4:(s5 + 4)] * (a^2 + sigma_2),
                5 * a * m[5:(s5 + 5)],
                m[6:(s5 + 6)]
            )
        )

        jac_b <- rbind(
            cbind(
                diag(5), matrix(0, 5, 5)
            ),
            cbind(
                matrix(0, 5, 5),
                matrix(
                    c(
                        b_L - b_H, b_M - b_H, p_L, p_M, 1 - p_L - p_M,
                        b_L^2 - b_H^2, b_M^2 - b_H^2, 2 * p_L * b_L, 2 * p_M * b_M, 2 * (1 - p_L - p_M) * b_H,
                        b_L^3 - b_H^3, b_M^3 - b_H^3, 3 * p_L * b_L^2, 3 * p_M * b_M^2, 3 * (1 - p_L - p_M) * b_H^2,
                        b_L^4 - b_H^4, b_M^4 - b_H^4, 4 * p_L * b_L^3, 4 * p_M * b_M^3, 4 * (1 - p_L - p_M) * b_H^3,
                        b_L^5 - b_H^5, b_M^5 - b_H^5, 5 * p_L * b_L^4, 5 * p_M * b_M^4, 5 * (1 - p_L - p_M) * b_H^4
                    ), 5, 5, byrow = TRUE
                )
            )
        )
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

    eval_g_ineq <- function(theta){
        return(
            list(
                "constraints" = c(theta[6] + theta[7] - 1,
                                  theta[8] - theta[9],
                                  theta[9] - theta[10],
                                  theta[8] - theta[10]),
                "jacobian"= rbind(c(rep(0, 5), 1, 1, 0, 0, 0),
                                  c(rep(0, 5), 0, 0, 1, -1, 0),
                                  c(rep(0, 5), 0, 0, 0, 1, -1),
                                  c(rep(0, 5), 0, 0, 1, 0, -1))
            )
        )
    }

    # OPTIMIZATION
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
                                lb = c(-Inf, 0, -Inf, 0, -Inf, 0, 0, -Inf, -Inf, -Inf),
                                ub = c(Inf, Inf, Inf, Inf, Inf, 1, 1, Inf, Inf, Inf),
                                eval_g_ineq = eval_g_ineq,
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution

    # Two-step GMM
    if (second_step) {
        h_mat <- h_moment_mat_fn(theta_hat)
        h_mean <- colMeans(h_mat)
        if (!check_inverse((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))) {
            warning("Moment function ill-posed.")
            W  <- weight_mat
        } else {
            W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
        }
        nlopt_sol <- nloptr::nloptr(theta_hat,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, 0, -Inf, 0, -Inf, 0, 0, -Inf, -Inf, -Inf),
                                ub = c(Inf, Inf, Inf, Inf, Inf, 1, 1, Inf, Inf, Inf),
                                eval_g_ineq = eval_g_ineq,
                                opts = opts
        )
        theta_hat <- nlopt_sol$solution
    }

    # test statistics
    h_mat <- h_moment_mat_fn(theta_hat)
    h_mean <- colMeans(h_mat)
    if (!check_inverse((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))) {
        warning("Moment function ill-posed.")
        W  <- weight_mat
    } else {
        W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
    }
    D <- jac_h_fn(theta_hat)

    # Compute the moments
    p_L_hat <- theta_hat[6]
    p_M_hat <- theta_hat[7]
    p_H_hat <- 1 - theta_hat[6] - theta_hat[7]
    b_L_hat <- theta_hat[8]
    b_M_hat <- theta_hat[9]
    b_H_hat <- theta_hat[10]

    Eb_1_hat <- p_L_hat * b_L_hat + p_M_hat * b_M_hat + p_H_hat * b_H_hat
    Eb_2_hat <- p_L_hat * b_L_hat^2 + p_M_hat * b_M_hat^2 + p_H_hat * b_H_hat^2
    Eb_3_hat <- p_L_hat * b_L_hat^3 + p_M_hat * b_M_hat^3 + p_H_hat * b_H_hat^3
    Eb_4_hat <- p_L_hat * b_L_hat^4 + p_M_hat * b_M_hat^4 + p_H_hat * b_H_hat^4
    Eb_5_hat <- p_L_hat * b_L_hat^5 + p_M_hat * b_M_hat^5 + p_H_hat * b_H_hat^5

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
            G_gamma <- rbind(t(x_mat[, 1:(s1 + 1)]) %*% z,
                             2 * t(mxy1_mat[, 1:(s2 + 1)]) %*% z,
                             3 * t(mxy2_mat[, 1:(s3 + 1)]) %*% z,
                             4 * t(mxy3_mat[, 1:(s4 + 1)]) %*% z,
                             5 * t(mxy4_mat[, 1:(s5 + 1)]) %*% z) / n
            meat_mat <- h_mat + t(G_gamma %*% L_mat %*% solve(Q_ww) %*% t(w * xi_hat))
            V_meat_mat <- t(meat_mat) %*% (meat_mat) / n
            sandwich_left <- solve(t(D) %*% W %*% D) %*% t(D) %*% W
            V_theta <- sandwich_left %*% V_meat_mat %*% t(sandwich_left) / n
        }
    }


    if(is.null(z)) {
        gamma = NULL
        gamma_se = NULL
    }

    list(
        theta_b = theta_hat[6:10],
        theta_m = c(Eb_1_hat, Eb_2_hat - Eb_1_hat^2, Eb_2_hat, Eb_3_hat, Eb_4_hat, Eb_5_hat),
        theta_m_gmm = res_temp$moment_hat,
        theta_u = theta_hat[2:5],
        theta_se = sqrt(diag(V_theta)),
        gamma = gamma_hat,
        gamma_se = gamma_se_hat,
        V_theta = V_theta
    )
}

# ----------First step estimation by solving equations----------
moment_est_3 <- function(x, y) {

    m1 <- mean(x)
    m2 <- mean(x^2)
    m3 <- mean(x^3)
    m4 <- mean(x^4)
    m5 <- mean(x^5)
    m6 <- mean(x^6)
    my1 <- mean(y)
    my2 <- mean(y^2)
    my3 <- mean(y^3)
    my4 <- mean(y^4)
    my5 <- mean(y^5)
    my1x1 <- mean(x * y)
    my2x2 <- mean(x^2 * y^2)
    my3x1 <- mean(x^1 * y^3)
    my4x2 <- mean(x^2 * y^4)
    my5x1 <- mean(x^1 * y^5)

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

    # Fourth Moment
    a_mat <- matrix(c(
        1, m4,
        m2, m6
    ), 2, 2, byrow = TRUE)
    b_vec <- c(
        my4 - 4 * m3 * a * b3 - 6 * m2 * (a^2 + sigma_2) * b2 - 4 * m1 * (a^3 + 3 * a * sigma_2 + sigma_3) * b1 - (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3),
        my4x2 - 4 * m5 * a * b3 - 6 * m4 * (a^2 + sigma_2) * b2 - 4 * m3 * (a^3 + 3 * a * sigma_2 + sigma_3) * b1 - (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3) * m2
    )
    res_temp <- solve(a_mat, b_vec)
    sigma_4 <- res_temp[1]
    b4  <- res_temp[2]

    # Fifth Moment
    a_mat <- matrix(c(
        1, m5,
        m1, m6
    ), 2, 2, byrow = TRUE)
    b_vec <- c(
        my5 - 5 * m4 * a * b4 - 10 * m3 * (a^2 + sigma_2) * b3 - 10 * m2 * (a^3 + 3 * a * sigma_2 + sigma_3) * b2 - 5 * m1 * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) * b1 - (a^5 + 10 * a^3 * sigma_2 + 10 * a^2 * sigma_3 + 5 * a * sigma_4),
        my5x1 - 5 * m5 * a * b4 - 10 * m4 * (a^2 + sigma_2) * b3 - 10 * m3 * (a^3 + 3 * a * sigma_2 + sigma_3) * b2 - 5 * m2 * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) * b1 - (a^5 + 10 * a^3 * sigma_2 + 10 * a^2 * sigma_3 + 5 * a * sigma_4) * m1
    )
    res_temp <- solve(a_mat, b_vec)
    sigma_5 <- res_temp[1]
    b5  <- res_temp[2]

    if (sigma_2 < 0) sigma_2 <- 0
    if (b2 < 0) b2 <- 0
    if (sigma_4 < 0) sigma_4 <- 0
    if (b4 < 0) b4 <- 0

    c(a, sigma_2, sigma_3, sigma_4, sigma_5, b1, b2, b3, b4, b5)
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
moment_est_gmm_3 <- function(x, y, z, s_max = 6, second_step = TRUE) {

    # with intercept, over-identification, GMM framework

    n <- length(x)

    s1 <- s_max - 1
    s2 <- s_max - 2
    s3 <- s_max - 3
    s4 <- s_max - 4
    s5 <- s_max - 5


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
    mxy4_mat <- y^4 * x_mat
    mxy5_mat <- y^5 * x_mat


    mxy1 <- colMeans(y * x_mat)
    mxy2 <- colMeans(y^2 * x_mat)
    mxy3 <- colMeans(y^3 * x_mat)
    mxy4 <- colMeans(y^4 * x_mat)
    mxy5 <- colMeans(y^5 * x_mat)
    m <- colMeans(x_mat)

    # return moment matrix
    h_moment_mat_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3), E(u_i^4), E(u_i^5), Eb_1, Eb_2, Eb_3, Eb_4, Eb_5)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        sigma_4 <- theta[4]
        sigma_5 <- theta[5]

        Eb_1 <- theta[6]
        Eb_2 <- theta[7]
        Eb_3 <- theta[8]
        Eb_4 <- theta[9]
        Eb_5 <- theta[10]

        h1 <-
            a * x_mat[, 1:(s1 + 1)] + Eb_1 * x_mat[, 2:(s1 + 2)] - mxy1_mat[, 1:(s1 + 1)]
        h2 <-
            (a^2 + sigma_2) * x_mat[, 1:(s2 + 1)] + 2 * a * Eb_1 * x_mat[, 2:(s2 + 2)] + Eb_2 * x_mat[, 3:(s2 + 3)] -
            mxy2_mat[, 1:(s2 + 1)]
        h3 <-
            (a^3 + 3 * a * sigma_2 + sigma_3) * x_mat[, 1:(s3 + 1)] + Eb_3 * x_mat[, 4:(s3 + 4)] +
            3 * x_mat[, 2:(s3 + 2)] * Eb_1 * (a^2 + sigma_2) + 3 * a * x_mat[, 3:(s3 + 3)] * Eb_2 - mxy3_mat[, 1:(s3 + 1)]

        h4 <-
            (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) * x_mat[, 1:(s4 + 1)] +
            4 * x_mat[, 2:(s4 + 2)] * Eb_1 * (a^3 + 3 * a * sigma_2 + sigma_3)  +
            6 * x_mat[, 3:(s4 + 3)] * Eb_2 * (a^2 + sigma_2) +
            4 * a * x_mat[, 4:(s4 + 4)] * Eb_3 +
            Eb_4 * x_mat[, 5:(s4 + 5)] - mxy4_mat[, 1:(s4 + 1)]

        h5 <-
            (a^5 + 10 * a^3 * sigma_2 + 10 * a^2 * sigma_3 + 5 * a * sigma_4 + sigma_5) * x_mat[, 1:(s5 + 1)] +
            5 * x_mat[, 2:(s5 + 2)] * Eb_1 * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) +
            10 * x_mat[, 3:(s5 + 3)] * Eb_2 * (a^3 + 3 * a * sigma_2 + sigma_3) +
            10 * x_mat[, 4:(s5 + 4)] * Eb_3 * (a^2 + sigma_2) +
            5 * a * x_mat[, 5:(s5 + 5)] * Eb_4 +
            Eb_5 * x_mat[, 6:(s5 + 6)] - mxy5_mat[, 1:(s5 + 1)]

        cbind(h1, h2, h3, h4, h5)
    }

    h_moment_fn <- function(theta) {
        # theta = (alpha, sigma^2, E(u_i^3), E(u_i^4), E(u_i^5), Eb_1, Eb_2, Eb_3, Eb_4, Eb_5)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        sigma_4 <- theta[4]
        sigma_5 <- theta[5]

        b1 <- theta[6]
        b2 <- theta[7]
        b3 <- theta[8]
        b4 <- theta[9]
        b5 <- theta[10]

        h1 <- m[2:(s1 + 2)] * b1 + m[1:(s1 + 1)] * a - mxy1[1:(s1 + 1)]
        h2 <- m[3:(s2 + 3)] * b2 + 2 * m[2:(s2 + 2)] * (a * b1) + m[1:(s2 + 1)] * (a^2 + sigma_2) - mxy2[1:(s2 + 1)]
        h3 <- m[4:(s3 + 4)] * b3 + 3 * m[3:(s3 + 3)] * b2 * a + 3 * m[2:(s3 + 2)] * b1 * (a^2 + sigma_2) +
            m[1:(s3 + 1)] * (a^3 + 3 * a * sigma_2 + sigma_3) - mxy3[1:(s3 + 1)]
        h4 <- m[5:(s4 + 5)] * b4 +
            4 * m[4:(s4 + 4)] * b3 * a +
            6 * m[3:(s4 + 3)] * b2 * (a^2 + sigma_2) +
            4 * m[2:(s4 + 2)] * b1 * (a^3 + 3 * a * sigma_2 + sigma_3) +
            m[1:(s4 + 1)] * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) - mxy4[1:(s4 + 1)]
        h5 <- m[6:(s5 + 6)] * b5 +
            5 * m[5:(s5 + 5)] * b4 * a +
            10 * m[4:(s5 + 4)] * b3 * (a^2 + sigma_2) +
            10 * m[3:(s5 + 3)] * b2 * (a^3 + 3 * a * sigma_2 + sigma_3) +
            5 * m[2:(s5 + 2)] * b1 * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4) +
            m[1:(s5 + 1)] * (a^5 + 10 * a^3 * sigma_2 + 10 * a^2 * sigma_3 + 5 * a * sigma_4 + sigma_5) - mxy5[1:(s5 + 1)]

        c(h1, h2, h3, h4, h5)

    }


    jac_h_fn <- function(theta) {

        # theta = (alpha, sigma^2, E(u_i^3), E(u_i^4), E(u_i^5), Eb_1, Eb_2, Eb_3, Eb_4, Eb_5)
        a <- theta[1]
        sigma_2 <- theta[2]
        sigma_3 <- theta[3]
        sigma_4 <- theta[4]
        sigma_5 <- theta[5]

        b1 <- theta[6]
        b2 <- theta[7]
        b3 <- theta[8]
        b4 <- theta[9]
        b5 <- theta[10]

        rbind(
            cbind(m[1:(s1 + 1)], 0, 0, 0, 0, m[2:(s1 + 2)], 0, 0, 0, 0),
            cbind(
                2 * a * m[1:(s2 + 1)] + 2 * m[2:(s2 + 2)] * b1,
                m[1:(s2 + 1)],
                0,
                0,
                0,
                2 * a * m[2:(s2 + 2)],
                m[3:(s2 + 3)],
                0,
                0,
                0
            ),
            cbind(
                3 * a^2 * m[1:(s3 + 1)] + 6 * a * m[2:(s3 + 2)] * b1 +
                    3 * m[3:(s3 + 3)] * b2 + 3 * sigma_2 * m[1:(s3 + 1)],
                3 * m[2:(s3 + 2)] * b1 + 3 * a * m[1:(s3 + 1)],
                m[1:(s3 + 1)],
                0,
                0,
                3 * m[2:(s3 + 2)] * (a^2 + sigma_2),
                3 * a * m[3:(s3 + 3)],
                m[4:(s3 + 4)],
                0,
                0
            ),
            cbind(
                (4 * a^3 + 12 * sigma_2 * a + 4 * sigma_3) * m[1:(s4 + 1)] +
                    12 * (a^2 + sigma_2)* m[2:(s4 + 2)] * b1 +
                    12 * a * m[3:(s4 + 3)] * b2 +
                    4 * m[4:(s4 + 4)] * b3,
                6 * a^2 * m[1:(s4 + 1)] +
                    12 * a * m[2:(s4 + 2)] * b1 +
                    6 * m[3:(s4 + 3)] * b2,
                4 * a * m[1:(s4 + 1)] + 4 * m[2:(s4 + 2)] * b1,
                m[1:(s4 + 1)],
                0,
                4 * m[2:(s4 + 2)] * (a^3 + 3 * a * sigma_2 + sigma_3),
                6 * m[3:(s4 + 3)] * (a^2 + sigma_2),
                4 * a * m[4:(s4 + 4)],
                m[5:(s4 + 5)],
                0
            ),
            cbind(
                (5 * a^4 + 30 * a^2 * sigma_2 + 20 * a * sigma_3 + 5 * sigma_4) * m[1:(s5 + 1)] +
                    5 * (4 * a^3 + 12 * a * sigma_2 + 4 * sigma_3) * m[2:(s5 + 2)] * b1 +
                    30 * (a^2 + sigma_2) * m[3:(s5 + 3)] * b2 +
                    20 * a * m[4:(s5 + 4)] * b3 +
                    5 * m[5:(s5 + 5)] * b4,
                10 * a^3 * m[1:(s5 + 1)] +
                    30 * a^2 * m[2:(s5 + 2)] * b1 +
                    30 * a * m[3:(s5 + 3)] * b2 +
                    10 * m[4:(s5 + 4)] * b3,
                10 * a^2 * m[1:(s5 + 1)] +
                    20 * a * m[2:(s5 + 2)] * b1 +
                    10 * m[3:(s5 + 3)] * b2,
                5 * a * m[1:(s5 + 1)] +
                    5 * m[2:(s5 + 2)] * b1,
                m[1:(s5 + 1)],
                5 * m[2:(s5 + 2)] * (a^4 + 6 * a^2 * sigma_2 + 4 * a * sigma_3 + sigma_4),
                10 * m[3:(s5 + 3)] * (a^ 3 + 3 * a * sigma_2 + sigma_3),
                10 * m[4:(s5 + 4)] * (a^2 + sigma_2),
                5 * a * m[5:(s5 + 5)],
                m[6:(s5 + 6)]
            )
        )
    }

    # Estimate the weighting matrix
    theta_init <- moment_est_3(x, y)
    h_mat <- h_moment_mat_fn(theta_init)
    h_mean <- colMeans(h_mat)
    if (!check_inverse((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))) {
        warning("Moment function ill-posed.")
        W  <- diag(length(h_mean))
    } else {
        W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
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
    opts <- list(
        "algorithm" = "NLOPT_LD_LBFGS",
        "xtol_rel" = 1.0e-8
    )
    nlopt_sol <- nloptr::nloptr(theta_init,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, 0, -Inf, 0, -Inf, -Inf, 0, -Inf, 0, -Inf),
                                ub = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution


    # Two-step GMM
    if(second_step) {
        h_mat <- h_moment_mat_fn(theta_hat)
        h_mean <- colMeans(h_mat)
        if (!check_inverse((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))) {
            warning("Moment function ill-posed.")
            W  <- diag(length(h_mean))
        } else {
            W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
        }
        nlopt_sol <- nloptr::nloptr(theta_hat,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, 0, -Inf, 0, -Inf, -Inf, 0, -Inf, 0, -Inf),
                                ub = c(Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf, Inf),
                                opts = opts
        )
        theta_hat <- nlopt_sol$solution
    }

    # test statistics
    h_mat <- h_moment_mat_fn(theta_hat)
    h_mean <- colMeans(h_mat)
    if (!check_inverse((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))) {
        warning("Moment function ill-posed.")
        W  <- diag(length(h_mean))
    } else {
        W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
    }
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
                             3 * t(mxy2_mat[, 1:(s3 + 1)]) %*% z,
                             4 * t(mxy3_mat[, 1:(s4 + 1)]) %*% z,
                             5 * t(mxy4_mat[, 1:(s5 + 1)]) %*% z) / n
            meat_mat <- h_mat + t(G_gamma %*% L_mat %*% solve(Q_ww) %*% t(w * xi_hat))
            V_meat_mat <- t(meat_mat) %*% (meat_mat) / n
            sandwich_left <- solve(t(D) %*% W %*% D) %*% t(D) %*% W
            V_theta <- sandwich_left %*% V_meat_mat %*% t(sandwich_left) / n
        }

        # kappa2 <- theta_hat[5] - theta_hat[4]^2
        # H_grad_vec <- t(c(0, 0, 0, -2 * theta_hat[4], 1, 0))
        # kappa2_se <- c(sqrt(H_grad_vec %*% V_theta %*% t(H_grad_vec)))
    }

    list(
        theta = theta_hat,
        theta_se = sqrt(diag(V_theta)),
        weight_mat = W
    )
}



# ----------Second step estimation----------
init_est_b_3 <- function(x,
                         y,
                         z,
                         s_max = 6) {

    moment_est_result <- moment_est_gmm_3(x, y, z, s_max)
    moment_est_result_parameter <- moment_est_result$theta

    # first_step_est: estimated parameters only
    # theta = (alpha, sigma^2, E(u_i^3), E(u_i^4), E(u_i^5), Eb_1, Eb_2, Eb_3, Eb_4, Eb_5)

    a <- moment_est_result_parameter[1]
    sigma_2 <- moment_est_result_parameter[2]
    sigma_3 <- moment_est_result_parameter[3]
    sigma_4 <- moment_est_result_parameter[4]
    sigma_5 <- moment_est_result_parameter[5]
    b1 <- moment_est_result_parameter[6]
    b2 <- moment_est_result_parameter[7]
    b3 <- moment_est_result_parameter[8]
    b4 <- moment_est_result_parameter[9]
    b5 <- moment_est_result_parameter[10]

    f_fn <- function(theta) {

        p_L <- theta[1]
        p_M <- theta[2]
        p_H <- 1 - theta[1] - theta[2]
        b_L <- theta[3]
        b_M <- theta[4]
        b_H <- theta[5]

        f_value <- c(
            p_L * b_L + p_M * b_M + p_H * b_H ,
            p_L * b_L^2 + p_M * b_M^2 + p_H * b_H^2,
            p_L * b_L^3 + p_M * b_M^3 + p_H * b_H^3,
            p_L * b_L^4 + p_M * b_M^4 + p_H * b_H^4,
            p_L * b_L^5 + p_M * b_M^5 + p_H * b_H^5
        ) - c(b1, b2, b3, b4, b5)
        f_value
    }
    jac_f_fn <- function(theta) {
        p_L <- theta[1]
        p_M <- theta[2]
        p_H <- 1 - theta[1] - theta[2]
        b_L <- theta[3]
        b_M <- theta[4]
        b_H <- theta[5]

        jac_mat <- matrix(
            c(
                b_L - b_H, b_M - b_H, p_L, p_M, 1 - p_L - p_M,
                b_L^2 - b_H^2, b_M^2 - b_H^2, 2 * p_L * b_L, 2 * p_M * b_M, 2 * (1 - p_L - p_M) * b_H,
                b_L^3 - b_H^3, b_M^3 - b_H^3, 3 * p_L * b_L^2, 3 * p_M * b_M^2, 3 * (1 - p_L - p_M) * b_H^2,
                b_L^4 - b_H^4, b_M^4 - b_H^4, 4 * p_L * b_L^3, 4 * p_M * b_M^3, 4 * (1 - p_L - p_M) * b_H^3,
                b_L^5 - b_H^5, b_M^5 - b_H^5, 5 * p_L * b_L^4, 5 * p_M * b_M^4, 5 * (1 - p_L - p_M) * b_H^4
            ), 5, 5, byrow = TRUE
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
    theta_start <- c(0.25, 0.5, b1 - gap, b1, b1 + gap)

    eval_g_ineq <- function(theta){
        return(
            list(
                "constraints" = c(theta[1] + theta[2] - 1,
                                  theta[3] - theta[4],
                                  theta[4] - theta[5],
                                  theta[3] - theta[5]),
                "jacobian"= rbind(c(1, 1, 0, 0, 0),
                                  c(0, 0, 1, -1, 0),
                                  c(0, 0, 0, 1, -1),
                                  c(0, 0, 1, 0, -1))
            )
        )
    }

    local_opts <- list("algorithm" = "NLOPT_LD_MMA",
                       "xtol_rel"  = 1.0e-10)
    opts <- list(
        "algorithm" = "NLOPT_LD_AUGLAG",
        "xtol_rel" = 1.0e-15,
        "maxeval" = 1e6,
        "local_opts" = local_opts
    )

    nlopt_sol <- nloptr::nloptr(
        theta_start,
        eval_f = g_fn,
        eval_grad_f = grad_g_fn,
        lb = c(0, 0, -Inf, -Inf, -Inf),
        ub = c(1, 1,  Inf, Inf, Inf),
        eval_g_ineq = eval_g_ineq,
        opts = opts
    )
    theta_hat <- nlopt_sol$solution

    list(
        theta_hat = theta_hat,
        moment_hat = c(b1, b2, b3, b4, b5), # From GMM
        moment_u_hat = c(a, sigma_2, sigma_3, sigma_4, sigma_5),
        weight_mat = moment_est_result$weight_mat
    )
}
