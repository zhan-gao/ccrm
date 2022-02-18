#' Estimation: heteroskedastic error, two categories
#'
#' @param x regressors (N-by-1)
#' @param y dependent variable (N-by-1)
#' @param s_max maximum order of moments of x
#' @param theta_init initial value for distribution parameters
#'
#' @return A list contains estimated coefficients and inferential statistics
#' \item{theta_hat}{Estimated coefficient: a, p, b_L, b_H, Eb1, Var_b}
#'
#' @export
#'

ccrm_est_hetero <- function(x,
                            y,
                            s_max,
                            theta_init = NULL
                            ) {

    # order of theta_init: a, p, b_L, b_H

    s1 <- s_max - 1
    s2 <- s_max - 2
    s3 <- s_max - 3

    if (is.null(theta_init)) {
        theta_temp <- init_est_hetero(x, y, s_max)
        theta_init <- theta_temp[1:4]
    }

    # construct x matrix and xy matrix with different order
    x_mat <- sapply(1:s_max, function(i) {
        x^i
    })
    m <- colMeans(x_mat)
    mxy1 <- colMeans(y * x_mat)
    mxy2 <- colMeans(y^2 * x_mat)
    mxy3 <- colMeans(y^3 * x_mat)

    my1 <- mean(y)
    my2 <- mean(y^2)
    my3 <- mean(y^3)

    h_moment_fn <- function(theta) {
        # theta = (a, p, b_L, b_H)
        # return a column vector stacking all moment conditions
        a <- theta[1]
        p <- theta[2]
        b_L <- theta[3]
        b_H <- theta[4]

        b1 <- p * b_L + (1 - p) * b_H
        b2 <- p * b_L^2 + (1 - p) * b_H^2
        b3 <- p * b_L^3 + (1 - p) * b_H^3

        h1 <- c(-my1 + a + m[1] * b1,
                -mxy1[1:s1] + m[1:s1] * a + m[2:s_max] * b1)
        h2 <- - mxy2[1:s2] + (
            m[1:s2] * my2 + b2 * (m[3:(s2 + 2)] - m[2] * m[1:s2]) + 2 * a * b1 * (m[2:(s2 + 1)] - m[1] * m[1:s2])
        )
        h3 <- - mxy3[1:s3] + (
            b1 * (3 * m[2:(s3+1)] * my2 - 3 * m[1] * m[1:s3] * my2) +
                a * (3 * my2 - 3 * m[1:s3] * my2) + m[1:s3] * my3 +
                a ^ 3 * (3 * m[1:s3] - 3) +
                b3 * (m[4:(s3 + 3)] - m[3] * m[1:s3]) -
                a * b2 * (3 * m[2] - 3 * m[3:(s3 + 2)]) -
                b1 * b2 * (3 * m[2] * m[2:(s3+1)] - 3 * m[1] * m[2] * m[1:s3]) -
                a * b1 ^ 2 * (-6 * m[1:s3] * m[1] ^ 2 + 6 * m[2:(s3+1)] * m[1]) -
                a ^ 2 * b1 * (6 * m[1] - 6 * m[1] * m[1:s3])
        )

        c(h1, h2, h3)
    }

    jac_h_fn <- function(theta) {

        # theta = (a, p, b_L, b_H)
        # return a Jacobian matrix qx4 where q is the number of moment conditions
        a <- theta[1]
        p <- theta[2]
        b_L <- theta[3]
        b_H <- theta[4]

        b1 <- p * b_L + (1 - p) * b_H
        b2 <- p * b_L^2 + (1 - p) * b_H^2
        b3 <- p * b_L^3 + (1 - p) * b_H^3

        m1 <- m[1]
        m2 <- m[2]
        m3 <- m[3]
        ms <- m[1:s3]
        ms1 <- m[2:(s3 + 1)]
        ms2 <- m[3:(s3 + 2)]
        ms3 <- m[4:(s3 + 3)]

        jac_m <- rbind(
            cbind(c(1, m[1:s1]), m, 0, 0),
            cbind(
                2 * b1 * (m[2:(s2 + 1)] - m[1] * m[1:s2]),
                2 * a  * (m[2:(s2 + 1)] - m[1] * m[1:s2]),
                m[3:(s2 + 2)] - m[2] * m[1:s2],
                0
            ),
            cbind(
                3*my2 - 3*b2*m2 + 3*b2*ms2 + 3*a^2*ms - ms*(3*my2 - (6*a + 6*b1*m1)*(a + b1*m1)) - (6*a + 6*b1*m1)*(a + b1*ms1) - 3*a^2 - 6*a*b1*m1 + 6*a*b1*ms1,
                3*a^2*ms1 + ms*(m1*(3*a^2 + 6*b1*m1*a - 3*my2 + 3*b2*m2) - 3*a^2*m1 + 6*a*m1*(a + b1*m1)) - ms1*(3*a^2 + 6*b1*m1*a - 3*my2 + 3*b2*m2) - 6*a*m1*(a + b1*ms1),
                3*a*ms2 - ms*(3*a*m2 - 3*m2*(a + b1*m1)) - 3*m2*(a + b1*ms1),
                ms3 - m3*ms
            )
        )

        jac_b <- matrix(c(1, 0, 0, 0,
                          0, b_L - b_H, p , 1-p,
                          0, b_L^2 - b_H^2, 2 * p * b_L , 2 * (1 - p) * b_H,
                          0, b_L^3 - b_H^3, 3 * p * b_L^2, 3 * (1 - p) * b_H^2),
                        4, 4, byrow = TRUE)

        jac_m %*% jac_b

    }


    # Estimate the weighting matrix
    # h_mat <- h_moment_mat_fn(theta_init)
    # h_mean <- colMeans(h_mat)
    # W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))

    q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(t(h_fn_value) %*% h_fn_value)
    }

    grad_q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(2 * t(jac_h_fn(theta)) %*% h_fn_value)
    }

    # OPTIMIZATION
    opts <- list(
        "algorithm" = "NLOPT_LD_LBFGS",
        "xtol_rel" = 1.0e-8
    )
    nlopt_sol <- nloptr::nloptr(theta_init,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, 0, -Inf, -Inf),
                                ub = c(Inf, 1, Inf, Inf),
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution

    p_hat <- theta_hat[2]
    b_L_hat <- theta_hat[3]
    b_H_hat <- theta_hat[4]
    Eb_1_hat <- p_hat * b_L_hat + (1 - p_hat) * b_H_hat
    Eb_2_hat <- p_hat * b_L_hat^2 + (1 - p_hat) * b_H_hat^2
    Eb_3_hat <- p_hat * b_L_hat^3 + (1 - p_hat) * b_H_hat^3

    # a, p, b_L, b_H, Eb1, Var_b
    c(theta_hat, Eb_1_hat, p_hat * (1 - p_hat) * (b_L_hat - b_H_hat)^2)
}

# ----------First step estimation by solving equations----------
moment_est_direct <- function(x, y) {

    # Passed a prelim test!

    m1 <- mean(x)
    m2 <- mean(x^2)
    m3 <- mean(x^3)
    m4 <- mean(x^4)
    my1 <- mean(y)
    my2 <- mean(y^2)
    my3 <- mean(y^3)
    mx2y2 <- mean(x^2 * y^2)
    mx1y3 <- mean(x * y^3)

    # OLS estimate of alpha and b1
    a_mat <- matrix(c(
        1, mean(x),
        mean(x), mean(x^2)
    ), 2, 2, byrow = TRUE)
    b_vec <- c(mean(y), mean(x * y))
    res_temp <- solve(a_mat, b_vec)
    a <- res_temp[1]
    b1 <- res_temp[2]

    # Estimate b2

    numer_2 <- (mx2y2 - m2 * my2) - 2 * (m3 - m1 * m2) * a * b1
    denom_2 <- m4 - m2^2
    b2 <- numer_2 / denom_2

    # Estimate b3
    denom_3 <- m4 - m1 * m3
    numer_3 <- mx1y3 - ( a*(3*my2 - 3*m1*my2) + m1*my3 + a^3*(3*m1 - 3) + b1*(- 3*my2*m1^2 + 3*m2*my2) - a*b2*(3*m2 - 3*m3) - a*b1^2*(- 6*m1^3 + 6*m2*m1) + b1*b2*(3*m1^2*m2 - 3*m2^2) - a^2*b1*(- 6*m1^2 + 6*m1))
    b3 <- numer_3 / denom_3

    if (b2 < 0) b2 <- 0
    c(a, b1, b2, b3)
}

# --------------------
#' First step estimation minimum distance
#'
#' @param x
#' @param y
#' @param s_max maximum order of moments of x
#'
#' @return theta_hat (a, b1, b2, b3)
#'
#'
moment_est_md <- function(x, y, s_max) {

    s1 <- s_max - 1
    s2 <- s_max - 2
    s3 <- s_max - 3

    # construct x matrix and xy matrix with different order
    x_mat <- sapply(1:s_max, function(i) {
        x^i
    })
    m <- colMeans(x_mat)
    mxy1 <- colMeans(y * x_mat)
    mxy2 <- colMeans(y^2 * x_mat)
    mxy3 <- colMeans(y^3 * x_mat)

    my1 <- mean(y)
    my2 <- mean(y^2)
    my3 <- mean(y^3)

    h_moment_fn <- function(theta) {
        # theta = (a, b1, b2, b3)
        # return a column vector stacking all moment conditions
        a <- theta[1]
        b1 <- theta[2]
        b2 <- theta[3]
        b3 <- theta[4]

        h1 <- c(-my1 + a + m[1] * b1,
                -mxy1[1:s1] + m[1:s1] * a + m[2:s_max] * b1)
        h2 <- - mxy2[1:s2] + (
            m[1:s2] * my2 + b2 * (m[3:(s2 + 2)] - m[2] * m[1:s2]) + 2 * a * b1 * (m[2:(s2 + 1)] - m[1] * m[1:s2])
        )
        h3 <- - mxy3[1:s3] + (
            b1 * (3 * m[2:(s3+1)] * my2 - 3 * m[1] * m[1:s3] * my2) +
            a * (3 * my2 - 3 * m[1:s3] * my2) + m[1:s3] * my3 +
            a ^ 3 * (3 * m[1:s3] - 3) +
            b3 * (m[4:(s3 + 3)] - m[3] * m[1:s3]) -
            a * b2 * (3 * m[2] - 3 * m[3:(s3 + 2)]) -
            b1 * b2 * (3 * m[2] * m[2:(s3+1)] - 3 * m[1] * m[2] * m[1:s3]) -
            a * b1 ^ 2 * (-6 * m[1:s3] * m[1] ^ 2 + 6 * m[2:(s3+1)] * m[1]) -
            a ^ 2 * b1 * (6 * m[1] - 6 * m[1] * m[1:s3])
        )

        c(h1, h2, h3)
    }

    jac_h_fn <- function(theta) {

        # theta = (alpha, b1, b2, b3)
        # return a Jacobian matrix qx4 where q is the number of moment conditions

        a <- theta[1]
        b1 <- theta[2]
        b2 <- theta[3]
        b3 <- theta[4]

        m1 <- m[1]
        m2 <- m[2]
        m3 <- m[3]
        ms <- m[1:s3]
        ms1 <- m[2:(s3 + 1)]
        ms2 <- m[3:(s3 + 2)]
        ms3 <- m[4:(s3 + 3)]

        rbind(
            cbind(c(1, m[1:s1]), m, 0, 0),
            cbind(
                2 * b1 * (m[2:(s2 + 1)] - m[1] * m[1:s2]),
                2 * a  * (m[2:(s2 + 1)] - m[1] * m[1:s2]),
                m[3:(s2 + 2)] - m[2] * m[1:s2],
                0
            ),
            cbind(
                3*my2 - 3*b2*m2 + 3*b2*ms2 + 3*a^2*ms - ms*(3*my2 - (6*a + 6*b1*m1)*(a + b1*m1)) - (6*a + 6*b1*m1)*(a + b1*ms1) - 3*a^2 - 6*a*b1*m1 + 6*a*b1*ms1,
                3*a^2*ms1 + ms*(m1*(3*a^2 + 6*b1*m1*a - 3*my2 + 3*b2*m2) - 3*a^2*m1 + 6*a*m1*(a + b1*m1)) - ms1*(3*a^2 + 6*b1*m1*a - 3*my2 + 3*b2*m2) - 6*a*m1*(a + b1*ms1),
                3*a*ms2 - ms*(3*a*m2 - 3*m2*(a + b1*m1)) - 3*m2*(a + b1*ms1),
                ms3 - m3*ms
            )
        )
    }

    # Estimate the weighting matrix
    # theta_init <- moment_est_hetero(x, y)
    # h_mat <- h_moment_mat_fn(theta_init)
    # h_mean <- colMeans(h_mat)
    # W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))

    theta_init <- moment_est_direct(x, y)

    q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(t(h_fn_value) %*% h_fn_value)
    }

    grad_q_fn <- function(theta) {
        h_fn_value <- h_moment_fn(theta)
        c(2 * t(jac_h_fn(theta)) %*% h_fn_value)
    }

    # OPTIMIZATION
    opts <- list(
        "algorithm" = "NLOPT_LD_LBFGS",
        "xtol_rel" = 1.0e-8
    )
    nlopt_sol <- nloptr::nloptr(theta_init,
                                eval_f = q_fn,
                                eval_grad_f = grad_q_fn,
                                lb = c(-Inf, -Inf, 0, -Inf),
                                ub = c(Inf, Inf, Inf, Inf),
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution

    return(theta_hat)
}

# ----------Second step estimation----------

init_est_hetero <- function(x,
                            y,
                            s_max,
                            method = "md") {

    if (method == "md") {
        moment_est_result_parameter <- moment_est_md(x, y, s_max)
    } else {
        moment_est_result_parameter <- moment_est_direct(x, y)
    }


    # first_step_est: estimated parameters only
    # a, b1, b2, b3
    a <- moment_est_result_parameter[1]
    b1 <- moment_est_result_parameter[2]
    b2 <- moment_est_result_parameter[3]
    b3 <- moment_est_result_parameter[4]
    var_b <- b2 - b1^2

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

    opts <- list(
        "algorithm" = "NLOPT_LD_LBFGS",
        "xtol_rel" = 1.0e-8
    )
    nlopt_sol <- nloptr::nloptr(theta_start,
                                eval_f = g_fn,
                                eval_grad_f = grad_g_fn,
                                lb = c(0, -Inf, -Inf),
                                ub = c(1, Inf, Inf),
                                opts = opts
    )
    theta_hat <- nlopt_sol$solution

    c(a, theta_hat, b1, b2, b3)
}
