#' Estimation: homoskedastic error, two categories
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-1)
#' @param theta_init initial value for distribution parameters
#' @param s_1
#' @param s_2 
#' @param s_3
#'
#' @return A list contains estimated coefficients and inferential statistics
#' \item{theta_hat}{Estimated coefficient}
#' \item{wald_stat}{Wald statistic}
#' \item{t_stat}{t statistic}
#' \item{V_theta}{variance}
#' 
#' @export 
ccrm_est <- function(x,
                     y,
                     theta_init = NULL,
                     s_1 = 3,
                     s_2 = 2,
                     s_3 = 1) {

    # theta_init: estimated parameters only
    # alpha, sigma_2, pi, b_l, b_h


    n <- length(x)
    s_max <- max(c(s_1 + 1, s_2 + 2, s_3 + 3))

    if (is.null(theta_init)) {
        theta_temp <- init_est(x, y, s_1, s_2, s_3)$theta
        theta_init <- theta_temp[1:5]
    }


    # return moment matrix
    h_moment_mat_fn <- function(theta) {
        # theta = (alpha, sigma^2,  pi, b_L, b_H)
        alpha <- theta[1]
        sigma_2 <- theta[2]
        p <- theta[3]
        b_L <- theta[4]
        b_H <- theta[5]

        Eb_1 <- p * b_L + (1 - p) * b_H
        Eb_2 <- p * b_L^2 + (1 - p) * b_H^2
        Eb_3 <- p * b_L^3 + (1 - p) * b_H^3

        # dim: max{s_1 + 2, s_2 + 3, s_3 + 4}
        x_s_mat <- sapply(0:s_max, function(i) {
            x^i
        })
        Ex_s_mat <- x_s_mat
        Ey_1x_s_mat <- y * x_s_mat
        Ey_2x_s_mat <- y^2 * x_s_mat
        Ey_3x_s_mat <- y^3 * x_s_mat

        h_1 <-
            alpha * Ex_s_mat[, 1:(s_1 + 1)] + Eb_1 * Ex_s_mat[, 2:(s_1 + 2)] - Ey_1x_s_mat[, 1:(s_1 + 1)]
        h_2 <-
            (alpha^2 + sigma_2) * Ex_s_mat[, 1:(s_2 + 1)] + 2 * alpha * Eb_1 * Ex_s_mat[, 2:(s_2 + 2)] + Eb_2 * Ex_s_mat[, 3:(s_2 + 3)] -
            Ey_2x_s_mat[, 1:(s_2 + 1)]
        h_3 <-
            (alpha^3 + 3 * alpha * sigma_2) * Ex_s_mat[, 1:(s_3 + 1)] + Eb_3 * Ex_s_mat[, 4:(s_3 + 4)] +
            3 * Ex_s_mat[, 2:(s_3 + 2)] * Eb_1 * (alpha^2 + sigma_2) + 3 * alpha * Ex_s_mat[, 3:(s_3 + 3)] * Eb_2 - Ey_3x_s_mat[, 1:(s_3 + 1)]

        cbind(h_1, h_2, h_3)
    }

    jac_h_fn <- function(theta) {
        # gamma = (alpha, sigma^2, Eb_1, Eb_2, Eb_3)
        alpha <- theta[1]
        sigma_2 <- theta[2]
        p <- theta[3]
        b_L <- theta[4]
        b_H <- theta[5]

        Eb_1 <- p * b_L + (1 - p) * b_H
        Eb_2 <- p * b_L^2 + (1 - p) * b_H^2
        Eb_3 <- p * b_L^3 + (1 - p) * b_H^3

        # dim: max(s_1 + 2, s_2 + 3, s_3 + 4)
        x_s_mat <- sapply(0:s_max, function(i) {
            x^i
        })
        Ex_s <- colMeans(x_s_mat)

        rbind(
            cbind(Ex_s[1:(s_1 + 1)], 0, Ex_s[2:(s_1 + 2)] * (b_L - b_H), Ex_s[2:(s_1 + 2)] * p, Ex_s[2:(s_1 + 2)] * (1 - p)),
            cbind(
                2 * alpha * Ex_s[1:(s_2 + 1)] + 2 * Ex_s[2:(s_2 + 2)] * Eb_1,
                Ex_s[1:(s_2 + 1)],
                2 * alpha * Ex_s[2:(s_2 + 2)] * (b_L - b_H) + Ex_s[3:(s_2 + 3)] * (b_L^2 - b_H^2),
                2 * alpha * Ex_s[2:(s_2 + 2)] * p + 2 * Ex_s[3:(s_2 + 3)] * p * b_L,
                2 * alpha * Ex_s[2:(s_2 + 2)] * (1 - p) + 2 * Ex_s[3:(s_2 + 3)] * (1 - p) * b_H
            ),
            cbind(
                3 * alpha^2 * Ex_s[1:(s_3 + 1)] + 6 * alpha * Ex_s[2:(s_3 + 2)] * Eb_1 +
                    3 * Ex_s[3:(s_3 + 3)] * Eb_2 + 3 * sigma_2 * Ex_s[1:(s_3 + 1)],
                3 * Ex_s[2:(s_3 + 2)] * Eb_1 + 3 * alpha * Ex_s[1:(s_3 + 1)],
                3 * Ex_s[2:(s_3 + 2)] * (alpha^2 + sigma_2) * (b_L - b_H) + 3 * alpha * Ex_s[3:(s_3 + 3)] * (b_L^2 - b_H^2) +
                    Ex_s[4:(s_3 + 4)] * (b_L^3 - b_H^3),
                3 * Ex_s[2:(s_3 + 2)] * (alpha^2 + sigma_2) * p + 6 * alpha * Ex_s[3:(s_3 + 3)] * p * b_L +
                    3 * Ex_s[4:(s_3 + 4)] * p * b_L^2,
                3 * Ex_s[2:(s_3 + 2)] * (alpha^2 + sigma_2) * (1 - p) + 6 * alpha * Ex_s[3:(s_3 + 3)] * (1 - p) * b_H +
                    3 * Ex_s[4:(s_3 + 4)] * (1 - p) * b_H^2
            )
        )
    }

    # Estimate the weighting matrix
    h_mat <- h_moment_mat_fn(theta_init)
    h_mean <- colMeans(h_mat)
    W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))

    q_fn <- function(theta) {
        h_fn_value <- colMeans(h_moment_mat_fn(theta))
        c(t(h_fn_value) %*% W %*% h_fn_value)
    }

    grad_q_fn <- function(theta) {
        h_fn_value <- colMeans(h_moment_mat_fn(theta))
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
        lb = c(-Inf, 0, 0, -Inf, -Inf),
        ub = c(Inf, Inf, 1, Inf, Inf),
        opts = opts
    )
    theta_hat <- nlopt_sol$solution

    p_hat <- theta_hat[3]
    b_L_hat <- theta_hat[4]
    b_H_hat <- theta_hat[5]
    Eb_1_hat <- p_hat * b_L_hat + (1 - p_hat) * b_H_hat
    Eb_2_hat <- p_hat * b_L_hat^2 + (1 - p_hat) * b_H_hat^2
    Eb_3_hat <- p_hat * b_L_hat^3 + (1 - p_hat) * b_H_hat^3

    # test statistics
    h_mat <- h_moment_mat_fn(theta_hat)
    h_mean <- colMeans(h_mat)
    W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
    D <- jac_h_fn(theta_hat)
    r_fn <- b_L_hat - b_H_hat

    # ERROR CATCHING
    # If pi ~ 0.5 and b_L ~ B_H, then D can be ill-posed (nearly rank deficient)
    if (abs(Matrix::rcond(t(D) %*% W %*% D)) < 1e-16) {
        V_theta <- NULL
        warning("Jacobian of moment function ill-posed.")
        wald_stat <- 0
        t_stat <- c(NA, NA)
    } else {
        V_theta <- solve(t(D) %*% W %*% D)
        grad_r_fn <- c(0, 0, 0, 1, -1)
        wald_stat <- n * (r_fn^2) / c(t(grad_r_fn) %*% V_theta %*% grad_r_fn)
        t_stat <- theta_hat[c(4, 5)] / sqrt(diag(V_theta)[c(4, 5)] / n)
    }

    list(
        theta = c(theta_hat, Eb_1_hat, p_hat * (1 - p_hat) * (b_L_hat - b_H_hat)^2, Eb_2_hat, Eb_3_hat),
        wald_stat = wald_stat,
        t_stat = t_stat,
        V_theta = V_theta
    )
}
 
# ----------First step estimation by solving equations----------
moment_est <- function(x, y) {
    n <- length(x)

    a_mat <- matrix(c(
        1, mean(x),
        mean(x), mean(x^2)
    ), 2, 2, byrow = TRUE)
    b_vec <- c(mean(y), mean(x * y))
    res_temp <- solve(a_mat, b_vec)
    alpha <- res_temp[1]
    Eb_1 <- res_temp[2]


    a_mat <- matrix(c(
        1, mean(x^2),
        mean(x^2), mean(x^4)
    ), 2, 2, byrow = TRUE)
    b_vec <- c(
        mean(y^2) - alpha^2 - 2 * alpha * mean(x) * Eb_1,
        mean((x^2) * (y^2)) - (alpha^2) * mean(x^2) - 2 * alpha * mean(x^3) * Eb_1
    )
    res_temp <- solve(a_mat, b_vec)
    sigma_2 <- res_temp[1]
    Eb_2 <- res_temp[2]

    Eb_3 <-
        (mean(x * (y^3)) - alpha^3 * mean(x) - 3 * alpha^2 * mean(x^2) * Eb_1 -
            3 * alpha * mean(x^3) * Eb_2 - 3 * alpha * sigma_2 -
            3 * mean(x^2) * Eb_1 * sigma_2
        ) / mean(x^4)

    if (sigma_2 < 0) sigma_2 <- 0
    if (Eb_2 < 0) Eb_2 <- 0

    c(alpha, sigma_2, Eb_1, Eb_2, Eb_3)
}

# ----------First step estimation GMM----------
moment_est_gmm <- function(x, y, s_1 = 3, s_2 = 2, s_3 = 1) {

    # with intercept, over-identification, GMM framework

    n <- length(x)
    s_max <- max(c(s_1 + 1, s_2 + 2, s_3 + 3))

    # return moment matrix
    h_moment_mat_fn <- function(theta) {
        # gamma = (alpha, sigma^2, Eb_1, Eb_2, Eb_3)
        alpha <- theta[1]
        sigma_2 <- theta[2]
        Eb_1 <- theta[3]
        Eb_2 <- theta[4]
        Eb_3 <- theta[5]

        # dim: max(s_1 + 2, s_2 + 3, s_3 + 4)
        x_s_mat <- sapply(0:s_max, function(i) {
            x^i
        })
        Ex_s_mat <- x_s_mat
        Ey_1x_s_mat <- y * x_s_mat
        Ey_2x_s_mat <- y^2 * x_s_mat
        Ey_3x_s_mat <- y^3 * x_s_mat

        h_1 <-
            alpha * Ex_s_mat[, 1:(s_1 + 1)] + Eb_1 * Ex_s_mat[, 2:(s_1 + 2)] - Ey_1x_s_mat[, 1:(s_1 + 1)]
        h_2 <-
            (alpha^2 + sigma_2) * Ex_s_mat[, 1:(s_2 + 1)] + 2 * alpha * Eb_1 * Ex_s_mat[, 2:(s_2 + 2)] + Eb_2 * Ex_s_mat[, 3:(s_2 + 3)] -
            Ey_2x_s_mat[, 1:(s_2 + 1)]
        h_3 <-
            (alpha^3 + 3 * alpha * sigma_2) * Ex_s_mat[, 1:(s_3 + 1)] + Eb_3 * Ex_s_mat[, 4:(s_3 + 4)] +
            3 * Ex_s_mat[, 2:(s_3 + 2)] * Eb_1 * (alpha^2 + sigma_2) + 3 * alpha * Ex_s_mat[, 3:(s_3 + 3)] * Eb_2 - Ey_3x_s_mat[, 1:(s_3 + 1)]

        cbind(h_1, h_2, h_3)
    }

    jac_h_fn <- function(theta) {
        # gamma = (alpha, sigma^2, Eb_1, Eb_2, Eb_3)
        alpha <- theta[1]
        sigma_2 <- theta[2]
        Eb_1 <- theta[3]
        Eb_2 <- theta[4]
        Eb_3 <- theta[5]

        # dim: max(s_1 + 2, s_2 + 3, s_3 + 4)
        x_s_mat <- sapply(0:s_max, function(i) {
            x^i
        })
        Ex_s <- colMeans(x_s_mat)

        rbind(
            cbind(Ex_s[1:(s_1 + 1)], 0, Ex_s[2:(s_1 + 2)], 0, 0),
            cbind(
                2 * alpha * Ex_s[1:(s_2 + 1)] + 2 * Ex_s[2:(s_2 + 2)] * Eb_1,
                Ex_s[1:(s_2 + 1)],
                2 * alpha * Ex_s[2:(s_2 + 2)],
                Ex_s[3:(s_2 + 3)],
                0
            ),
            cbind(
                3 * alpha^2 * Ex_s[1:(s_3 + 1)] + 6 * alpha * Ex_s[2:(s_3 + 2)] * Eb_1 +
                    3 * Ex_s[3:(s_3 + 3)] * Eb_2 + 3 * sigma_2 * Ex_s[1:(s_3 + 1)],
                3 * Ex_s[2:(s_3 + 2)] * Eb_1 + 3 * alpha * Ex_s[1:(s_3 + 1)],
                3 * Ex_s[2:(s_3 + 2)] * (alpha^2 + sigma_2),
                3 * alpha * Ex_s[3:(s_3 + 3)],
                Ex_s[4:(s_3 + 4)]
            )
        )
    }

    # Estimate the weighting matrix
    theta_init <- moment_est(x, y)
    h_mat <- h_moment_mat_fn(theta_init)
    h_mean <- colMeans(h_mat)
    W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))

    q_fn <- function(theta) {
        h_fn_value <- colMeans(h_moment_mat_fn(theta))
        c(t(h_fn_value) %*% W %*% h_fn_value)
    }

    grad_q_fn <- function(theta) {
        h_fn_value <- colMeans(h_moment_mat_fn(theta))
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
        lb = c(-Inf, 0, -Inf, 0, -Inf),
        ub = c(Inf, Inf, Inf, Inf, Inf),
        opts = opts
    )
    theta_hat <- nlopt_sol$solution


    # test statistics
    h_mat <- h_moment_mat_fn(theta_hat)
    h_mean <- colMeans(h_mat)
    W <- solve((t(h_mat) %*% h_mat / n) - (h_mean %*% t(h_mean)))
    D <- jac_h_fn(theta_hat)
    V_theta <- solve(t(D) %*% W %*% D)

    r_fn <- theta_hat[4] - theta_hat[3]^2
    grad_r_fn <- c(0, 0, -2 * theta_hat[3], 1, 0)
    wald_stat <- n * (r_fn^2) / c(t(grad_r_fn) %*% V_theta %*% grad_r_fn)
    list(
        theta = theta_hat,
        wald_stat = wald_stat,
        V_theta = V_theta,
        var_b = theta_hat[4] - theta_hat[3]^2
    )
}

# ----------Second step estimation----------
initial_value_gap <- function(Eb_1, Eb_2) {

    # Using the estimated variance
    # If it's a valid estimate, return the s.d.
    # If not, return a random number between 0 and 1

    var_hat <- Eb_2 - Eb_1^2
    if (var_hat > 1e-4) {
        return(sqrt(var_hat))
    } else {
        return(runif(1))
    }
}

init_est <- function(x,
                     y,
                     s_1 = 3,
                     s_2 = 2,
                     s_3 = 1) {
    moment_est_result <- moment_est_gmm(x, y, s_1, s_2, s_3)
    moment_est_result_parameter <- moment_est_result$theta

    # first_step_est: estimated parameters only
    # alpha, sigma_2, Eb_1, Eb_2, Eb_3
    alpha <- moment_est_result_parameter[1]
    sigma_2 <- moment_est_result_parameter[2]
    Eb_1 <- moment_est_result_parameter[3]
    Eb_2 <- moment_est_result_parameter[4]
    Eb_3 <- moment_est_result_parameter[5]
    var_b <- Eb_2 - Eb_1^2


    f_fn <- function(theta) {
        f_value <- c(
            theta[1] * theta[2] + (1 - theta[1]) * theta[3],
            theta[1] * theta[2]^2 + (1 - theta[1]) * theta[3]^2,
            theta[1] * theta[2]^3 + (1 - theta[1]) * theta[3]^3
        ) - c(Eb_1, Eb_2, Eb_3)
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

    gap <- initial_value_gap(Eb_1, Eb_2)
    theta_start <- c(0.5, Eb_1 - gap, Eb_1 + gap)

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

    list(
        theta = c(alpha, sigma_2, theta_hat, Eb_1, var_b, Eb_2, Eb_3),
        wald_stat = moment_est_result$wald_stat
    )
}
