# -----------------------------------------------------------------------------
# Estimation of categorical coefficient regression models
#   with ONE random coefficient (p_x = 1) that has TWO categories (K = 2)
# -----------------------------------------------------------------------------

#' GMM Estimator of distributional parameters of a K = 2 random coefficient model
#'
#' This function estimates the distribution of beta_i using GMM. The moment conditions are
#' listed in Section 3 of Gao and Pesaran (2023)
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-p_x)
#' @param z control variables (N-by-p_z)
#' @param theta_init Initial estimates of distributional parameters
#' @param s_max Maximum order of moments used in estimation
#' @param remove_intercept whether subtract estimated intercept term in calculation of y_tilde
#' @param iter_gmm Whether to use iterative GMM (If FALSE, use OW-GMM with initial estimates)
#' @param seed seed for random number generator
#'
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
#'
ccrm_est_2 <- function(x, y, z = NULL, theta_init = NULL, s_max = 4, remove_intercept = TRUE, iter_gmm = TRUE, seed = 2023) {

    n <- length(x)

    # Check the support of x
    x_supp <- length(unique(x))
    if(x_supp == 1) {
        stop("x is homogeneous. ||Supp(x)||_0 >= 2 is required for identification.")
    }

    # Order of moments of x
    s1 <- min(x_supp - 1, s_max - 1)
    s2 <- min(x_supp - 1, s_max - 2)
    s3 <- min(x_supp - 1, s_max - 3)

    # Initial estimates
    gmm_res <- moment_est_gmm(x, y, z, s_max, remove_intercept, iter_gmm)
    if(gmm_res$kappa2 > 1) {
        cat(
            "Warning: the GMM estimate of the variance of beta is negative. \n",
            "THe test-statistic for homogeneity is ", gmm_res$kappa2_tstat, "\n",
            "The estimation will proceed, however, the results are not expected to be informative.\n"
        )
    }
    init_est_b_res <- init_est_b(x, y, z, s_max, remove_intercept, iter_gmm, gmm_res, seed)
    if(is.null(theta_init)) {
        theta_init <- c(gmm_res$theta[1:3], init_est_b_res$theta)
        weight_mat <- gmm_res$weight_mat
    }

    # OLS estimate and replace gamma by gamma_hat
    if (!is.null(z)) {
        y <- gmm_res$init_est$ols_res$y_tilde
        xi_hat <- gmm_res$init_est$ols_res$xi
        w <- cbind(1, x, z)
        Q_ww <- t(w) %*% w / n
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
    if(iter_gmm) {
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
        theta_b_se <- rep(NaN, 3)
        theta_m_se <- rep(NaN, 3)
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

        theta_b_se <- sqrt(diag(V_theta))[4:6]

        H_grad <- rbind(
            c(0, 0, 0, b_L_hat - b_H_hat, p_hat, 1 - p_hat),
            c(0, 0, 0, b_L_hat^2 - b_H_hat^2, 2 * p_hat * b_L_hat, 2 * (1 - p_hat) * b_H_hat),
            c(0, 0, 0, b_L_hat^3 - b_H_hat^3, 3 * p_hat * b_L_hat^2, 3 * (1 - p_hat) * b_H_hat^2)
        )
        theta_m_se <- sqrt(diag(H_grad %*% V_theta %*% t(H_grad)))
    }

    if(is.null(z)) {
        gamma_hat  <- NULL
        gamma_se_hat <- NULL
    } else {
        gamma_hat <- gmm_res$init_est$ols_res$phi[-2]
        gamma_se_hat <- gmm_res$init_est$ols_res$se[-2]
    }

    return(
        list(
            theta_b = theta_hat[4:6],
            theta_b_se = theta_b_se,
            theta_m = c(Eb_1_hat, Eb_2_hat, Eb_3_hat),
            theta_m_se = theta_m_se,
            theta_m_gmm = gmm_res$m_beta,
            theta_m_gmm_se = gmm_res$m_beta_se,
            gamma = gamma_hat,
            gamma_se = gamma_se_hat,
            V_theta = V_theta,
            gmm_res = gmm_res
        )
    )
}

# ------------------------------------------------------------------------------
# First step estimation by solving equations
# ------------------------------------------------------------------------------
#' Direct moment estimator using the identification equations
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-p_x)
#' @param z control variables (N-by-p_z)
#' @param remove_intercept whether subtract estimated intercept term in calculation of y_tilde
#'
#' @return A list contains
#' \item{para_est}{estimated (a, sigma2, sigma3, b1, b2, b3) }
#' \item{remove_intercept}{Save the parameter remove_intercept for future reference. = NA if z is not provided.}
#'
#'
moment_est_init <- function(x, y, z = NULL, remove_intercept = TRUE) {

    if (!is.null(z)) {
        ols_res <- ols_est(x, y, z, remove_intercept)
        y <- ols_res$y_tilde
    } else {
        ols_res <- NULL
    }

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
    if (!is.null(z)) {
        a <- ifelse(remove_intercept, 0, ols_res$phi[1])
        b1 <- ols_res$phi[2]
    } else {
        a_mat <- matrix(c(
            1, m1,
            m1, m2
        ), 2, 2, byrow = TRUE)
        b_vec <- c(my1, my1x1)
        res_temp <- solve(a_mat, b_vec)
        a <- res_temp[1]
        b1 <- res_temp[2]
    }

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

    list(
        para_est = c(a, sigma_2, sigma_3, b1, b2, b3),
        ols_res = ols_res
    )
}



# ------------------------------------------------------------------------------
# GMM Estimation of moments of beta_i
# ------------------------------------------------------------------------------

#' GMM Estimator of moments of beta_i, intercept and moments of u_i
#'
#' This function estimates the moments of beta_i using GMM. The moment conditions are
#' listed in Section 3 of Gao and Pesaran (2023)
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-p_x)
#' @param z control variables (N-by-p_z)
#' @param s_max Maximum order of moments used in estimation
#' @param remove_intercept whether subtract estimated intercept term in calculation of y_tilde
#' @param iter_gmm Whether to use iterative GMM (If FALSE, use OW-GMM with initial estimates)
#'
#' @return A list contains
#'  \item{m_beta}{estimated (b1, b2, b3)}
#'  \item{m_beta_se}{standard errors of m_beta}
#'  \item{theta}{estimated (a, sigma2, sigma3, b1, b2, b3)}
#'  \item{theta_se}{standard errors of theta}
#'  \item{V_theta}{estimated variance-covariance matrix of theta_hat}
#'  \item{weight_mat}{weight matrix used in estimation}
#'  \item{kappa2}{estimated kappa2 (b1^2 / b2)}
#'  \item{kappa2_se}{standard errors of kappa2}
#'  \item{kappa2_tstat}{t-statistic of kappa2}
#'  \item{init_est}{init_res}
#'
#' @import nloptr
#'
#' @export
#'
moment_est_gmm <- function(x, y, z = NULL, s_max = 4, remove_intercept = TRUE, iter_gmm = TRUE) {


    # sample size
    n <- length(x)

    # Check the support of x
    x_supp <- length(unique(x))
    if(x_supp == 1) {
        stop("x is homogeneous. ||Supp(x)||_0 >= 2 is required for identification.")
    }

    # Order of moments of x
    s1 <- min(x_supp - 1, s_max - 1)
    s2 <- min(x_supp - 1, s_max - 2)
    s3 <- min(x_supp - 1, s_max - 3)

    # OLS estimate and replace gamma by gamma_hat
    if (!is.null(z)) {
        init_res <- moment_est_init(x, y, z, remove_intercept)
        y <- init_res$ols_res$y_tilde
        xi_hat <- init_res$ols_res$xi
        w <- cbind(1, x, z)
        Q_ww <- t(w) %*% w / n
    } else {
        init_res <- moment_est_init(x, y)
    }
    theta_init <- init_res$para_est

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

    # Optimization
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
    if(iter_gmm) {
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

        kappa2 <- theta_hat[4]^2 / theta_hat[5]
        H_grad_vec <- t(c(
            0, 0, 0,
            2 * theta_hat[4] / theta_hat[5], - (theta_hat[4]^2) / (theta_hat[5]^2), 0
        ))
        kappa2_se <- c(sqrt(H_grad_vec %*% V_theta %*% t(H_grad_vec)))
        kappa2_tstat <- (kappa2 - 1) / kappa2_se
    }

    if(kappa2 > 1) {
        warning("Estimated variance of beta_i is negative...")
    }

    list(
        m_beta = theta_hat[4:6],
        m_beta_se = sqrt(diag(V_theta))[4:6],
        theta = theta_hat,
        theta_se = sqrt(diag(V_theta)),
        V_theta = V_theta,
        weight_mat = W,
        kappa2 = kappa2,
        kappa2_se = kappa2_se,
        kappa2_tstat = kappa2_tstat,
        init_est = init_res
    )
}


# ------------------------------------------------------------------------------
# Second step estimation
# ------------------------------------------------------------------------------

#' Initial estimation of the distributional parameters of \eqn{\beta_i}
#'
#' @param x regressor (N-by-1)
#' @param y dependent variable (N-by-p_x)
#' @param z control variables (N-by-p_z)
#' @param s_max Maximum order of moments used in estimation
#' @param remove_intercept whether subtract estimated intercept term in calculation of y_tilde
#' @param iter_gmm Whether to use iterative GMM (If FALSE, use OW-GMM with initial estimates)
#' @param gmm_res GMM estimation results from moment_est_gmm()
#' @param seed seed for random number generator to control the randomly generated result
#'
#' @return A list contains
#'  \item{theta_hat}{estimated (pi, b_L, b_H) based on least squares optimization }
#'  \item{sol_status}{whether estimated var(\beta_i) > 0 and we can solve for (pi, b_L, b_H)}
#'  \item{moment_gmm_res}{A list contains the results of GMM moment estimation}
#'
#' @import nloptr
#'
#' @export
#'
#'

init_est_b <- function(x, y, z = NULL, s_max = 4, remove_intercept = TRUE, iter_gmm = TRUE, gmm_res = NULL, seed = 2023) {

    if(!is.null(seed)) set.seed(seed)

    if(is.null(gmm_res)) {
        gmm_res <- moment_est_gmm(x, y, z, s_max, remove_intercept, iter_gmm)
    }
    theta_gmm <- gmm_res$theta

    # first_step_est: estimated parameters only
    # alpha, sigma_2, Eb_1, Eb_2, b3
    a <- theta_gmm[1]
    sigma_2 <- theta_gmm[2]
    sigma_3 <- theta_gmm[3]
    b1 <- theta_gmm[4]
    b2 <- theta_gmm[5]
    b3 <- theta_gmm[6]

    # Estimation
    if (b2 - b1^2 > 0) {
        # Direct estimation based on identification conditions
        b_plus <- (b3 - b1 * b2) / (b2 - b1^2)
        b_prod <- (b1 * b3 - b2^2) / (b2 - b1^2)
        b_L <- (b_plus - sqrt(b_plus^2 - 4 * b_prod)) / 2
        b_H <-  (b_plus + sqrt(b_plus^2 - 4 * b_prod)) / 2
        p <- (b_H - b1) / (b_H - b_L)
        theta_hat <- c(p, b_L, b_H)
        sol_status <-  1
    } else {
        warning("Estimated variance of beta_i is negative. A random solution is reported.")
        theta_hat <- c(
            runif(1),
            b1 - runif(1, 0.5, 1.5) * gap,
            b1 + runif(1, 0.5, 1.5) * gap
        )
        sol_status <-  0
    }

    list(
        theta_hat = theta_hat,
        sol_status = sol_status,
        moment_gmm_res = gmm_res
    )
}
