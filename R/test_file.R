# generate_data_hetero_cov <- function (n, # sample size
#                                       a, # intercept
#                                       b,
#                                       p,
#                                       gamma, # deterministic coefficient
#                                       chisq_df = 2) {
#
#     # K = 2
#
#     group_id <- as.integer(runif(n) > p) + 1
#     b0 <- b[group_id]
#
#     x <- (rchisq(n, df = chisq_df) - chisq_df) / sqrt(2 * chisq_df)
#     z_1 <- x + rnorm(n)
#     z_2 <- z_1 + rnorm(n)
#     z <- cbind(z_1, z_2)
#     sigma_i <- 0.5 * (1 + rchisq(n, df = 1))
#     e <- rnorm(n, mean = 0, sd = sqrt(sigma_i))
#
#     y <- a + x * b0 + as.numeric(z %*% gamma) + e
#
#     return(cbind(y, x, z))
# }
#
# R <- 1000
# res_mat <- array(0, c(R, 3, 3))
# for(r in 1:R) {
#
#     d <- generate_data_hetero_cov(1000, 0.5, c(1, 2), 0.5, c(1,1))
#     x <- d[, 'x']
#     y <- d[, 'y']
#     z <- d[, -c(1, 2)]
#
#     res_mat[r, ,1] <- ccrm_est_K2(x, y, z)$theta_b
#     res_mat[r, ,2] <- ccrm_est_K2(x, y, z, second_step = TRUE, weight_update = FALSE)$theta_b
#     res_mat[r, ,3] <- ccrm_est_K2(x, y, z, second_step = FALSE)$theta_b
# }
#
# sqrt(colMeans((res_mat[, , 1] - matrix(c(0.5, 1, 2), R, 3, byrow = TRUE))^2))
# sqrt(colMeans((res_mat[, , 2] - matrix(c(0.5, 1, 2), R, 3, byrow = TRUE))^2))
# sqrt(colMeans((res_mat[, , 3] - matrix(c(0.5, 1, 2), R, 3, byrow = TRUE))^2))

generate_data_vanilla_3 <- function(n, a, b, p, chisq_df = 2) {
    # K = 3

    sample_ind <- runif(n)
    group_id <- rowSums(
        sapply(p, function(x){
            as.integer(sample_ind > x)
        })
    )
    b0 <- b[group_id]

    x <- (rchisq(n, df = chisq_df) - chisq_df) / sqrt(2 * chisq_df)
    sigma_i <- 0.5 * (1 + rchisq(n, df = 1))
    u <- rnorm(n, mean = 0, sd = sqrt(sigma_i))

    y <- a + x * b0 + u

    return(cbind(y, x))
}

set.seed(10000)
n <- 100000
a <- 0.5
b <- c(1, 2, 3)
p <- c(0, 0.25, 0.75)
d <- generate_data_vanilla_3(n, a, b, p)
x <- d[, 'x']
y <- d[, 'y']
prop <- c(0.25, 0.5, 0.25)
sapply(1:5, function(x){
    sum(prop * (b^x))
})
res <- moment_est_3(x, y)
res

# weight_mat <- moment_est_gmm_3(x, y, NULL)$weight_mat
# res_gmm[6:10]
# moment_est_result_parameter <- res_gmm
# b1 <- moment_est_result_parameter[6]
# b2 <- moment_est_result_parameter[7]
# b3 <- moment_est_result_parameter[8]
# b4 <- moment_est_result_parameter[9]
# b5 <- moment_est_result_parameter[10]
# #
# b1 <- 2
# b2 <- 4.5
# b3 <- 11
# b4 <- 28.5
# b5 <- 77

# init_est_b_3(x, y, NULL)[6:11]

# M <- matrix(
#     c(1, b1, b2, b1, b2, b3, b2, b3, b4),
#     3, 3, byrow = TRUE
# )
# #
# D <- diag(c(1, -1, 1))
# #
# b_star <- solve(M %*% D) %*% c(b3, b4, b5) * c(-1, 1, -1)
# zz <-  polyroot(c(b_star, 1))
# zz <- Re(zz)
# #
# #
# solve(matrix(
#     c(
#         zz[1], zz[2], zz[3], zz[1]^2, zz[2]^2, zz[3]^2, zz[1]^3, zz[2]^3, zz[3]^3
#     ), 3, 3, byrow = TRUE
# ), c(b1, b2, b3))
#

#
#
# generate_data_vanilla <- function (n, # sample size
#                                    a, # intercept
#                                    b,
#                                    gamma, # deterministic coefficient
#                                    p,
#                                    dgp_para = 2) {
#
#     chisq_df <- dgp_para
#
#     group_id <- as.integer(runif(n) > p) + 1
#     b0 <- b[group_id]
#
#     x <- (rchisq(n, df = chisq_df) - chisq_df) / sqrt(2 * chisq_df)
#     z_1 <- x + rnorm(n)
#     z_2 <- z_1 + rnorm(n)
#     z <- cbind(z_1, z_2)
#     sigma_i <- 0.5 * (1 + rchisq(n, df = 1))
#     u <- rnorm(n, mean = 0, sd = sqrt(sigma_i))
#
#     y <- a + x * b0 + as.numeric(z %*% gamma) + u
#
#     return(cbind(y, x, z))
# }


# #
# set.seed(110)
# s_max = 3
# n = 100000
# a = 0.25
# b = c(1, 2)
# p = 0.5
# gamma = c(1, 1)
# data = generate_data_hetero_cov(n, a, b, p, gamma)
# y = data[, 1]
# x = data[, 2]
# z = data[, -c(1, 2)]

# coef_hat_ols <- lsfit(cbind(x,z), y_true)$coef
# gamma_hat <- coef_hat_ols[-(1:2)]
# y <- y_true - as.numeric(z %*% as.matrix(gamma_hat))

#
# # theta_m <- c(0.25, 1, 0, 1.5, 2.5, 4.5)
# # theta_b <- c(0.25, 1, 0, 0.5, 1, 2)
# theta_init_m <- moment_est(x, y)
# print(theta_init_m)
# theta_init <- init_est_b(x, y, s_max)
# print(theta_init)
# theta_hat <-  ccrm_est_K2(x, y, NULL, NULL, s_max)$theta
# print(theta_hat)
#
#
# theta_init_m_hetero <- moment_est_direct(x, y)
# print(theta_init_m)
# theta_init_hetero <- init_est_hetero(x, y, s_max)
# print(theta_init_hetero)
# theta_hat_hetero <- ccrm_est_hetero(x, y, 4)
# print(theta_hat_hetero)

# compute_u3 <- function(x, y) {
#     m1 <- mean(x)
#     m2 <- mean(x^2)
#     m3 <- mean(x^3)
#     m4 <- mean(x^4)
#     my1 <- mean(y)
#     my2 <- mean(y^2)
#     my3 <- mean(y^3)
#     my1x1 <- mean(x * y)
#     my2x2 <- mean(x^2 * y^2)
#     my3x1 <- mean(x^1 * y^3)
#     b1 <- 1.5
#     b2 <- 2.5
#     b3 <- 4.5
#     sigma_2  <- 1
#     a <- 0.25

#     u2 <- my2 - a^2 - 2 * a * m1 * b1 - m2*b2
#     u3 <- my3   - 3 * m2 * a * b2 - 3 * m1 * (a^2 + sigma_2) * b1 - (a^3 + 3 * a * sigma_2) - m3 * b3
#     c(u2, u3)

# }

# We want to test whether the new moment conditions work.






# setwd("C:/Users/zhang/Dropbox/Research/Categorical_Coef_Model/Numerical")
# source("dgp.R")
#
#
# # Sigma_mat <- matrix(0.25, 3, 3)
# # diag(Sigma_mat) <- 1
# # zz <- mvtnorm::rmvnorm(1000, sigma = Sigma_mat)
# # z <- zz^2
# # z_chisq <- rchisq(1000, df = 1)
#
#
# n = 10000 # sample size
# a = 0.25 # intercept
# b = c(1,2,3) # (b_1, b_2, ..., b_K)
# prob = c(0.3,0.3) # {p_1, p_2, ..., p_K} \in [0, 1]^K
# gamma = c(1,2) # deterministic coefficient
#
# K = 3
# d <- generate_data_hetero_cov(n, a, b, prob, gamma)
# y <- d[,1]
# x <- d[,2]
# z <- d[,-c(1,2)]
#
# res_init <- init_moment_est(x, y, z, K)
# theta_init <- c(res_init$gamma, res_init$Eb_hat)
#
# # -------------------------------
# # Eb_1_0 <- sum((1:3) * c(0.3,0.3,0.4)) #2.1
# # Eb_2_0 <- sum((1:3)^2 * c(0.3,0.3,0.4)) #5.1
# # Eb_3_0 <- sum((1:3)^3 * c(0.3,0.3,0.4)) #13.5
# # Eb_4_0 <- sum((1:3)^4 * c(0.3,0.3,0.4)) #37.5
# # Eb_5_0 <- sum((1:3)^5 * c(0.3,0.3,0.4)) #107.1
# # ------------------------------
#
# # d <- cbind(y,x,z)
#
# # -----------------------------------
# # # a loop to test init_est
# # eb = matrix(0, 500, 5)
# # for(r in 1:500){
# #
# #   d <- generate_data_hetero_cov(n, a, b, prob, gamma)
# #   y <- d[,1]
# #   x <- d[,2]
# #   z <- d[,-c(1,2)]
# #
# #   res_init <- init_moment_est(x, y, z, K)
# #   eb[r, ] = res_init$Eb_hat
# # }
#
#
# # ----------------------
# # start gmm
# t_0 <- Sys.time()
# res_gmm <- gmm::gmm(h_moment_mat_fn, d,
#                     t0 = theta_init,
#                     type = "twoStep")
# print(Sys.time() - t_0)
#
#
#
# # --------------------
#
# # ---------- try with new gmm
#
#
#
# multinom_gen <- function(r, exclude = TRUE) {
#
#   # Get all combinations of q_1, q_2, q_3
#   if (r == 2) {
#     multinom_order <- lapply(list(c(2,0,0), c(1,1,0)), function(l){allPerm(initMC(l))})
#   } else {
#     multinom_order <- lapply(genComp(r, 3, TRUE), function(l){allPerm(initMC(l))})
#   }
#   multinom_order <- do.call(rbind, multinom_order)
#
#   if (exclude) {
#     multinom_order <- multinom_order[(multinom_order[, 3] != 1), ]
#   }
#   # Compute multinomial coefficients
#   multinom_coef <- apply(multinom_order, 1, function(l){multinom(l, counts = TRUE)})
#
#
#   list(order = multinom_order,
#        coef = multinom_coef)
# }
