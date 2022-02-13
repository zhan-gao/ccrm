generate_data_hetero_cov <- function (n, # sample size
                                      a, # intercept
                                      b,
                                      p,
                                      gamma, # deterministic coefficient
                                      chisq_df = 2) {

    # K = 2

    group_id <- as.integer(runif(n) > p) + 1
    b0 <- b[group_id]

    x <- (rchisq(n, df = chisq_df) - chisq_df) / sqrt(2 * chisq_df)
    z_1 <- x + rnorm(n)
    z_2 <- z_1 + rnorm(n)
    z <- cbind(z_1, z_2)
    sigma_i <- 0.5 * (1 + rchisq(n, df = 1))
    e <- rnorm(n, mean = 0, sd = sqrt(sigma_i))

    y <- a + x * b0 + as.numeric(z %*% gamma) + e

    return(cbind(y, x, z))
}

n = 10000
a = 0.25
b = c(1, 2)
p = 0.5
gamma = c(1, 1)
data = generate_data_hetero_cov(n, a, b, p, gamma)
y_true = data[, 1]
x = data[, 2]
z = data[, -c(1, 2)]

coef_hat_ols <- lsfit(cbind(x,z), y_true)$coef
gamma_hat <- coef_hat_ols[-(1:2)]
y <- y_true - as.numeric(z %*% as.matrix(gamma_hat))


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
